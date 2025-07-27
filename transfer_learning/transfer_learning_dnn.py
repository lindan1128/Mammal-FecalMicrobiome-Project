import os
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from tqdm import tqdm

# 1. Data loading
DATA_DIR = ''

profile_files = {
    'cow1': 'df_sub_core_cow1_corrected.csv',
    'cow2': 'df_sub_core_cow2_corrected.csv',
    'pig1': 'df_sub_core_pig1_corrected.csv',
    'pig2': 'df_sub_core_pig2_corrected.csv',
    'dog1': 'df_sub_core_dog1_corrected.csv',
    'edd_singh': 'df_sub_core_edd_singh_corrected.csv',
    'ibd_alm': 'df_sub_core_ibd_alm_corrected.csv',
    'ibd_engstrand': 'df_sub_core_ibd_engstrand_corrected.csv',
    'ibd_gevers': 'df_sub_core_ibd_gevers_corrected.csv',
    'ibd_huttenhower': 'df_sub_core_ibd_huttenhower_corrected.csv',
    'cdi_schubert': 'df_sub_core_cdi_schubert_corrected.csv',
    'cdi_vincent': 'df_sub_core_cdi_vincent_corrected.csv',
    'cdi_youngster': 'df_sub_core_cdi_youngster_corrected.csv'
}

label_files = {
    k: f'group_{k}.csv' for k in profile_files
}

profiles = {k: pd.read_csv(os.path.join(DATA_DIR, v), index_col=0).T for k, v in profile_files.items()}
labels = {k: pd.read_csv(os.path.join(DATA_DIR, label_files[k]))['label'].values for k in profile_files}

# 2. Host grouping
host_map = {
    'cow': ['cow1', 'cow2'],
    'pig': ['pig1', 'pig2'],
    'dog': ['dog1'],
    'human': ['edd_singh', 'ibd_alm', 'ibd_engstrand', 'ibd_gevers', 'ibd_huttenhower', 'cdi_schubert', 'cdi_vincent', 'cdi_youngster']
}

# 3. DNN model definitions
class SimpleDNN(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 80),
            nn.ReLU(),
            nn.Linear(80, 40),
            nn.ReLU(),
            nn.Linear(40, 20),
            nn.ReLU(),
            nn.Linear(20, 1)
        )
    def forward(self, x):
        return self.net(x).squeeze(-1)

class AttentionMLP(nn.Module):
    def __init__(self, input_dim, hidden_dim=48, dropout_rate=0.3):
        super().__init__()
        self.hidden_dim = hidden_dim
        
        # Feature projection layer
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        
        # Attention mechanism
        self.attention = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.Tanh(),
            nn.Linear(hidden_dim//2, 1)
        )
        
        # Feedforward network
        self.feed_forward = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.ReLU(),
            nn.Dropout(dropout_rate)
        )
        
        # Output layer
        self.output = nn.Linear(hidden_dim//2, 1)
        
    def forward(self, x):
        # Project to hidden dimension
        x = self.input_proj(x)  # [batch_size, hidden_dim]
        
        # Calculate attention weights
        attention_weights = torch.softmax(self.attention(x), dim=0)  # [batch_size, 1]
        
        # Apply attention
        attended_x = x * attention_weights  # [batch_size, hidden_dim]
        
        # Feedforward network
        ff_output = self.feed_forward(attended_x)
        
        # Output prediction
        return self.output(ff_output).squeeze(-1)

class AutoencoderClassifier(nn.Module):
    def __init__(self, input_dim, latent_dim=16, dropout_rate=0.3):
        super().__init__()
        # Encoder - smaller hidden layers
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 40),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(40, 28),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(28, latent_dim)
        )
        
        # Decoder - correspondingly reduced
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 28),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(28, 40),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(40, input_dim)
        )
        
        # Classifier - simplified
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, 8),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(8, 1)
        )
    
    def forward(self, x, return_reconstruction=False):
        encoded = self.encoder(x)
        
        if return_reconstruction:
            decoded = self.decoder(encoded)
            classification = self.classifier(encoded).squeeze(-1)
            return classification, decoded
        else:
            return self.classifier(encoded).squeeze(-1)

# Function to count model parameters
def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def print_model_sizes(input_dim=50):
    """Print parameter count comparison of three models"""
    print(f"\nModel Parameter Count Comparison (input_dim={input_dim}):")
    print("="*50)
    
    models = {
        'SimpleDNN': SimpleDNN(input_dim),
        'AttentionMLP': AttentionMLP(input_dim), 
        'AutoencoderClassifier': AutoencoderClassifier(input_dim)
    }
    
    for name, model in models.items():
        param_count = count_parameters(model)
        print(f"{name:20}: {param_count:,} parameters")
    
    print("="*50 + "\n")

# 4. Training and evaluation functions
def train_model(model, loader, criterion, optimizer, device, is_autoencoder=False):
    model.train()
    for x, y in loader:
        x, y = x.to(device), y.to(device)
        optimizer.zero_grad()
        
        if is_autoencoder and hasattr(model, 'decoder'):
            # AutoencoderClassifier training
            classification, reconstruction = model(x, return_reconstruction=True)
            class_loss = criterion(classification, y)
            recon_loss = nn.MSELoss()(reconstruction, x)
            loss = class_loss + 0.1 * recon_loss  # Balance classification and reconstruction loss
        else:
            # Standard model training
            out = model(x)
            loss = criterion(out, y)
        
        loss.backward()
        optimizer.step()

def eval_model(model, loader, device):
    model.eval()
    y_true, y_pred, y_prob = [], [], []
    with torch.no_grad():
        for x, y in loader:
            x = x.to(device)
            out = model(x)
            prob = torch.sigmoid(out).cpu().numpy()
            pred = (prob > 0.5).astype(int)
            y_true.extend(y.cpu().numpy())
            y_pred.extend(pred)
            y_prob.extend(prob)
    return np.array(y_true), np.array(y_pred), np.array(y_prob)

# 5. Metrics calculation
def calc_metrics(y_true, y_pred, y_prob):
    auc = roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else np.nan
    acc = accuracy_score(y_true, y_pred)
    sen = recall_score(y_true, y_pred, pos_label=1)
    spc = recall_score(y_true, y_pred, pos_label=0)
    return {'auc': auc, 'acc': acc, 'sen': sen, 'spc': spc}

# 6. Main loop
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
n_repeat = 20
n_folds = 3
batch_size = 32
epochs = 50

# Model selection configuration
MODEL_TYPES = {
    'SimpleDNN': SimpleDNN,
    'AttentionMLP': AttentionMLP,
    'AutoencoderClassifier': AutoencoderClassifier
}

# Select models to test (can test multiple simultaneously)
models_to_test = ['SimpleDNN', 'AttentionMLP', 'AutoencoderClassifier']  # Can modify here to select models

all_results = []

# Print model parameter count comparison
print_model_sizes(input_dim=50)  # Assume 50 features

for source_host in tqdm(host_map, desc="Source host"):
    source_keys = host_map[source_host]
    X_source = np.concatenate([profiles[k].values for k in source_keys])
    y_source = np.concatenate([labels[k] for k in source_keys])
    y_source = (y_source == 'sick').astype(int)
    target_hosts = [h for h in host_map if h != source_host]
    for target_host in tqdm(target_hosts, desc=f"Target host ({source_host})", leave=False):
        for target_key in tqdm(host_map[target_host], desc=f"Target dataset ({target_host})", leave=False):
            X_target = profiles[target_key].values
            y_target = (labels[target_key] == 'sick').astype(int)
            for model_name in tqdm(models_to_test, desc=f"Model ({target_key})", leave=False):
                model_class = MODEL_TYPES[model_name]
                for strategy in tqdm(['Regional', 'Transfer', 'Independent'], desc=f"Strategy ({model_name})", leave=False):
                    for repeat in tqdm(range(n_repeat), desc=f"Repeat ({strategy})", leave=False):
                        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=repeat)
                        for fold, (train_idx, test_idx) in enumerate(skf.split(X_target, y_target)):
                            # Check if it's an autoencoder model
                            is_autoencoder = model_name == 'AutoencoderClassifier'
                            
                            # Data preparation
                            if strategy == 'Regional':
                                # Train only on source domain, test on target domain
                                model = model_class(X_source.shape[1]).to(device)
                                optimizer = optim.Adam(model.parameters(), lr=1e-3)
                                criterion = nn.BCEWithLogitsLoss()
                                ds = TensorDataset(torch.tensor(X_source, dtype=torch.float32), torch.tensor(y_source, dtype=torch.float32))
                                loader = DataLoader(ds, batch_size=batch_size, shuffle=True)
                                for ep in range(epochs):
                                    train_model(model, loader, criterion, optimizer, device, is_autoencoder)
                                # Testing
                                test_ds = TensorDataset(torch.tensor(X_target[test_idx], dtype=torch.float32), torch.tensor(y_target[test_idx], dtype=torch.float32))
                                test_loader = DataLoader(test_ds, batch_size=batch_size)
                                y_true, y_pred, y_prob = eval_model(model, test_loader, device)
                            elif strategy == 'Transfer':
                                # First train on source domain, then fine-tune on target domain training set
                                model = model_class(X_source.shape[1]).to(device)
                                optimizer = optim.Adam(model.parameters(), lr=1e-3)
                                criterion = nn.BCEWithLogitsLoss()
                                ds = TensorDataset(torch.tensor(X_source, dtype=torch.float32), torch.tensor(y_source, dtype=torch.float32))
                                loader = DataLoader(ds, batch_size=batch_size, shuffle=True)
                                for ep in range(epochs):
                                    train_model(model, loader, criterion, optimizer, device, is_autoencoder)
                                # Fine-tuning
                                finetune_ds = TensorDataset(torch.tensor(X_target[train_idx], dtype=torch.float32), torch.tensor(y_target[train_idx], dtype=torch.float32))
                                finetune_loader = DataLoader(finetune_ds, batch_size=batch_size, shuffle=True)
                                for ep in range(epochs):
                                    train_model(model, finetune_loader, criterion, optimizer, device, is_autoencoder)
                                # Testing
                                test_ds = TensorDataset(torch.tensor(X_target[test_idx], dtype=torch.float32), torch.tensor(y_target[test_idx], dtype=torch.float32))
                                test_loader = DataLoader(test_ds, batch_size=batch_size)
                                y_true, y_pred, y_prob = eval_model(model, test_loader, device)
                            elif strategy == 'Independent':
                                # Train only on target domain training set
                                model = model_class(X_target.shape[1]).to(device)
                                optimizer = optim.Adam(model.parameters(), lr=1e-3)
                                criterion = nn.BCEWithLogitsLoss()
                                train_ds = TensorDataset(torch.tensor(X_target[train_idx], dtype=torch.float32), torch.tensor(y_target[train_idx], dtype=torch.float32))
                                train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
                                for ep in range(epochs):
                                    train_model(model, train_loader, criterion, optimizer, device, is_autoencoder)
                                # Testing
                                test_ds = TensorDataset(torch.tensor(X_target[test_idx], dtype=torch.float32), torch.tensor(y_target[test_idx], dtype=torch.float32))
                                test_loader = DataLoader(test_ds, batch_size=batch_size)
                                y_true, y_pred, y_prob = eval_model(model, test_loader, device)
                            # Record results
                            metrics = calc_metrics(y_true, y_pred, y_prob)
                            all_results.append({
                                'source_host': source_host,
                                'target_host': target_host,
                                'target_dataset': target_key,
                                'model': model_name,
                                'strategy': strategy,
                                'repeat': repeat,
                                'fold': fold,
                                **metrics
                            })

# 7. Save all results
results_df = pd.DataFrame(all_results)
results_df.to_csv('transfer_learning_results_with_attention.csv', index=False)

# 8. Print average metrics
print("="*80)
print("TRANSFER LEARNING RESULTS COMPARISON")
print("="*80)

# Overall performance grouped by model and strategy
overall_summary = results_df.groupby(['model', 'strategy'])[['auc', 'acc', 'sen', 'spc']].mean().reset_index()
print("\nOverall Performance by Model and Strategy:")
print(overall_summary.round(4))

# Model comparison (average across all strategies)
model_comparison = results_df.groupby('model')[['auc', 'acc', 'sen', 'spc']].mean().reset_index()
print("\nModel Comparison (Average across all strategies):")
print(model_comparison.round(4))

# Detailed comparison by strategy
print("\n" + "-"*60)
for strategy in ['Regional', 'Transfer', 'Independent']:
    print(f"\n{strategy} Strategy - Model Comparison:")
    strategy_results = results_df[results_df['strategy'] == strategy].groupby('model')[['auc', 'acc', 'sen', 'spc']].mean()
    print(strategy_results.round(4))

# Statistical significance comparison (simple performance difference analysis)
print("\n" + "-"*60)
print("Performance Improvement of AttentionMLP over SimpleDNN:")
for strategy in ['Regional', 'Transfer', 'Independent']:
    simple_auc = results_df[(results_df['model'] == 'SimpleDNN') & (results_df['strategy'] == strategy)]['auc'].mean()
    attention_auc = results_df[(results_df['model'] == 'AttentionMLP') & (results_df['strategy'] == strategy)]['auc'].mean()
    improvement = attention_auc - simple_auc
    print(f"{strategy}: {improvement:+.4f} AUC improvement")

print(f"\nDetailed results saved to: transfer_learning_results_with_attention.csv") 