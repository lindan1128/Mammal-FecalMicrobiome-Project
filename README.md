# Mammal-FecalMicrobiome-Project

### Project structure
	Main/
	├── scripts/                      # Codes for analysis 
	│   ├── otu_batch_correction.r    # Covariates controlling 
	│   ├── host_specific_cooccurrence_network.r      # Host-specific co-occurrence network construction
	│   ├── bootstrapped_cooccurrence_network.r       # Bootstrapped co-occurrence network construction
	│ 	├── maaslin2.r                # Association analysis
    │ 	├── rf_prediction.r           # Random forest-based classification
    │ 	├── lasso_prediction.r        # Lasso-based classification
    │ 	├── alpha_diversity.r         # Alpha diversity
    │ 	├── permanova.r               # PERMANOVA diversity
   	├── 16s/                         # Codes for 16s data processing for each study 
   	├── mgs/                         # Codes for mgs data processing (assembly-free)
   	├── microbiomeHD/                # Human validation datasets extracted from microbiomeHD 
	├── Supplemental_Table/          # Folder for supplemental tables
	├── Supplemental_Figure/         # Folder for supplemental figures
	├── README.md                    # Readme file
	
