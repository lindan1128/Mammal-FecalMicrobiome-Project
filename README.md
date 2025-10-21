### Shared Fecal Microbiome Characteristics of Common Diseases in Domesticated Mammals



### Project structure
	Main/
   	├── data_process/                  # Codes for 16s data, WGS data  processing for each study 
	├── Rscripts/                      # Codes for analysis 
	│   ├── diversity.r                # alpha & beta diversity 
	│   ├── differential_microbes.r    # differential genera identification	
	│   ├── differential_microbes_batch_corrected.r       # differential genera identification with MMUPHin and ComBat correction
	│ 	├── maaslin2.r                 # MaasLin2 analysis
    │ 	├── aldex2.r                   # Aldex2 analysis
    │ 	├── ancombc2.r                 # AncomBC2 analysis
    │ 	├── netmoss2.r                 # Netmoss2 analysis
	├── machine_learning/              # Codes for machine learning
	│   ├── single_disease_prediction.py  # random forest, lasso regression and neural network for disease prediction
	│   ├── multi_diseases_prediction.py  # multi-class random forest and neural network for multi-disease prediction
	│   ├── transfer_learning.py       # transfer learning for cross-study generalizability evaluation
	├── README.md                      # Readme file
	
	
