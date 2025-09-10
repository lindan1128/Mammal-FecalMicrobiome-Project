# Mammal-FecalMicrobiome-Project

### Project structure
	Main/
	├── scripts/                      # Codes for analysis 
	│   ├── diversity.r               # alpha & beta diversity 
	│   ├── differential_microbes.r      # differential genera identification	
    │   ├── differential_microbes_batch_corrected.r       # differential genera identification with MMUPHin and limma correction
	│ 	├── maaslin2.r                # MaasLin2 analysis
    │ 	├── aldex2.r                  # Aldex2 analysis
    │ 	├── ancombc2.r                # AncomBC2 analysis
    │ 	├── netmoss2.r                # Netmoss2 analysis
    │ 	├── random_forest.r           # random forest modeling
    │ 	├── lasso.r                   # lasso regression modeling
    │ 	├── network.r                 # neural network modeling
	├── transfer_learning/            # Codes for transfer_learning
	│   ├── transfer_single_source.py  # transfer learning with single source study
	│   ├── transfer_multiple_sources.py  # transfer learning with multiple source studies
   	├── data_process/                         # Codes for 16s data, WGS data  processing for each study 
	├── Supplemental_Table/          # Folder for supplemental tables
	├── Supplemental_Figure/         # Folder for supplemental figures
	├── README.md                    # Readme file
	
