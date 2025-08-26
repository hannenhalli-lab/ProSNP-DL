Tool to Assess Features of Oncogenesis (OncoFeats)

Necessary files (included in the folder)
- tcgaAnalysisModules.R -  Main script with all functions and analyses
- GeneListsForFisherTest.Rda 
- mutationFrequencyListOfcancers.Rda  
- TCGAclinicalDataFile.Rda 
- TCGA_tpms.RDS (needs to be downloaded from https://drive.google.com/file/d/1H6q8OS_vEeRBY7qA0AeZ4dngZglmf9t3/view?usp=drive_link)
- gtex_tpms.rds (needs to be downloaded from https://drive.google.com/file/d/1Yu0PpsIhs6kJ5NOdK1tzw4MQhFEgByaf/view?usp=drive_link) 
- tcgaCNVdata.rds

Necessary Biowulf Modules to load
- R
	- to load R type out “module load R”

Input file specifications:
- The Input file should be a two column file
- Column 1 (Gene) - the gene can be in any of the following formats
	- ENSEMBL
	- ENTREZID
	- PFAM
	- REFSEQ
	- SYMBOL
	- UNIGENE
	- UNIPROT

- Column 2 (GeneType) - GeneTypes are arbitrarily assigned by the user, they are just predefined gene groups, there can be as many or as few as you want. 

Ex. Input - 6 genes with 3 labels
	
Gene    GeneType
TPM1    regen
MMP3    wound
LDLR    wound
ACTC1   stress
LECT2   regen
MYH11   stress 

Instructions for running the script 
- This code is meant to run on the terminal, initiate an R terminal session by typing ‘R’ 
- In your R terminal session type ‘source(tcgaAnalysisModules.R)’ to load necessary R libraries, code, and datasets.
- There are only 2 functions you need to run in order to do your analysis
	- loadTCGAdatasets()
		-  This function allows you to load your specific genes and select the cancer type and tissue type context for all analyses
		-  The function will give you prompts to choose from that makes navigating the analysis a little easier
		-  Once Everything is loaded, the function will tell you the percentage of genes that fail to map. This will not hinder any subsequent analysis
		**** When the code gives you options for available input methods, currently only “local” data is avalable
	- analyzeTCGAdata()
		- This function is meant to be run AFTER loadTCGAdatasets(), otherwise it will prompt an error
		- This function will let you choose what analysis you want to do with your set of genes and then output a file into your current directory.
		- Current supported analyses are
			1. Gene Expression
			2. Mutations
			3. Copy Number Variations
			4. Survival Analysis
			5. Fisher Test

**Since all output files are printed to the current directory, please make sure to move your files to your own directory so that the shared folder isn’t cluttered
   We will add a feature soon so that you can specify where you want all your output files to go


