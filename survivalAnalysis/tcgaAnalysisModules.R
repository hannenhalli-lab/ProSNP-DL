#source('FunctionsForTCGA.R')
generateGeneFile = F
if (generateGeneFile) {
	library("msigdbr")
	database = "msigdb"
	if (database == "msigdb") {
	  geneSet <- msigdbr(species = "human", category = "C6") #subcategory = c("IMMUNESIGDB")
	  t2g <- geneSet %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
	  t2n <- geneSet %>% dplyr::distinct(gs_name, gs_id) %>% as.data.frame()
  
	  downRegulated <- grep("DN", t2g$gs_name)
	  upRegulated <- grep("UP", t2g$gs_name)
  
	  upregulated <- data.frame(GeneName = t2g[upRegulated, 'gene_symbol'], GeneType = "Upregulated") %>% unique()
	  downregulated <- data.frame(GeneName = t2g[downRegulated, 'gene_symbol'], GeneType = "Downregulated") %>% unique()
  
	  uniqueDown <- downregulated[!(downregulated$GeneName %in% upregulated$GeneName), ]
	  uniqueUP <- upregulated[!(upregulated$GeneName %in% downregulated$GeneName), ]
	  InternalInputFile <- rbind(uniqueDown[1:100, ], uniqueUP[1:100, ]) %>% data.frame()
	}
	geneFile <- InternalInputFile
}


library(ensembldb, quietly = TRUE); library(EnsDb.Hsapiens.v86, quietly = TRUE); EnsdB <- EnsDb.Hsapiens.v86
Transcripts <- transcripts(EnsdB, columns = c("seq_name", "gene_id", "gene_name", "tx_id", "tx_biotype", "gene_biotype"), return.type = "DataFrame")
					#listColumns(EnsdB , "tx")), filter = TxBiotypeFilter("nonsense_mediated_decay"),
GenesDf <- unique(data.frame(GeneID = Transcripts$gene_id, GeneName = Transcripts$gene_name, GeneType = Transcripts$gene_biotype))

# #####################################################SNV analysis#########################################################

promptFunction <- function(vectorList) {
	userInput <- NA
	while(is.na(userInput)) {
		userInput <- readline(prompt = cat(blue$bold("Select one of the number : ")))
		userInput <- as.integer(userInput)
	}
	Sys.sleep(1)
	cat(red$bold(glue("\n\nyou have selected -> {vectorList[userInput]}\n")), "\n")
	return(userInput)
}

messageFunction <- function(vectorList, class, level) {
	message(glue("Avialable {class} {level} include:"))
	cat(red$bold(vectorList), sep = "\n")
}

outputFilePop <- function(fileNameString) {
	message(sprintf("writing output to: %s\n", fileNameString))
}

combineFiles <- function(File1, File2) {
	commonGenes <- rownames(File1)[rownames(File1) %in% rownames(File2)]
	File1 <- File1[rownames(File1) %in% commonGenes, ] %>% .[order(rownames(.)), ]
	File2 <- File2[rownames(File2) %in% commonGenes, ] %>% .[order(rownames(.)), ]
	combinedFile <- cbind(File1, File2) %>% data.frame()
	return(combinedFile)
}

message("loading datasets...")
tcgaExpression <- readRDS('TCGA_tpms.RDS')
GtexExpression <- readRDS('gtex_tpms.rds')
tcgaCNV <- readRDS('tcgaCNVdata.rds')
load(file = 'TCGAclinicalDataFile.Rda')
load('GeneListsForFisherTest.Rda')
load('mutationFrequencyListOfcancers.Rda')


supportedIDs <- c("ENSEMBL", "ENTREZID", "PFAM", "REFSEQ", "SYMBOL", "UNIGENE", "UNIPROT")
loadTCGAdatasets <- function(inputMethods = inputTypes, clinicalDf = clinicalData, tcgaList = tcgaExpression, gtexList = GtexExpression) {
	messageFunction(vectorList = inputMethods, class = "input", level = "methods")
	inputMethod <- promptFunction(inputMethods)
	#message(glue("Avialable {class} {level} include:"))
	#return(inputMethod)
	if (inputMethod == 2) {
	  
	  message("please enter your data frame containing gene names and gene labels\n")
	  geneFilePath <- NA
	  while(is.na(geneFilePath)) {
	    geneFilePath <- readline(prompt = cat(blue$bold("enter the full path to tab delimited gene file without quotes: ")))
		if (geneFilePath == "") {
			geneFilePath <- NA
		} else {
		    geneFilePath <- as.character(geneFilePath)
		    geneFile <- read.table(geneFilePath, sep = "\t", header = T)
			colnames(geneFile) <- c("Gene", "GeneType")
			message(glue("supported gene ids are:"))
			cat(red$bold(supportedIDs), sep = "\n")
			
			geneID <- NA
			while (is.na(geneID)) {
				geneID <- readline(prompt = cat(blue$bold("which gene ids your file has (please type one from the display above): "))) %>% toupper
				if(geneID == "") {
					geneID <- NA
				} else {
					geneID <- geneID %>% as.character()
					if (geneID == "ENSEMBL") {
						colnames(geneFile) <- c("ENSEMBL", "GeneType")
					} else {
						convertedIDs <- idConversion(dataForTest = geneFile[ ,1], from = geneID, to = "ENSEMBL", dataBase = orgDb, removeNA = T, returnType = "dataFrame")
						geneFile <- merge(geneFile, convertedIDs, by.x = colnames(geneFile)[1], by.y = colnames(convertedIDs)[1]) %>% subset(., select = c("ENSEMBL", "GeneType"))
					}
				}
			}
			
		}
		#print(geneFilePath)
	  }
   
	  cancerType <- NA
	  while(is.na(cancerType)) {
	    cancer_names <- paste(gsub("TCGA-","",names(tcgaList)), collapse=" | ")
	    cancerType <- readline(prompt = cat(blue$bold(paste("please enter the cancer type: ",cancer_names, " "))))
		if (cancerType == "") {
			cancerType <- NA
		} else {
		    clinical <- clinicalDf %>% subset(., type == cancerType) #[clinicalData$type %in% cancerType, ]
		    cancerType <- glue("TCGA-{cancerType}")
		    cancerExpression <- tcgaList[[cancerType]]
			cancerSamples <- TCGAquery_SampleTypes(names(cancerExpression), typesample = c("TP"))
			cancerExpression <- cancerExpression[, names(cancerExpression) %in% cancerSamples] %>% .[ , order(names(.))]
			colnames(cancerExpression) <- gsub("-\\w+$", "", colnames(cancerExpression) )
			cancerExpression <- cancerExpression[ ,!duplicated(colnames(cancerExpression))]
			colnames(cancerExpression) <- gsub("\\.", "-", colnames(cancerExpression))
			rownames(cancerExpression) <- gsub("\\..+", "", rownames(cancerExpression))
			cancerCNV <- tcgaCNV[[cancerType]]
			
			mutationNames <- names(mutationList)
			if (sum(grepl(cancerType, names(mutationList))) == 1) {
				cancerMutation <- mutationList[[cancerType]]
			}else {
				cancerMutation <- "mutation data not present for this cancer type"
			}
		}
	  }
	  sampleGroups <- NA
	  while (is.na(sampleGroups)) {
	  	sampleGroups <- readline(prompt = cat(blue$bold("do you have sample groups (yes/no): "))) %>% tolower
		if (sampleGroups == "") {
			 sampleGroups <- NA
		} else {
			if (sampleGroups == "yes") {
			    sampleGroupFile <- readline(prompt = cat(blue$bold("enter the full path to tab delimited sample group file without quotes: ")))
			    sampleGroupFile <- as.character(sampleGroupFile)
			    sampleGroupFile <- read.table(sampleGroupFile, sep = "\t", header = T)
				
			} else {
				sampleGroupFile <- "No sample groupings were provided"
			}
		}
	  }
	  tissueType <- NA; 
	  while (is.na(tissueType)) {
		 tissue_names <- paste(names(gtexList), collapse=" | ")	
	  	 tissueType <- readline(prompt = cat(blue$bold(paste("please enter one of the gtex tissue from the list (Brain, Liver etc, select 000 if not relevant): ",tissue_names, " "))))
		 if (tissueType == "") {
		 	tissueType <- NA
		 } else {
			 if(tissueType == "000") {	
				 normalExpression <- "No tissue was selected"
			 }else {
			 	normalExpression <- gtexList[[tissueType]]
				rownames(normalExpression) <- gsub("\\..+", "", rownames(normalExpression))
			 }
		 } 
	  }
	  inputList <<- list(cancerType = cancerType, tissueType = tissueType, geneFile = geneFile, sampleGroupFile = sampleGroupFile, clinical = clinical, 
		  	cancerMutation = cancerMutation, cancerCNV = cancerCNV, cancerExpression = cancerExpression, normalExpression = normalExpression)
	  #return()
	}
}

inputTypes <- c("1. download data", "2. local data")

# inputList <- loadTCGAdatasets(inputMethods = inputMethods)

library(crayon)
library(glue)


analysisTypes <- c("1. Gene Expression", "2. Mutations", "3. Copy Number Variations", "4. Survival Analysis", "5. Fisher Test")

analyzeTCGAdata <- function(data = inputList) {
	
	cancerExpression <- data[['cancerExpression']]
	normalExpression <- data[['normalExpression']]
	clinical <- data[['clinical']]
	geneFile <- data[['geneFile']]
	cancerType <- data[['cancerType']]
	tissueType <- data[['tissueType']]
	sampleGroupFile <- data[['sampleGroupFile']]
	cancerCNV <- data[['cancerCNV']]
	cancerMutation <- data[['cancerMutation']]
	#print(dim(cancerCNV))
	
	messageFunction(vectorList = analysisTypes, class = "analysis", level = "types")
	analysisType <- promptFunction(analysisTypes)
	
	if (analysisType == 1) {
		
		analysisMethods <- c("1. Correlation anlysis", "2. Fold change analysis")
					
		
		messageFunction(vectorList = analysisMethods, class = "analysis", level = "methods")
		analysisMethod <- promptFunction(analysisMethods)
		

		if (analysisMethod == 1) {
			optionVector <- c("1. Correlation coefficient in normal samples", "2. Correlation coefficient in cancer samples", "3. Correlation coefficient in both cases")
			messageFunction(vectorList = optionVector, class = "analysis", level = "options")
			analysisOption <- promptFunction(optionVector)
			if (analysisOption == 1) {
				fileNameString <- glue("{tissueType}.correlation.txt")
				correlationMatrix <- correlationFunction(normalExpression, geneFile = geneFile)
				message(glue("writing correlation matrix to {fileNameString}"))
				write.table(correlationMatrix, file = fileNameString, sep = "\t", quote = F)
			}
			if (analysisOption == 2) {
				fileNameString <- glue("{cancerType}.correlation.txt")
				correlationMatrix <- correlationFunction(cancerExpression, geneFile = geneFile)
				message(glue("writing correlation matrix to {fileNameString}"))
				write.table(correlationMatrix, file = fileNameString, sep = "\t", quote = F)
			}
			if (analysisOption == 3) {
				fileNameString1 <- glue("{tissueType}.correlation.txt")
				fileNameString2 <- glue("{cancerType}.correlation.txt")
				
				normalCorrelation <- correlationFunction(normalExpression, geneFile = geneFile)
				cancerCorrelation <- correlationFunction(cancerExpression, geneFile = geneFile)
				
				message(glue("writing correlation matrix to {fileNameString1}"))
				write.table(normalCorrelation, file = fileNameString1, sep = "\t", quote = F)
				message(glue("writing correlation matrix to {fileNameString2}"))
				write.table(cancerCorrelation, file = fileNameString2, sep = "\t", quote = F)
			}
		}
		if (analysisMethod == 2) {
			optionVector <- c("1. Normal vs. cancer", "2. Early vs. Late", "3. Custom Sample groups")
			messageFunction(vectorList = optionVector, class = "analysis", level = "options")
			analysisOption <- promptFunction(optionVector)
			
			if (is.data.frame(normalExpression) & is.data.frame(cancerExpression)) {
				combinedFile <- combineFiles(normalExpression, cancerExpression)
			} else {
				message("For fold change analysis, you must select a cancer as well as a relevant normal tissue...cell the function again")
			}
			  #combinedFile <- processExpressionFile(combinedFile, Transcripts = Transcripts, replicates = F)
			
			if (analysisOption == 1) {
				fileNameString <- glue("ExpressionSummaryFor{cancerType}andGtex-{tissueType}.txt")
				sampleGroupFile1 <- data.frame(sampleName = names(normalExpression), sampleType = "normal")
				sampleGroupFile2 <- data.frame(sampleName = names(cancerExpression), sampleType = "cancer")
				sampleGroupFile <- rbind(sampleGroupFile1, sampleGroupFile2) %>% data.frame()
			}
			
			if (analysisOption == 2) {
				fileNameString <- glue("ExpressionSummaryFor{cancerType}andGtex-{tissueType}.EarlyLate.txt")
				cat(red$bold("generating early vs late sample groupings...."))
				sampleGroupFile <- tcgaStaging(cancerType = cancerType, samples = colnames(cancerExpression), getSamples = F)
				sampleTable <- table(sampleGroupFile$sampleType)
				print(sampleTable)
				normalGroups <- data.frame(sampleName = colnames(normalExpression), sampleType = "normal")
				sampleGroupFile <- rbind(sampleGroupFile, normalGroups) %>% data.frame
			} 
			if (analysisOption == 3) { 
				fileNameString <- glue("ExpressionSummaryFor{cancerType}andGtex-{tissueType}.customSampleGroups.txt")
				if (is.data.frame(sampleGroupFile)) {
					
				} else {
					message("please provide the path to sample group file")
					sampleGroupFile <- NA
					while(is.na(sampleGroupFile)) {
					    sampleGroupFile <- readline(prompt = cat(blue$bold("enter the full path to tab delimited sample group file without quotes: ")))
						if (sampleGroupFile == "") {
							sampleGroupFile <- NA
						} else {
						    sampleGroupFile <- as.character(sampleGroupFile)
						    sampleGroupFile <- read.table(sampleGroupFile, sep = "\t", header = T)
						}
					}
				   
				}
				normalGroups <- data.frame(sampleName = colnames(normalExpression), sampleType = "normal")
				sampleGroupFile <- rbind(sampleGroupFile, normalGroups) %>% data.frame
			}
			expressionSummaryList <- expressionSummaryFunction(combinedExpressionFile = combinedFile, sampleGroupFile)
			expressionSummaryDf <- do.call("cbind", expressionSummaryList)
            
            outputMethods <- c("1. All Genes LFC output", "2. Input File LFC output")
					
		
            messageFunction(vectorList = outputMethods, class = "analysis", level = "methods")
            outputMethods <- promptFunction(outputMethods)
            
            if (outputMethods == 1) {
                outputfileNameString <- paste("AllGenes_",fileNameString,sep="")
                outputFilePop(outputfileNameString)
                write.table(expressionSummaryDf, file = outputfileNameString, sep = "\t", quote = F)
            } else {
                outputfileNameString <- paste("GeneSet_",fileNameString,sep="")
                outputFilePop(outputfileNameString)
                write.table(expressionSummaryDf[rownames(expressionSummaryDf) %in% inputList[["geneFile"]]$ENSEMBL,],
                            file = outputfileNameString, sep = "\t", quote = F)
            }
			
		}
	}
	if (analysisType == 2) {
		#cat(red$bold("Not available yet\n********\n"))
		fileNameString <- glue("mutationSummaryFor{cancerType}.txt")
		mutationDf <- merge(geneFile, cancerMutation, by.x = "ENSEMBL", by.y = "Gene", all.x = T)
		
		outputFilePop(fileNameString)
		write.table(mutationDf, file = fileNameString, sep = "\t", quote = F)
	}
	if (analysisType == 3) {
		fileNameString <- glue("CNVSummaryFor{cancerType}.txt")
		analysisMethods <- c("1. Summarize genome-wide CNVs", "2. Summarize CNVs for geneList")
		messageFunction(vectorList = analysisMethods, class = "analysis", level = "methods")
		analysisMethod <- promptFunction(analysisMethods)
		
		cnvSummaryList <- summarizeCNVs(cancerCNV)
		cnvSummaryDf <- bind_rows(cnvSummaryList, .id = "Gene")
		cNames <- colnames(cnvSummaryDf)
		
		colnames(cnvSummaryDf)[cNames == "-1"] <- "Lost"
		colnames(cnvSummaryDf)[cNames == "0"] <- "Unchanged"
		colnames(cnvSummaryDf)[cNames == "1"] <- "Gained"
		
		if (analysisMethod == 1) {
			outputFilePop(fileNameString)
			write.table(cnvSummaryDf, file = fileNameString, sep = "\t", quote = F)
		}
		if (analysisMethod == 2) {
			cnvSummaryDf <- merge(cnvSummaryDf, geneFile, by.x = "Gene", by.y = "ENSEMBL", all.x = T)
			cnvSummaryDf$GeneType <- as.character(cnvSummaryDf$GeneType)
			cnvSummaryDf$GeneType[is.na(cnvSummaryDf$GeneType)] <- "Other"
			outputFilePop(fileNameString)
			write.table(cnvSummaryDf, file = fileNameString, sep = "\t", quote = F)
		}
	}
	if (analysisType == 4) {
		analysisMethods <- c("1. Cox regression with median of genesets", "2. Cox regression with mean of genesets", 
				"3. Cox regression for each gene in geneset", "4. Log Rank")
		messageFunction(vectorList = analysisMethods, class = "analysis", level = "methods")
		analysisMethod <- promptFunction(analysisMethods)

		geneList <- geneFile[ ,1] %>% as.character()
		
		subExpression <- log(cancerExpression + 1) %>% data.frame(check.names = F)

		if(analysisMethod == 1){
		  clinical <- data[['clinical']]
		  fileNameString <- "coxRegressionUsingMedianExpressionOfGeneSets.txt"
		  subData <- merge(subExpression, geneFile, by.x = 0, by.y = "ENSEMBL") %>% .[ ,-1]
		  medianExpression <- subData %>% group_by(GeneType) %>% summarise_all(., median, na.rm = T) %>% 
		          data.frame(check.names = F) %>% set_rownames(.[ ,1]) %>% .[ ,-1] %>% t() %>% data.frame
		  dataToTest <- merge(clinical, medianExpression, by.x = "samples", by.y = 0)
		  SurvivalAnalysisList <- list()
		  for (colIndex in 9:ncol(dataToTest)) {
			  geneLabel <- colnames(dataToTest)[colIndex]
			  HR <- NA; Pval <- NA;
			  res <- list(HR = HR, Pval = Pval)
			    #dataToTest$BinnedExpression <- GetBins(dataToTest$cancerExpression, binStep = binStep)
			tryCatch({
					cox.base = coxph(Surv(time,status) ~ dataToTest[[colIndex]] + age, data = dataToTest)
					HR = summary(cox.base)$coefficients['dataToTest[[colIndex]]',2]
					Pval = summary(cox.base)$coefficients['dataToTest[[colIndex]]',5]
					res <- list(HR = HR, Pval = Pval)
			},error=function(x){})
			SurvivalAnalysisList[[geneLabel]] <- res
		  }
		  SurvivalAnalysisDf <- bind_rows(SurvivalAnalysisList, .id = "GeneType")
		  outputFilePop(fileNameString)
		  write.table(SurvivalAnalysisDf, file = fileNameString, sep = "\t", quote = F, row.names = F)
		}
		if(analysisMethod == 2){
			fileNameString <- "coxRegressionUsingMeanExpressionOfGeneSets.txt"
			subData <- merge(cancerExpression, geneFile, by.x = 0, by.y = "ENSEMBL") %>% .[ ,-1]
			meanExpression <- subData %>% group_by(GeneType) %>% summarise_all(., mean, na.rm = T) %>% 
			data.frame(check.names = F) %>% set_rownames(.[ ,1]) %>% .[ ,-1] %>% t() %>% data.frame
			dataToTest <- merge(clinical, meanExpression, by.x = "samples", by.y = 0)
			SurvivalAnalysisList <- list()
			for (colIndex in 9:ncol(dataToTest)) {
				geneLabel <- colnames(dataToTest)[colIndex]
				HR <- NA; Pval <- NA;
				res <- list(HR = HR, Pval = Pval)
				#dataToTest$BinnedExpression <- GetBins(dataToTest$cancerExpression, binStep = binStep)
				tryCatch({
				  cox.base = coxph(Surv(time,status) ~ dataToTest[[colIndex]] + age, data = dataToTest)
				  HR = summary(cox.base)$coefficients['dataToTest[[colIndex]]',2]
				  Pval = summary(cox.base)$coefficients['dataToTest[[colIndex]]',5]
				  res <- list(HR = HR, Pval = Pval)
				},error=function(x){})
				SurvivalAnalysisList[[geneLabel]] <- res
			}
			SurvivalAnalysisDf <- bind_rows(SurvivalAnalysisList, .id = "GeneType")
			outputFilePop(fileNameString)
			write.table(SurvivalAnalysisDf, file = fileNameString, sep = "\t",  quote = F, row.names = F)
		}
		if(analysisMethod == 3){
			fileNameString <- "SurvivalUsingcoxForGeneSet.txt"
			SurvivalAnalysisList <- list()
			SurvivalAnalysisList <- lapply(geneList, function(geneName) survivalFunction(geneName, expressionData = cancerExpression, clinicalData = clinical))
			SurvivalAnalysisDf <- do.call("rbind", SurvivalAnalysisList) %>% set_rownames(geneList)
			outputFilePop(fileNameString)
			write.table(SurvivalAnalysisDf, file = fileNameString, sep = "\t", quote = F, row.names = T)
		}
		if(analysisMethod == 4){
			fileNameString <- "SurvivalUsinglogRankForGeneSet.txt"
			SurvivalAnalysisList <- list()
			SurvivalAnalysisList <- lapply(geneList, function(geneName) survivalFunction(geneName, expressionData = cancerExpression, clinicalData = clinical, test = "logRank"))
			SurvivalAnalysisDf <- processLogrankOutput(SurvivalAnalysisList, geneLabels = geneList)
			outputFilePop(fileNameString)
			write.table(SurvivalAnalysisDf, file = fileNameString, sep = "\t", quote = F, row.names = F)
		}
	}
	if (analysisType == 5) {
		splitVector <- geneFile$GeneType
		geneList <- split(geneFile, splitVector) %>% lapply(., function(x) x[ ,1])
		geneTypes <- names(geneList)
		
		analysisMethods <- c("1. Fisher test against Hallmark gene sets", "2. Fisher test against cancerSEA gene sets")
		
		messageFunction(vectorList = analysisMethods, class = "analysis", level = "methods")
		analysisMethod <- promptFunction(analysisMethods)
		
		
		universeList <- GenesDf %>% subset(., GeneType == "protein_coding", select = 'GeneID') %>% unlist() %>% unique() %>% as.character()
		message("using all protein coding genes as background...")
		
		
		if (analysisMethod == 1) {
			fileNameString <- glue("ResultOfFisherTestwithHallmarks.txt")
			enrichmentResult <- lapply(hallMarkList, function(x) lapply(geneList, function(y) fisherTest(list1 = x, list2 = y, universe = universeList)))
			enrichmentResult <- lapply(enrichmentResult, function(x) bind_rows(x, .id = "geneType")) #%>% data.frame()
			enrichmentResult <- bind_rows(enrichmentResult, .id = "signature")
		}
		if (analysisMethod == 2) {
			fileNameString <- glue("ResultOfFisherTestwithcancerSEA.txt")
			enrichmentResult <- lapply(cancerSeaList, function(x) lapply(geneList, function(y) fisherTest(list1 = x, list2 = y, universe = universeList)))
			enrichmentResult <- lapply(enrichmentResult, function(x) bind_rows(x, .id = "geneType")) #%>% data.frame()
			enrichmentResult <- bind_rows(enrichmentResult, .id = "signature")
		}
		outputFilePop(fileNameString)
		write.table(enrichmentResult, file = fileNameString, sep = "\t", quote = F, row.names = F)
	}
}





################################Functions used for analysis##########################
library(dplyr); library(magrittr); 
library(preprocessCore);
library(glue); library(TCGAbiolinks); 
library(survival); library(OneR)
library(gplots)

library(org.Hs.eg.db)
library(clusterProfiler)
orgDb <- "org.Hs.eg.db"

idConversion <- function(dataForTest, from, to, dataBase = orgDb, removeNA = T, returnType = "dataFrame") {
	geneList <- dataForTest #names(dataForTest)
	idsDf <- bitr(geneList, fromType = from, toType = to, OrgDb = dataBase, drop = F)
	idsDf <- idsDf[!duplicated(idsDf[[from]]),] ###randomly retain unique ids from one to multiple comversion cases
	#
	if (removeNA) {
		idsDf <- idsDf[!(is.na(idsDf[[to]])), ]
		message("the cases where gene id conversion was not have been removed")
	}
	if (returnType == "dataFrame") {
		return(idsDf)
	} else {
		geneList <- idsDf[[to]]
		return(geneList)
	}
}


processLogrankOutput <- function (survResults, geneLabels = geneList) {
	lowMedian <- lapply(survResults, function(x) x[["groupMedians"]][2]) %>% unlist()
	highMedian <- lapply(survResults, function(x) x[["groupMedians"]][1]) %>% unlist()
	groupMediansData <- data.frame(lowMedian = lowMedian, highMedian = highMedian)
	medianDifference <- groupMediansData[ ,2] - groupMediansData[ ,1]
	indexes <- !is.na(medianDifference)
	groupMediansData <- groupMediansData[indexes, ]

	survResults <- survResults[indexes]
	geneLabels <- geneLabels[indexes]
	
	SurvFitData <- lapply(survResults, function(x) x[["fit"]])
	survDiffData <- lapply(survResults, function(x) x[["sd"]])
	
	#save(SurvFitData, SurvFitData, file = 'LogrankResultsForSurvivalWithDrosophilaLikeSetup.Rda')
	Pvalue <- lapply(survDiffData, function(x) (1 - pchisq(x$chisq, length(x$n) - 1))) %>% unlist() 
	EffectSize <- lapply(survDiffData, function(x) {ggVec = (x$obs / x$exp); ggVec[2] / ggVec[1]}) %>% unlist()
	#Direction <- lapply(survDiffData, function(x) {ggVec = (x$obs / x$exp); ggVec[2] / ggVec[1]}) %>% unlist()
	
	survivalResults <- data.frame(Gene = geneLabels, EffectSize, Pvalue) %>% cbind(., groupMediansData)
}

summarizeCNVs <- function(cnvData = cancerCNV) {
	cnvSummary <- apply(cnvData, 1, function(x) table(x)) #%>% data.matrix() %>% t())
	return(cnvSummary)
}


survivalFunction <- function(geneName, clinicalData, expressionData, test = "cox", binStep = 10, binned = F) {
	GeneIndex <- which(rownames(expressionData) == geneName)
	expressionLevel <- expressionData[GeneIndex, ] ##### the expression level should be named vector
	#HR <- NA; Pval <- NA;
	res <- NA
	expressionLevel <- unlist(expressionLevel)
	dataToTest <- merge(clinicalData, expressionLevel, by.x = "samples", by.y = 0)
	
	names(dataToTest) <- c(names(clinicalData),  "geneExpression")
	dataToTest$age = scale(dataToTest$age,center = T,scale = T)
	if (test == "logRank") {
		fit <- NA; sd <- NA; groupMedians <- c(NA, NA) 
		medianLevel <- median(dataToTest$geneExpression, na.rm = T, na.action = na.pass)
		dataToTest$BinnedExpression <- ifelse(dataToTest$geneExpression >= medianLevel, "HighExpression", "LowExpression")
		groupMedians <- dataToTest %>% group_by(BinnedExpression) %>% summarise(median = median(geneExpression, na.rm = T, na.action = na.pass))
		groupMedians <- na.omit(groupMedians)
		groupMedians <- groupMedians[ ,'median'] %>% unlist() %>% set_names(unlist(groupMedians[ ,'BinnedExpression']))
		res <- list(fit = fit, sd = sd, groupMedians = groupMedians)
		tryCatch({
			fit <- survfit(Surv(time, status) ~ BinnedExpression, data = dataToTest)
			sd <- survdiff(Surv(time, status) ~ BinnedExpression, data = dataToTest)
			res <- list(fit = fit, sd = sd, groupMedians = groupMedians)
	    },error=function(x){})
		return(res)
	}
	
	if (test == "cox") {
		HR <- NA; Pval <- NA;
		medianExpression = median(expressionLevel, na.rm = T, na.action = na.pass)
		meanExpression = mean(expressionLevel, na.rm = T, na.action = na.pass)
		res <- list(HR = HR, Pval = Pval)
		#dataToTest$BinnedExpression <- GetBins(dataToTest$geneExpression, binStep = binStep)
		tryCatch({
			if (binned) {
				dataToTest$BinnedExpression <- GetBins(dataToTest$geneExpression, binStep = binStep)
				cox.base = coxph(Surv(time,status) ~ BinnedExpression + age, data = dataToTest)
				HR = summary(cox.base)$coefficients['BinnedExpression',2]
				Pval = summary(cox.base)$coefficients['BinnedExpression',5]
			} else {
				cox.base = coxph(Surv(time,status) ~ geneExpression + age, data = dataToTest)
				HR = summary(cox.base)$coefficients['geneExpression',2]
				Pval = summary(cox.base)$coefficients['geneExpression',5]
			}
			
			res <- list(HR = HR, Pval = Pval, medianExpression = medianExpression, meanExpression = meanExpression)
		},error=function(x){})
		return(res)
	}
}

# ###sample names are of format. TCGA-XXX-XXX-XXX
# ###cancer type should be TCGA acronyms (for eg. LUAD, GBM, OV, BRCA, etc)

tcgaStaging <- function(cancerType, samples = NA, getSamples = T) {
  stageInformation <- "Not defined"
  if (grepl("-", cancerType)) {
  }else {
  	cancerType <- glue("TCGA-{cancerType}")
  }
  if (getSamples) {
	  query <- GDCquery(project = c(cancerType),  data.category = "Sequencing Reads",  experimental.strategy = "RNA-Seq")  
	  samples <- getResults(query, cols = c("file_name","cases", "file_id", "sample_type", "is_ffpe")) %>% subset(., is_ffpe == F, select = cases) %>% unlist()
	   samples <-  gsub("-\\w+-\\w+-\\w+-\\w+$", "", samples)
  } else {
  	 samples <- gsub("\\.", "-", samples)
  }
  #clinical <- GDCquery_clinic(paste("TCGA-",cancerType, sep = ""), type = "clinical")
  clinical <- GDCquery_clinic(cancerType, type = "clinical")
  if (is.null(clinical$ajcc_pathologic_stage)) {
    stageInformation = "Not defined"
    if (cancerType == "OV") {
      clinical <- data.frame(Barcode = clinical$submitter_id, Stage = clinical$figo_stage)#, MetaStasis = clinical$ajcc_pathologic_m)
      clinical$Stage <- gsub("[A-HJ-UW-Z]$", "", clinical$Stage)
      if (length(table(clinical$Stage)) > 0) {
        stageInformation = "defined"
      }
    }
  } else {
    clinical <- data.frame(Barcode = clinical$submitter_id, Stage = clinical$ajcc_pathologic_stage)#, MetaStasis = clinical$ajcc_pathologic_m)
    clinical$Stage <- gsub("[A-HJ-UW-Z]$", "", clinical$Stage)
    if (length(table(clinical$Stage)) > 0) {
      stageInformation = "defined"
    }
  }
  
  #samples <- gsub("-\\w+$", "", samples)
  #clinical <- clinical[match(samples, clinical$Barcode), ]
  samples <- data.frame(samples)
  clinical <- merge(clinical, samples, by.x = "Barcode", by.y = "samples")
  early = c("Stage I", "Stage II"); late = c("Stage III", "Stage IV")
  earlySamples <- clinical[clinical$Stage %in% early, 'Barcode'] %>% as.character()
  lateSamples <- clinical[clinical$Stage %in% late, 'Barcode'] %>% as.character()
  stageDf <- rbind(data.frame(sampleName = earlySamples, sampleType = "Early Cancer"), data.frame(sampleName = lateSamples, sampleType = "Late Cancer"))
  
  return(stageDf)
}


processExpressionFile <- function(dataFrame, Transcripts, replicates = F, normalize = T) {
  if (normalize) {
    print("performing quantile normalization")
  } else {
    print("processing without quantile normalization")
  }
  rownames(dataFrame) <- gsub("\\.\\d+", "", rownames(dataFrame))
  matches <- match(rownames(dataFrame), Transcripts$tx_id)
  Genes <- Transcripts[matches, 'gene_name']
  dataFrame <- cbind("GeneName" = Genes, dataFrame)
  GeneLevelTPM <- aggregate(.~GeneName, data = dataFrame, sum, na.rm = T, na.action = na.pass) 
  if (normalize) {
    GeneLevelTPM <- GeneLevelTPM[ ,-1] %>% as.matrix() %>% normalize.quantiles() %>% 
      set_rownames(GeneLevelTPM$GeneName) %>% 
      set_colnames(colnames(GeneLevelTPM)[-1]) %>% data.frame(check.names = F)
  }else {
    GeneLevelTPM <- GeneLevelTPM[ ,-1] %>% set_rownames(GeneLevelTPM$GeneName)
  }
  if (replicates) {
    print("you have passed replicates, taking averages")
    GeneLevelTPM <- GeneLevelTPM %>% t() %>% data.frame() %>%
      mutate(Stage = gsub("_Rep\\d+", "", colnames(GeneLevelTPM))) %>% 
      aggregate(.~Stage, data = ., mean, na.rm = T, na.action = na.pass) %>% 
      set_rownames(.$Stage) %>% .[,-1] %>% t() %>% data.frame() %>%
      set_rownames(rownames(GeneLevelTPM))
  }	else {
    print("no replicates are assumed in samples")
  }
  
  #GeneLevelTPM <- GeneLevelTPM[ ,mixedsort(names(GeneLevelTPM))]
  return(GeneLevelTPM)
}


expressionSummaryFunction <- function(combinedExpressionFile, sampleGroupFile) {
	combinedExpressionFile <- combinedExpressionFile %>% data.matrix() %>% normalize.quantiles() %>% data.frame() %>% 
				set_rownames(rownames(combinedExpressionFile)) %>% set_colnames(colnames(combinedExpressionFile))
				
	sampleList <- split(sampleGroupFile, sampleGroupFile$sampleType)
	sampleGroupFile$sampleName <- gsub("\\.", "-", sampleGroupFile$sampleName)
	colnames(combinedExpressionFile) <- gsub("\\.", "-", colnames(combinedExpressionFile))
	sampleList <- split(sampleGroupFile, sampleGroupFile$sampleType)
	expressionSummaryList <- list()
	for (sampleType in names(sampleList)) {
		samplesToUse <- sampleList[[sampleType]][ ,1] %>% as.character()
		subExpression <- combinedExpressionFile[ ,colnames(combinedExpressionFile) %in% samplesToUse]
		medianExpression <- apply(subExpression, 1, median, na.rm = T, na.action = na.pass)
		meanExpression <- apply(subExpression, 1, mean, na.rm = T, na.action = na.pass)
		tempDf <- cbind(meanExpression, medianExpression) %>% set_colnames(c(glue("mean.{sampleType}"), glue("median.{sampleType}")))
		expressionSummaryList[[sampleType]] <- tempDf
	}
	#foldChangeDf <- bind_rows(foldChangeList, .id = "sampleType")
	return(expressionSummaryList)
}

correlationFunction <- function(expressionFile, geneFile) {
	expressionFile <- expressionFile %>% data.matrix() %>% normalize.quantiles() %>% data.frame() %>% 
				set_rownames(rownames(expressionFile)) %>% set_colnames(colnames(expressionFile))
	#expressionFile <- log(expressionFile + 1)
	splitVector <- geneFile$GeneType
	geneList <- split(geneFile, splitVector)
	geneTypes <- names(geneList)
	subExpression1 <- expressionFile[rownames(expressionFile) %in% geneList[[geneTypes[1]]][ ,1], ] %>% t()
	subExpression2 <- expressionFile[rownames(expressionFile) %in% geneList[[geneTypes[2]]][ ,1], ] %>% t()
	
	subExpression1 <- log2(subExpression1 + 1)
	subExpression2 <- log2(subExpression2 + 1)
	
	correlationList <- list()
	correlationMatrix <- cor(subExpression1, subExpression2)
	return(correlationMatrix)
}


####do a fisher test####
fisherTest <- function(list1, list2, universe) {
	list1 <- list1[list1 %in% universe]
	list2 <- list2[list2 %in% universe]
	
	universe <- length(universe)
	c1 <- sum(list1 %in% list2)
	c2 <- length(list1) - c1
	c3 <- length(list2) - c1
	c4 <- universe - (c3 + length(list1))
	
	fisherMatrix <- matrix(c(c1,c2,c3,c4), nrow = 2)
	fisherTest <- fisher.test(fisherMatrix)
	enrichment <- fisherTest$estimate
	pvalue <- fisherTest$p.value
	
	return(c(enrichment = enrichment, pvalue = pvalue))
}

# ##
