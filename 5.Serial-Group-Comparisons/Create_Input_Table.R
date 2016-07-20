# This R script prepares two files (All TAXA and All OTUS) as inputs for the Serial-Group-Comparisons Script.

#************************
# A total of 4 files is required for merging. 
# 1. A file containing the alpha-diversity measures.
# 2. A file with the normalized relative abundances of OTUs across samples.
# 3. A file with the relative abundances of existing taxonomic groups.
# 4. A file with the categorical grouping of the samples and any additional meta-data.
#******************************

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of this script as the working folder (e.g. D:/imngs-toolbox/Rhea/differential-abundances")
#' Note: the path is denoted by forward slash "/"
setwd("D:/imngs-toolbox/Rhea/5.Serial-Group-Comparisons")  #<--- CHANGE ACCORDINGLY

#' Please give the name of the file with alpha-diversity measures
alpha <- "alpha-diversity.tab";                         #<--- CHANGE ACCORDINGLY

#' Please give the name of the file with OTUs relative abundance
RelativeAbundanceOTUs <- "OTUs_Table-norm-rel.tab";     #<--- CHANGE ACCORDINGLY

#' Please give the name of the file with relative abundances of different taxonomic levels
TaxanomyAll <- "tax.summary.all.tab";                   #<--- CHANGE ACCORDINGLY

#' Please give the name of the meta file with sample groups and additional metadata variables if available
MetaFile <- "mapping_file.tab";                         #<--- CHANGE ACCORDINGLY

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ###### 
##################################################################################

###################       Load all required libraries     ########################

# Check if required package GUniFrac is already installed
if (!require("compare")) {

# Package is not already installed and will be installed automatically
  install.packages("compare")
}
library(compare)

###################            Read input table              ####################

# Reading all files
alpha <- read.table(file=alpha,header=TRUE,sep="\t",row.names=1,check.names = F) 
RelativeAbundanceOTUs <-as.data.frame(t(read.table(file=RelativeAbundanceOTUs,header=TRUE,sep="\t",row.names = 1,check.names = F)))
MetaFile <- read.table(file=MetaFile,header=TRUE,sep="\t",comment.char = "",row.names = 1,check.names = F)
TaxanomyAll <- read.table(file=TaxanomyAll,header=TRUE,sep="\t",row.names=NULL,check.names = F)
ColnameTo_assign<- TaxanomyAll[,1]
TaxanomyAll[,1] <- NULL
TaxanomyAll <- as.data.frame(t(TaxanomyAll))
colnames(TaxanomyAll) <- ColnameTo_assign

######################          MAIN PROGRAM                #####################

# Preparing TAXA table:
combine_taxa<- cbind.data.frame(MetaFile[rownames(TaxanomyAll),],alpha[rownames(TaxanomyAll),],TaxanomyAll) # merging Meta+Alpha+Taxa based on same order row.
combine_taxa$SampleID <- row.names(combine_taxa)
combine_taxa<-combine_taxa[,c(ncol(combine_taxa),1:(ncol(combine_taxa)-1))]# replacing columnn at the beginging.

# Prepare OTU table:
combine_OTUs<- cbind.data.frame(MetaFile[rownames(RelativeAbundanceOTUs),],alpha[rownames(RelativeAbundanceOTUs),],RelativeAbundanceOTUs) # merging Meta+Alpha+OTUs based on same order row.
combine_OTUs$SampleID <- row.names(combine_OTUs)
colnames(combine_OTUs)
combine_OTUs<-combine_OTUs[,c(ncol(combine_OTUs),1:(ncol(combine_OTUs)-1))]# replacing columnn at the beginging.

#################################################################################
######                        Write Output Files                           ######
#################################################################################

# Writing tables
if(compareIgnoreOrder(row.names(TaxanomyAll),row.names((MetaFile)))$result & compareIgnoreOrder(row.names(alpha),row.names((MetaFile)))$result
   & compareIgnoreOrder(row.names(RelativeAbundanceOTUs),row.names((MetaFile)))$result){
  write.table(combine_taxa,file="TaxaCombined.tab",sep="\t",row.names=FALSE)
  write.table(combine_OTUs,file="OTUsCombined.tab",sep="\t",row.names=FALSE)
  print("Files are combined successfully")
}else{ 
  print("ATTENTION !!!! Sample names differ across files. Script aborted.")
  print("Please ensure that identical sample names are used.")
}

#################################################################################
######                           End of Script                             ######
#################################################################################