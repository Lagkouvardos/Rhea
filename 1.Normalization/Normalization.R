#' Version 1.6
#' This script was last modified on 11/02/2016
#' Script Task: Normalize OTU-tables
#' Author: Ilias Lagkouvardos
#'
#' Normalize abundance values of the input OTU table
#' Calculate relative abundances for all OTUs based on normalized values
#' 
#' Input: Please enter following parameters
#' 1. Set the path to the directory where the file is stored
#' 2. Write the name of the OTU-table of interest in quotes
#' 
#' Output: The script generates four tab-delimited files
#' 1. Normalized counts without taxonomy information
#' 2. Normalized counts with taxonomy information
#' 3. Normalized relative abundances without taxonomy information
#' 4. Normalized relative abundances with taxonomy information
#' 
#' Concept:
#' Normalization via division by the sum of sequences in a given sample
#' and mulitplication by the minimum sum across all samples
#' It is used instead of the classic rarefactioning approach
#' to prevent confounding effects of subsampling sequences with possible consequences on diversity
#' 
#' Note:
#' Files are stored in the current folder 
#' If a file is needed for downstream analysis, it is also automatically added to the appropriate folder
#' under the condition that original folder structure of Rhea is maintained.


##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/normalize/)
#' Note: the path is denoted by forward slash "/"
setwd("D:/Git/imngs-toolbox/Rhea/1.Normalization") #<--- CHANGE ACCORDINGLY

#' Please give the file name of the original OTU-table with taxonomic classification 
file_name<-"OTUs-Table.tab"  #<--- CHANGE ACCORDINGLY

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ###### 
##################################################################################

# Load the tab-delimited file containing the values to be be checked (rownames in the first column)
otu_table <-  read.table (file_name,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")

# Save taxonomy information in vector
taxonomy <- as.vector(otu_table$taxonomy)

# Delete column with taxonomy information in dataframe
otu_table$taxonomy <- NULL

# Calculate the minimum sum of all columns/samples
min_sum <- min(colSums(otu_table))

# Divide each value by the sum of the sample and multiply by the minimal sample sum
norm_otu_table <- t(min_sum * t(otu_table) / colSums(otu_table))

# Calculate relative abundances for all OTUs over all samples
# Divide each value by the sum of the sample and multiply by 100
rel_otu_table <- t(100 * t(otu_table) / colSums(otu_table))

# Re-insert the taxonomy information in normalized counts table
norm_otu_table_tax <- cbind(norm_otu_table,taxonomy)

# Reinsert the taxonomy information in relative abundance table
rel_otu_table_tax <- cbind(rel_otu_table,taxonomy)

#################################################################################
######                        Write Output Files                           ######
#################################################################################

# Write the normalized table in a file and copy in directories alpha-diversity and beta-diversity if existing
write.table(norm_otu_table, "OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(norm_otu_table, "../2.Alpha-Diversity/OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
suppressWarnings (try(write.table(norm_otu_table, "../3.Beta-Diversity/OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))

# Write the normalized table with taxonomy in a file
write.table(norm_otu_table_tax, "OTUs_Table-norm-tax.tab", sep = "\t",col.names = NA, quote = FALSE)

# Write the normalized relative abundance table in a file and copy in directory Serial-Group-Comparisons if existing
write.table(rel_otu_table, "OTUs_Table-norm-rel.tab", sep = "\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table, "../5.Serial-Group-Comparisons/OTUs_Table-norm-rel.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))

# Write the normalized relative abundance with taxonomy table in a file and copy in directory Taxonomic-Binning if existing
write.table(rel_otu_table_tax, "OTUs_Table-norm-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table_tax, "../4.Taxonomic-Binning/OTUs_Table-norm-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))

#################################################################################
######                           End of Script                             ######
#################################################################################