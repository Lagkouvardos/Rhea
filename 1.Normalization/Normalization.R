#' Version 2.0
#' This script was last modified on 20/01/2020
#' Script Task: Normalize OTU-tables
#' Author: Ilias Lagkouvardos
#' Contributions by: Thomas Clavel, Sandra Reitmeier
#'
#' Normalize abundance values of the input OTU table
#' Calculate relative abundances for all OTUs based on normalized values
#' Calculate rarefaction curves to help estimate sufficiency of sequencing depth for each sample
#'
#' Input: Please enter following parameters
#' 1. Set the path to the directory where the file is stored
#' 2. Write the name of the OTU-table of interest in quotes
#'
#' The script generates five tab-delimited files and one pdf file
#' 1. Normalized counts with taxonomy information
#' 2. Normalized counts without taxonomy information
#' 3. Normalized relative abundances with taxonomy information
#' 4. Normalized relative abundances without taxonomy information
#' 5. Rarefaction curves for all samples and the most undersampled ones (default 5 cases) as PDF
#' 6. Slope of the Rarefaction curve expressed as the number of species per 100 reads
#'
#' Concept:
#' The default method followed is normalization via division by the sum of sequences in a given sample
#' and multiplication by the minimum sum across all samples. It is used instead of the classic rarefaction approach
#' to avoid unnecessary variation due to the random subsampling and loss of information due to rounding.
#' The option of random subsampling is still available for normalization if deemed necessary by users (see Set parameters section below).
#' Rarefaction curves are showing species richness with respect to sequencing depth (number of reads).
#' The depth of sequencing may be considered too low for those samples with a rarefaction curve that does not reach a plateau at the available maximum number of reads.
#' This indicates that additional, less abundant species are probably in the sample
#' but were not detected at the available depth of sequencing.
#' The terminal slope of the curve for each sample is documented in the tab-delimited file
#' as the number of species added in richness by the last 100 reads

#' Note:
#' Files are stored in the current folder
#' If a file is needed for downstream analysis, it is also automatically added to the appropriate next folder
#' under the condition that the original folder structure of Rhea is maintained.
#' If the primary evaluation of sufficient sequencing depth led to subsequent removal of samples from the analysis,
#' the normalization script should be re-run with the new OTU table.

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/normalize/)
#' Note: the path is denoted by forward slash "/"
setwd("D:/path/to/Rhea/1.Normalization")      #<--- CHANGE ACCORDINGLY

#' Please give the file name of the original OTU-table with taxonomic classification
file_name <- "OTUs-Table.tab"                   #<--- CHANGE ACCORDINGLY

#' Please select the normalisation method
#' 0 = No random subsampling, no rounding
#' 1 = Random subsampling with rounding
method <- 0                                   #<--- CHANGE ACCORDINGLY

#' Pease select the normalization level used
#' 0 = Minimum sampling size
#' 1 = Fixed value (e.g. 1000)
level <- 0                                    #<--- CHANGE ACCORDINGLY

#' Please choose the value at which all samples will be normalized. (Only used if level selected is 1)
normCutoff <- 1000

#' Please choose the number of samples with the steepest rarefaction curves to be selectively plotted
#' The default number of samples presented separately is 5
labelCutoff <- 5                              #<--- CHANGE ACCORDINGLY

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ######
##################################################################################

###################       Load all required libraries     ########################

# Check if required packages are already installed, and install if missing
packages <- c("GUniFrac","vegan")

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos = "http://cloud.r-project.org/")
  }
}

# Applying the installation on the list of packages
lapply(packages, InsPack)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))

###################       Read all required input files      ####################

# Load the tab-delimited file containing the values to be be checked (rownames in the first column)
otu_table <-  read.table (file_name,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")


# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

####################       Normalize OTU Table          ###################


# Save taxonomy information in vector
taxonomy <- as.vector(otu_table$taxonomy)

# Delete column with taxonomy information in dataframe
otu_table$taxonomy <- NULL



if (level == 0) {
  # Calculate the minimum sum of all columns/samples
  min_sum <- min(colSums(otu_table))
} else {
  # The minimum size is set to a fixed reference level
  min_sum <- normCutoff
}


if (method == 0) {
  # Divide each value by the sum of the sample and multiply by the minimal sample sum
  norm_otu_table <- t(min_sum * t(otu_table) / colSums(otu_table))
} else {
  # Rarefy the OTU table to an equal sequencing depth
  norm_otu_table <- Rarefy(t(otu_table),depth = min_sum)
  norm_otu_table <- t(as.data.frame(norm_otu_table$otu.tab.rff))
}

# Calculate relative abundances for all OTUs over all samples
# Divide each value by the sum of the sample and multiply by 100
rel_otu_table <- t(100 * t(otu_table) / colSums(otu_table))


# Re-insert the taxonomy information in normalized counts table
norm_otu_table_tax <- cbind(norm_otu_table,taxonomy)

# Re-insert the taxonomy information in relative abundance table
rel_otu_table_tax <- cbind(rel_otu_table,taxonomy)

################################################################################
# Generate a two-sided pdf with a rarefaction curve for all samples and a curve
pdf(file = "RarefactionCurve.pdf")

# Plot the rarefaction curve for all samples
rarefactionCurve <- rarecurve(data.frame(t(otu_table)),
                              step = 20,
                              col = "black",
                              lty = "solid",
                              label = F,
                              xlab = "Number of Reads",
                              ylab = "Number of Species",
                              main = "Rarefaction Curves of All Samples")

# Generate empty vectors for the analysis of the rarefaction curve
slope = vector()
SampleID = vector()
angle <- c()

# Iterate through all samples
for (i in seq_along(rarefactionCurve)) {
  # If the sequencing depth is greater than 100, the slope of the line that passes between the last and last-100 count is calculated
  richness <- ifelse(length(rarefactionCurve[[i]]) > 6, 
                     (rarefactionCurve[[i]][length(rarefactionCurve[[i]])] - rarefactionCurve[[i]][length(rarefactionCurve[[i]])-5])/(attr(rarefactionCurve[[i]], "Subsample")[length(rarefactionCurve[[i]])]-attr(rarefactionCurve[[i]], "Subsample")[length(rarefactionCurve[[i]])-5]) , 1000)
  angle[i] <- ifelse(richness!=1000, atan(richness)*180/pi, NA)
  slope <- c(slope,richness)
  SampleID <- c(SampleID,as.character(names(otu_table)[i]))
}

# Generate the output table for rarefaction curve
curvedf <- cbind(SampleID, slope, angle)
ordered_vector <- order(as.numeric(curvedf[,2]), decreasing = TRUE)
curvedf <- curvedf[order(as.numeric(curvedf[,2]), decreasing = TRUE),]

# Generates a graph with all samples
# Underestimated cases are shown in red
for (i in 1:labelCutoff) {
  N <- attr(rarefactionCurve[[ordered_vector[i]]], "Subsample")
  lines(N, rarefactionCurve[[ordered_vector[i]]],col="red")
}

# Determine the plotting width and height
Nmax <- sapply(rarefactionCurve[ordered_vector[1:labelCutoff]], function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefactionCurve[ordered_vector[1:labelCutoff]], max)

# Creates an empty plot for rarefaction curves of underestimated cases
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Reads",
     ylab = "Number of Species", type = "n", main=paste(labelCutoff,"- most undersampled cases"))

for (i in 1:labelCutoff) {
  N <- attr(rarefactionCurve[[ordered_vector[i]]], "Subsample")
  lines(N, rarefactionCurve[[ordered_vector[i]]],col="red")
  text(max(attr(rarefactionCurve[[ordered_vector[i]]],"Subsample")), max(rarefactionCurve[[ordered_vector[i]]]), curvedf[i,1],cex=0.6)
}

dev.off()

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

# Write the rarefaction table
write.table(curvedf, "RarefactionCurve.tab", sep ="\t", quote = FALSE, row.names = FALSE)


# Error message
if(!flag) { stop("
                 It was not possible to install all required R libraries properly.
                 Please check the installation of all required libraries manually.\n
                 Required libaries: GUniFrac")
}

#################################################################################
######                           End of Script                             ######
#################################################################################
