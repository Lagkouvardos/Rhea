#' Version 1.0
#' This script was last modified on 11/02/2016
#' Script: Taxonomic Binning
#'
#' Provides an overview of sample-specific relative abundances for all taxonomic levels
#'
#' Input:
#' 1. Set the path to the directory where the file is stored
#' 2. Write the name of the examined OTU file
#'
#' Output: The script is generating seven tab-delimited files
#' 1. Relative taxonomic abundance at the kingdom level for each sample
#' 2. Relative taxonomic abundance at the phyla level for each sample
#' 3. Relative taxonomic abundance at the class level for each sample
#' 4. Relative taxonomic abundance at the order level for each sample
#' 5. Relative taxonomic abundance at the family level for each sample
#' 6. Relative taxonomic abundance at the genera level for each sample
#' 7. Relative taxonomic abundance at all taxonomic levels for each sample
#' 
#' Graphical Output:
#' 8. Distribution of taxonomic relative abundances across all taxonomic groups for all samples
#'
#' Concept:
#' Taxonomic information is split into the different taxonomic levels for each OTU
#' Relative abundance values of OTUs belonging to the same taxonomic group are summed up
#' Individual taxonomic composition for each sample is generated
#'
#' Note:
#' If taxonomic information at a specific level is missing, the entry is replaced by
#' the last available taxonomic level including the prefix "unknown_"

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/numtax/)
#' Note: the path is denoted by forward slash "/"
setwd("D:/path/to/Rhea/4.Taxonomic-Binning") #<--- CHANGE ACCORDINGLY

#' Please give the file name of the OTU-table containing relative abundances and taxonomic classification 
otu_file<-"OTUs_Table-norm-rel-tax.tab"  #<--- CHANGE ACCORDINGLY


######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ###### 
##################################################################################

###################            Read input table              ####################

# Load the tab-delimited file containing the abundances and taxonomic information to be checked (rownames in the first column)
otu_table <-  read.table (otu_file,
                          check.names = FALSE,
                          header = TRUE,
                          dec = ".",
                          sep = "\t",
                          row.names = 1,
                          comment.char = "")

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Create a dataframe with a number of rows identical to the number of OTUs in the dataset
taxonomy <- otu_table[,dim(otu_table)[2]]

# Test if the taxonomy column is in the correct format (delimited by semicolon)
if(all(grepl("(?:[^;]*;){6}", taxonomy))==FALSE) {

#Send error message if taxonomy is not in the right format
  stop("Wrong number of taxonomic classes\n

Taxonomic levels have to be separated by semicolons (six in total). 
IMPORTANT: if taxonomic information at any level is missing, the semicolons are still needed:\n
       
      e.g.Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;
      e.g.Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;;")
} else { 

# Delete the taxonomy row from the OTU-table
otuFile <- otu_table[,c(1:dim(otu_table)[2] - 1)]

# Initialize empty dataframe
taxonomy_new <- NULL

for (i in 1:dim(otu_table)[1]) {
  # Split taxonomic information in its taxonomic classes
  # Kingdom - Phylum - Class - Family - Order - Genus
  splitTax <- strsplit(x = as.character(taxonomy[i]), ";")
  
  # Save the position where the first empty string (sequence of characters) occurs
  value <- which(splitTax[[1]] == "")[1]
  
  # Save the last known taxa information
  lastTaxa = splitTax[[1]][value - 1]
  
  # Replace all empty values by the last taxa information and the prefix "unkown_"
  splitTax <-replace(splitTax[[1]],splitTax[[1]] == "",paste("unknown_",lastTaxa))
 
  # Write new taxonomic information in the dataframe
  taxonomy_new[i] <- list(splitTax)
}

# Adjust dataframe with modified taxonomic information
taxonomy_new <- t(as.data.frame(taxonomy_new))
row.names(taxonomy_new) <- row.names(otuFile)

# Add level information to all taxonomies
# For taxonomies related to kingdom level
taxonomy_new[,1] <- sub("^","k__",taxonomy_new[,1])

# For taxonomies related to phylum level
taxonomy_new[,2] <- sub("^","p__",taxonomy_new[,2])

# For taxonomies related to class level
taxonomy_new[,3] <- sub("^","c__",taxonomy_new[,3])

# For taxonomies related to order level
taxonomy_new[,4] <- sub("^","o__",taxonomy_new[,4])

# For taxonomies related to family level
taxonomy_new[,5] <- sub("^","f__",taxonomy_new[,5])

# For taxonomies related to genus level
taxonomy_new[,6] <- sub("^","g__",taxonomy_new[,6])

#################################################################################

# Create list with taxonomic information for each taxonomy level
class_list <-
  list(
    unique(taxonomy_new[,1]),unique(taxonomy_new[,2]),
    unique(taxonomy_new[,3]), unique(taxonomy_new[,4]),
    unique(taxonomy_new[,5]),unique(taxonomy_new[,6])
  )

# Clone the created list for further processing
sample_list <- class_list
list_length <- NULL

# Iterate through all six taxonomy levels
for (a in 1:6) {
  
  lis <- lapply(class_list[a], lapply, length)
  names(lis)<-lapply(class_list[a],length)
  
  # Individual number of taxonomies for each taxonomic level
  num_taxa <- as.integer(names(lis))
  list_length[a] <- num_taxa
  
  # Iterate through taxonomic class specific taxonomies
  for (b  in 1:num_taxa) {
    
    # Initialize list with the value zero for all taxonomies
    sample_list[[a]][[b]] <- list(rep.int(0,dim(otuFile)[2]))
    
  }
}

#################################################################################
#################################################################################
# Save relative abundances of all samples for each taxonomy

# Iterate through all OTUs
for (i in 1:dim(otu_table)[1]) {
  
  # Iterate through all taxonomic levels
  for (m in 1:6) {
    
    # List of m-th taxonomies of i-th taxonomic levels (e.g. m = Kingdom, i = 4th OTU -> Clostridiales)
    taxa_in_list <- list(taxonomy_new[i,])[[1]][m]
    
    # Record the current position in a list
    position <- which(class_list[[m]] == taxa_in_list)
    
    # All rows with taxonomic information of n-th sample
    matrix <- data.matrix(otuFile)
    sub_sample_tax <-(subset(matrix,taxonomy_new[,m] == taxa_in_list))
    
    # Get the actual value out of the list (initialized with zero)
    temp <- unlist(sample_list[[m]][[position]])
    
    # Calculate the summed up relative abundances for the particular taxonomic class for n-th sample
    temp <- colSums(sub_sample_tax)
    
    # Replace values by new summed values
    sample_list[[m]][[position]] <- list(temp)
    
  }
}

#################################################################################
######                         Write output                                ######
#################################################################################

# Generate tables for each taxonomic class

##Kingdom table
# Create table with taxonomic information (kingdom level)
kingdom <-  matrix(unlist(sample_list[[1]]),nrow = dim(otuFile)[2],ncol = list_length[1],dimnames = list(names(otuFile),unlist(class_list[[1]])))
kingdom <- (t(kingdom))

##Phylum table
# Create table with taxonomic information (phylum level)
phyla <-matrix(unlist(sample_list[[2]]),nrow = dim(otuFile)[2],ncol = list_length[2],dimnames = list(names(otuFile),unlist(class_list[[2]])))
phyla <- (t(phyla))

# Order table according to taxonomic name (descending)
phyla <- phyla[order(row.names(phyla)),]

## Class table
# Create table with taxonomic information (class level)
classes <- matrix(unlist(sample_list[[3]]), nrow = dim(otuFile)[2], ncol = list_length[3], dimnames = list(names(otuFile),unlist(class_list[[3]])))
classes <- (t(classes))

# Order dataframe according to taxonomic name (descending)
classes <- classes[order(row.names(classes)),]

## Orders
# create table with taxonomic information (Order)
orders <-matrix(unlist(sample_list[[4]]),nrow = dim(otuFile)[2],ncol = list_length[4],dimnames = list(names(otuFile),unlist(class_list[[4]])))
orders <- (t(orders))

# Order dataframe according to taxonomic name (descending)
orders <- orders[order(row.names(orders)),]

## Family table
# Create table with taxonomic information (family level)
families <-matrix(unlist(sample_list[[5]]),nrow = dim(otuFile)[2],ncol = list_length[5],dimnames = list(names(otuFile),unlist(class_list[[5]])))
families <- (t(families))

# Order dataframe according to taxonomic name (descending)
families <- families[order(row.names(families)),]

## Genus level
# Create table with taxonomic information (generum level)
genera <- matrix(unlist(sample_list[[6]]),nrow = dim(otuFile)[2],ncol = list_length[6],dimnames = list(names(otuFile),unlist(class_list[[6]])))
genera <- (t(genera))

# Order dataframe according to taxonomic name (descending)
genera <- genera[order(row.names(genera)),]

# Merge all dataframes
tax_summary <-rbind.data.frame(kingdom,phyla,classes,orders,families,genera)

# Identify duplicates and remove them
tax_summary <- tax_summary[!duplicated(row.names(tax_summary)),]

################################################################################
######                        Write Output Files                           ######
#################################################################################

# Create a directory 
dir.create("Taxonomic-Binning")

# Set path for all outputs to the new directory
setwd("Taxonomic-Binning")

# Write output files for taxonomic composition of every sample
write.table(kingdom,"0.Kingdom.all.tab",sep = "\t",col.names = NA)
write.table(phyla,"1.Phyla.all.tab",sep = "\t",col.names = NA)
write.table(classes,"2.Classes.all.tab",sep = "\t",col.names = NA)
write.table(orders,"3.Orders.all.tab",sep = "\t",col.names = NA)
write.table(families,"4.Families.all.tab",sep = "\t",col.names = NA)
write.table(genera,"5.Genera.all.tab",sep = "\t",col.names = NA)
write.table(tax_summary,"tax.summary.all.tab",sep = "\t",col.names = NA)
suppressWarnings (try(write.table(tax_summary, "../../5.Serial-Group-Comparisons/tax.summary.all.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))

#################################################################################
######                        Write Graphical Output                       ######
#################################################################################

pdf("taxonomic-overview.pdf")
par(xpd=T, mar=par()$mar+c(0,0,0,9))

#Kingdom
#k_col=distinctColorPalette(dim(kingdom)[1])
k_col=rainbow(dim(kingdom)[1])
k_col=sample(k_col)
barplot(kingdom,col=k_col, cex.names=0.5, ylab="cumulative relative abundance (%)", las=2, main="Taxonomic binning at Kingdom level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(kingdom)),cex=0.7,col = rev(k_col),pch = 16,pt.cex = 1.2)

#Phyla
#p_col=distinctColorPalette(dim(phyla)[1])
p_col=rainbow(dim(phyla)[1])
p_col=sample(p_col)
barplot(phyla,col=p_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Phyla level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(phyla)),cex=0.7,col = rev(p_col),pch = 16,pt.cex = 1.2)

#Classes
c_col=rainbow(dim(classes)[1])
c_col=sample(c_col)
barplot(classes,col=c_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Class level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(classes)),cex=0.7,col = rev(c_col),pch = 16,pt.cex = 1.2)

#Orders
o_col=rainbow(dim(orders)[1])
o_col=sample(o_col)
barplot(orders,col=o_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Order level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(orders)),cex=0.7,col = rev(o_col),pch = 16,pt.cex = 1.2)

#Families
f_col=rainbow(dim(families)[1])
f_col=sample(f_col)
barplot(families,col=f_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Family level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(families)),cex=0.7,col = rev(f_col),pch = 16,pt.cex = 1.2)

#Genera
g_col=rainbow(dim(genera)[1])
g_col=sample(g_col)
barplot(genera,col=g_col, cex.names=0.5,ylab="cumulative relative abundance (%)",las=2, main="Taxonomic binning at Genus level")
legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(genera)),cex=0.7,col = rev(g_col),pch = 16,pt.cex = 1.2)

dev.off()
}

#################################################################################
######                           End of Script                             ######
#################################################################################
