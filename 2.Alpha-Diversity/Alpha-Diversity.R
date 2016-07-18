#' Version 1.3
#' Last modified on 22/02/2016
#' Script Task: Calculate alpha-diversity
#' Author: Ilias Lagkouvardos
#'
#' For meaningful species richness, use of normalized sequence counts is expected.
#' For richness calculation, only OTUs that are above 0.5 normalized counts are considered.

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the script as working folder (e.g D:/MyStudy/NGS/Rhea/alpha-diversity)
#' Note: the path is denoted by forward slash "/"
setwd("D:/imngs-toolbox/Rhea/2.Alpha-Diversity")    #<--- CHANGE ACCORDINGLY

#' Please give the file name of the normalized OTU-table (without taxonomic classification)
file_name <- "OTUs_Table-norm.tab"                #<--- CHANGE ACCORDINGLY

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                        Diversity Functions                           ###### 
##################################################################################

# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count=sum(x[x>0.5]^0)
  return(count)
}

# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total=sum(x)
  se=-sum(x[x>0]/total*log(x[x>0]/total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total=sum(x)
  se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
  return(se)
}

# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
{
  total=sum(x)
  si=sum((x[x>0]/total)^2)
  return(si)
}

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
  total=sum(x)
  si=round(1/sum((x[x>0]/total)^2),digits =2)
  return(si)
}

##################################################################################
######                             Main Script                              ###### 
##################################################################################

# Read a normalized OTU-table without taxonomy  
otu_table <- read.table (file_name, 
                       check.names = FALSE, 
                       header=TRUE, 
                       dec=".", 
                       sep = "\t",
                       row.names = 1)

# Order and transpose OTU-table
my_otu_table <- otu_table[,order(names(otu_table))] 
my_otu_table <-data.frame(t(my_otu_table))

# Apply diversity functions to table
otus_div_stats<-data.frame(my_otu_table[,0])
otus_div_stats$Richness<-apply(my_otu_table,1,Species.richness)
otus_div_stats$Shannon<-apply(my_otu_table,1,Shannon.entropy)
otus_div_stats$Shannon.effective<-apply(my_otu_table,1,Shannon.effective)
otus_div_stats$Simpson<-apply(my_otu_table,1,Simpson.concentration)
otus_div_stats$Simpson.effective<-apply(my_otu_table,1,Simpson.effective)

# Write the results in a file and copy in directory "Serial-Group-Comparisons" if existing
write.table(otus_div_stats, "alpha-diversity.tab", sep="\t", col.names=NA, quote=FALSE)
suppressWarnings (try(write.table(otus_div_stats[c(1,3,5)], "../5.Serial-Group-Comparisons/alpha-diversity.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))

##################################################################################
######                          End of Script                               ###### 
##################################################################################
