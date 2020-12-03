#' Version 1.2
#' Last modified on 25/02/2016
#' Script: Correlation
#' Author: Ilias Lagkouvardos
#'
#' Calculate correlations between continuous meta-variables and taxonomic variables
#'
#' The script requires three obligatory actions from users (below, L41): 
#' 1. Set the path to the directory where the present script is stored
#' 2. Set the name of the input table, containing information on both continuous and taxonomic variables (created manually)
#' Note:  Row names are referring to samples
#'        The first columns correspond to the dependent continuous variables (e.g. metabolite concentrations)
#'        to be correlated with taxonomic variables
#'        These columns are followed by the taxonomic variables (i.e. normalized relative abundances for as many taxa as wanted)
#' 3. Starting position (column) of the taxonomic variables (e.g. start of OTU relative abundances) 
#'
#' Optional Input (below, L55):
#' 1. Significance cutoff
#' 2. Correlation analysis within-variable categories (e.g. OTUs against OTUs)
#' 3. Handling of missing values
#'
#' The script generates six tables (1-6) and two pdf files (7 and 8):
#' 1. Scaled log ratio transformed data
#' 2. Correlation matrix (Pearson)
#' 3. The asymptotic p-value
#' 4. Number of observations used for each pair of variables
#' 5. Pairwise information about correlation coefficients, p-value, number of observations, and corrected p-value
#' 6. Corrected table with pairwise information (p-value < 0.05)
#' 7. Graphical illustration of correlation matrix (with and without significance cutoff)
#' 8. Pairwise linear model for log-transformed data
#'
#' Concept:
#' Calculate Pearson correlation coefficients
#' Significant pairs are determined by p-value (default cutoff is 0.05)
#' Graphical output is generated if absolute correlation value is >0.5
#'
#' Note:
#' It is important to take into account the number of observations for the interpretation of results!

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the present script as the working folder (e.g. D:/studyname/NGS-Data/Rhea/correlation/)
#' Note: the path is denoted by forward slash "/"
setwd("D:/path/to/Rhea/6.Correlations")                     #<--- CHANGE ACCORDINGLY !!!

#' Please give the file name of the table containing the variables for analysis
input_file <-"OTUsCombined_Corr_input_table.tab"              #<--- CHANGE ACCORDINGLY !!!

#' Please give the position where the taxonomic variables (OTUs or taxonomic groups) start!!
#' IMPORTANT: Since the first column in the input file will be used as row names, we do not count it!
otu_variables_start <- 10                                        #<--- CHANGE ACCORDINGLY !!!

#################################################################################
#########         Optional parameters in this section                       #####
#################################################################################
# Unless users follow specific purposes and know exactly what to test, default values are recommended.

# Set the cutoff for significance
# Possible values are any real number between 0 and 1 (default is 0.05)
signf_cutoff <- 0.05

# Calculate correlation among taxonomic variables
# If selected, this will result in many additional tests for correlation among the taxonomic data
# Possible parameters are 1 or 0 (default is 0):
# 1 = calculate correlations within OTUs or taxa
# 0 = NO test within taxonomic variables
includeTax <- 0

# Calculate correlation among meta-variables
# If selected, this will result in many additional tests for correlation among the meta-data
# Possible parameters are 1 or 0 (default is 0):
# 1 = calculate correlations within meta-variables
# 0 = NO test within meta-variables
includeMeta <- 0

# Handling of missing values for meta-variables
# Possible parameters are 1 or 0 (default is 0):
# 1 = missing values are filled with the mean for the corresponding variable
# 0 = NO imputation (replacing missing data with substituted values)
fill_NA <- 0

# Treat zeros in taxonomic variables as missing values
# Possible parameters are 1 or 0 (default is 1):
# 1 = Consider taxonomic zeros as missing values
# 0 = Keep zeros for the calculation of correlations
replace_zeros <- 1

# Set a cutoff for the minimum number of values (prevalence) for a given taxonomic variable to be considered for calculation
# OTUs or taxa with prevalences below the cutoff are excluded from the analysis because considered as non-relevant in the study (very incidental occurancies)
# An OTU is considered present if the value is NOT missing or zero
# This filter reduces the number of tests for poorly supported correlations
# Possible values are any real number between 0 and 1 (default is 0.3, i.e. 30 % of samples must have a value for the given variable)
prevalence_exclusion <- 0.3

# Set a cutoff for the minimal number of pairs observations required for calculation of correlations
# The decision is subjective and depends on the total number of samples
# This filter reduces the number of tests for poorly supported correlations
# Possible values: any positive integer between 0 and the total number of samples (default is 4)
min_pair_support <- 4

# Set a significance cutoff for graphical output
# Only correlations with an uncorrected p-value less than the cuttoff will be plotted
# Possible values: any real number between 0 and 1 (default is 0.05)
plot_pval_cutoff <- 0.05

# Set a correlation coefficient cutoff for graphical output
# Only correlations above the cuttof (absolute value) will be ploted
# Possible values: any real number between 0 and 1 (default is 0.5)
plot_corr_cutoff <- 0.5

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ######
##################################################################################

###################       Load all required libraries     ########################

# Check if required packages are already installed, and install if missing
packages <-c("Hmisc","corrplot") 

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack)
  } 
}

# Applying the installation on the list of packages
lapply(packages, InsPack)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))

###################            Read input table              ####################
# Load the tab-delimited file containing the values to be checked (rownames in the first column)
my_data <-
  read.table (
    file = input_file,
    check.names = FALSE,
    header = TRUE,
    dec = ".",
    sep = "\t",
    row.names = 1,
    comment.char = ""
  )

# Clean table from empty lines
my_data <- my_data[!apply(is.na(my_data) | my_data=="",1,all),]
####################            Functions                  #####################

# Function for filling missing values with the mean of the column
fill_NA.mean <- function(vec)
{
  # Calculate mean value of each column (excluding missing values)
  m <- mean(vec, na.rm = TRUE)
  # Replace missing values with mean
  vec[is.na(vec)] <- m
  # Return the new input data frame
  return(vec)
}

# Function to logarithmically normalized OTU values
log_ratio <- function(data)
{
  # Compute the logarithmus
  log_data <- log(data)
  # Calculate exponential function of column-wise mean values for finite log transformed data
  gm <- exp(mean(log_data[is.finite(log_data)]))
  # Compute the logarithmus
  log_gm <- log(gm)
  # Take the difference of both log-transformed datasets
  data <- log_data - log_gm
  # Return the new OTU table
  return(data)
}

##################### END of FUNCTIONS ###########################

######################  MAIN PROGRAM #####################
my_data <- as.data.frame(apply(my_data,2,as.numeric))

first_OTU <- colnames(my_data)[otu_variables_start]

# Split the meta and taxonomic parts of the table
# Choose the continuous scaled variables
my_meta_data <- my_data[1:otu_variables_start - 1]

# Choose the taxonomic variables
my_otu_data <- my_data[otu_variables_start:dim(my_data)[2]]

# Process the meta measurements according to selection
if (fill_NA == 0) {
  # Do not do anything, just rename the file
  my_meta_fixed =  my_meta_data
}

if (fill_NA == 1) {
  # Fill non-zero missing meta-values with the mean of the column (optional)
  # Apply the previously implemented function 'fill_Na.mean' to the meta-data subset
  my_meta_fixed = apply(my_meta_data, 2, fill_NA.mean)
}

# The maximal number of absence
prevalence_cutoff <- dim(my_otu_data)[1] - (prevalence_exclusion*dim(my_otu_data)[1])

# Count how many missing values are found for each OTU 
na_count <-sapply(my_otu_data, function(y) sum(length(which(is.na(y)))))

# Count how many zeros are found for each OTU 
zero_count <-sapply(my_otu_data, function(y) sum(length(which(y==0))))
prevalence_count <- na_count + zero_count

# A new OTU-table is generated, where the number of missing values is below the set cutoff
my_otu_data <- my_otu_data[, prevalence_count <= prevalence_cutoff ]

# If the parameter is set, zeros are replaced with missing values
if (replace_zeros == 1) {
  my_otu_data[my_otu_data==0] <- NA
}

# Replace zeros with 0.0001 to avoid infinite number when calculating logarithmus (log(0)=-INF)
my_otu_data[my_otu_data==0] <- 0.0001

# Transform compositional data by log ratio transformation
my_otu_fixed = apply(my_otu_data, 2, log_ratio)

# Merge the meta- and OTU data in one table
transformed_data <- cbind(my_meta_fixed, my_otu_fixed)

# Centre and scale the values
my_scaled_data <- scale(transformed_data, center = TRUE, scale = TRUE)

# Calculate all pairwise correlations using Pearson correlation method
my_rcorr <- rcorr(as.matrix(my_scaled_data), type = "pearson")

# Generate vector with variable names
var_names <- row.names(my_rcorr$r)

# Depending on which parameters were set at the beginning, one query type is selected
# In each query type, three matrices are generated: p-value matrix, correlation matrix, support matrix
# All possibles pairs are saved in a vector (pairs)
if(includeTax==1 & includeMeta==0){

  # Correlation among OTUs and NO correlation among meta-variables
  row_names <- var_names[c(otu_variables_start:dim(my_rcorr$r)[1])]
  col_names <- var_names
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r[c(otu_variables_start:dim(my_rcorr$r)[1]),]
  my_pvl_matrix <-my_rcorr$P[c(otu_variables_start:dim(my_rcorr$P)[1]),]
  my_num_matrix <- my_rcorr$n[c(otu_variables_start:dim(my_rcorr$n)[1]),]
  
  # Set variable for plotting
  diagonale=0
  
} else if(includeTax==1 & includeMeta==1){
  
  # Correlation among OTUs and correlation among meta-variables
    row_names <-var_names
  col_names <- var_names
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r
  my_pvl_matrix <-my_rcorr$P
  my_num_matrix <- my_rcorr$n
  # Set variable for plotting
  diagonale=0
  
} else if (includeTax==0 & includeMeta==1) {
  # NO correlation among OTUs and correlation among meta-variables
  
  row_names <-var_names[c(1:(otu_variables_start - 1))]
  col_names <- var_names
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)),]
  my_pvl_matrix <-my_rcorr$P[c(1:(otu_variables_start - 1)),]
  my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)),]
  # Set variable for plotting
  diagonale=1
  
} else {
  # NO correlation among OTUs and NO correlation among meta-variables
  
  row_names <- var_names[c(1:(otu_variables_start - 1))]
  col_names <- var_names[otu_variables_start:dim(my_rcorr$r)[1]]
  pairs <-expand.grid(row_names, col_names)
  my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$r)[1])]
  my_pvl_matrix <-my_rcorr$P[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$P)[1])]
  my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$n)[1])]
  # Set variable for plotting
  diagonale=1
  
}

# Select the corresponding p-value for each pair 
p_vector <- as.vector(my_pvl_matrix)

# Select the corresponding correlation coefficient for each pair 
c_vector <- as.vector(my_cor_matrix)

# Select the corresponding number of observations for each pair 
n_vector <- as.vector(my_num_matrix)

# Generate matrix with the pairwise comparisons
my_pairs <-
  matrix(ncol = 5,
         c(
           as.character(pairs[, 2]),
           as.character(pairs[, 1]),
           c_vector,
           p_vector,
           n_vector
         ))

# Delete all pairs with insufficient number of pairs
my_pairs <- subset(my_pairs, as.numeric(my_pairs[,5]) > min_pair_support)

# Adjust p-value for multiple testing using the Benjamin-Hochberg method
pVal_BH <- round(p.adjust(my_pairs[,4], method = "BH"), 4)

# Add the corrected p-value in the table 
my_pairs <- cbind(my_pairs,as.numeric(pVal_BH))

# Remove similar pairs (values along the diagonal)
my_pairs <- my_pairs[!as.character(my_pairs[, 1]) == as.character(my_pairs[, 2]),]

# Remove duplicate pairs
my_pairs <- my_pairs[!duplicated(my_pairs[, 3]), ]

# Created matrix columns represent correlation coefficients, p-values, number of observations, and corrected p-values
# Rows represent the pairs
matrix_names <- list(c(rep("",times=dim(my_pairs)[1])),
                     c(
                       "variable1",
                       "variable2",
                       "correlation",
                       "pValue",
                       "support",
                       "Corrected"
                     ))

dimnames(my_pairs) <- matrix_names

# Create subset of pairs with significant p-values
my_pairs_cutoff <- my_pairs[as.numeric(my_pairs[, 4]) <= signf_cutoff, ]

# Convert to matrix
my_pairs_cutoff <- matrix(my_pairs_cutoff,ncol=6,dimnames = list(c(rep("",times=dim(my_pairs_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))

# Create subset of significant pairs with strong correlation (above 0.5)
my_pairs_cutoff_corr <- my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= 0.5, ]

# Remove columns containing no information
my_cor_matrix <- my_cor_matrix[, colSums(is.na(my_cor_matrix)) != nrow(my_cor_matrix)]

# Missing values in the correlation matrix are set to zero
my_cor_matrix[is.na(my_cor_matrix)] <- 0
#################################################################################
######                        Generate Graphs                              ######
#################################################################################
# Take current path in one variable to store results in seperate folders in further steps
OriginalPath <- getwd()

# Take the name of the inputfile to name the folder
prefix = paste(strsplit(input_file,"[.]")[[1]][1],sep="_")

# Make a directory name with inputfile name and date
newdir <- paste(prefix,Sys.Date(), sep = "_")

# Create a directory 
dir.create(newdir)

# Set path for all outputs to the new directory
setwd(newdir)
                    
# Save visualized correlation matrix in "corrplot.pdf"
pdf(file = "corrplot.pdf")

# Visualization of all correlations between meta-variables and OTUs
corrplot(
  matrix(data=na.omit(my_cor_matrix),nrow=dim(as.data.frame(row_names))[1], ncol = dim(as.data.frame(col_names))[1],dimnames=list(row_names,col_names)),
  tl.col = "black",
  tl.srt = 65,
  tl.cex = 0.6,
  cl.cex = 0.5,
  diag=diagonale
)

dev.off()

# Check if significance value for graphical output was modified by the user
if (plot_pval_cutoff != signf_cutoff | plot_corr_cutoff != 0.5) {
  # Generate a new matrix with the signficance cutoff 
  my_pairs_cutoff <- my_pairs[as.numeric(my_pairs[, 4]) <= plot_pval_cutoff, ]
  
  # Extract all significant pairs with the set correlation cutoff
  corr_pval_cutoff <- my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= plot_corr_cutoff, ]
  corr_pval_cutoff <- matrix(corr_pval_cutoff,ncol=6, dimnames=list(c(rep("",times=dim(corr_pval_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))
} else  {
  # If the significance cutoff is 0.05 an the correlation cutoff is 0.5 (for plotting)
  # Take the previously generated matrix
  corr_pval_cutoff <- my_pairs_cutoff_corr
  corr_pval_cutoff <- matrix(corr_pval_cutoff,ncol=6, dimnames=list(c(rep("",times=dim(corr_pval_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))
 }

# Save linearized transformed correlations of significant pairs in "linear_sign_pairs.pdf"
pdf(file = "linear_sign_pairs.pdf",family="sans",fillOddEven=TRUE)

# Iterate through all significant pairs
for (i in 1:dim(corr_pval_cutoff)[1]) {
  # Save log-scaled transformed values of the first variable of the pair
  x_df <- transformed_data[names(transformed_data) %in% corr_pval_cutoff[i, 2]]
  
  # Save as numerical vector
  x <- x_df[, 1]
  
  # Save log-scaled transformed values of the second variable of the pair
  y_df <- transformed_data[names(transformed_data) %in% corr_pval_cutoff[i, 1]]
  
  # Save as numerical vector
  y <- y_df[, 1]
  
  # Create a linear model for the pair (excluding missing values)
  clm <- lm(y ~ x, na.action = na.exclude)
  
  # Determine the number of steps for the generation of the confidence intervals
  steps <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/1000
  
  # Sequence of more finely/evenly spaced data than original one for the calculation of the confidence interval
  newx <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), steps)
  
  # Predict confidence interval based on the generated linear model
  # Prediction intervals are calculated based on the residuals of the regression equation 
  # Prediction intervals account for the variability around the mean response inherent in any prediction
  # It represents the range where a single new observation is likely to fall
  a <- predict(clm, newdata = data.frame(x = newx), interval = "confidence")
  
  # Graphic display of all log-scaled transformed values of the ith-pair
  plot(
    x,
    y,
    ylab="",
    xlab="", 
    cex.axis = 0.75 ,
    xaxt = 'n',
    yaxt = 'n',
    xlim=c(min(x,na.rm = TRUE),max(x, na.rm = TRUE)),
   xaxs="i"
  )
  title(xlab = names(x_df),line=0.5,font.lab=2,cex.lab=1.4)
  title(ylab = names(y_df),line=0.5,font.lab=2,cex.lab=1.4)
  
  # Draw the confidence interval around the fitted line
  polygon(c(newx,rev(newx)),c(a[,2],rev(a[,3])),col="grey91",border=TRUE, lty="dashed")
  
  # Samples are shown as dots
  points(x,y)
  
  # Draw linear regression line
  abline(clm, lwd=2)
  
  # Take calculated pairwise p-value from rcorr
  pvalue_text <- paste("P-value:", round(as.numeric(corr_pval_cutoff[i, 4]), 4), sep = "")
  
  # Take corrected pairwise p-value
  pvalue_corr_text <- paste("Adj. p-value:", as.numeric(corr_pval_cutoff[i, 6]), sep = "")
  
  # Take calculated pairwise correlation coefficient from rcorr
  corr_text <- paste("Pearson's r:", round(as.numeric(corr_pval_cutoff[i, 3]), 4), sep = "")
  
  # Take calculated pairwise number of observations from rcorr
  support_text <- paste("supported by ", round(as.numeric(corr_pval_cutoff[i, 5]), 4)," observations", sep = "")
  
  # Show correlation coefficient and p-value in the plot
  mtext(corr_text, side = 3, line = 2)
  mtext(pvalue_text, side = 3, line = 1)
  mtext(pvalue_corr_text, side = 3, line = 0)
  mtext(support_text, side = 1, line = 2)
  
  abline(par("usr")[3],0)
  segments(par("usr")[1],a[1,2],par("usr")[1],a[1,3])
  abline(par("usr")[4],0)
  
}

dev.off()

#################################################################################
######                        Write Output Files                           ######
#################################################################################
# Take current path in one variable to store results in seperate folders in further steps
OriginalPath <- getwd()

# Take the name of the inputfile to name the folder
prefix = paste(strsplit(input_file,"[.]")[[1]][1],sep="_")

# Make a directory name with inputfile name and date
newdir <- paste(prefix,Sys.Date(), sep = "_")

# Create a directory 
dir.create(newdir)

# Set path for all outputs to the new directory
setwd(newdir)
                    
# Write the log-scale transformed table
write.table(my_scaled_data,"transformed.tab",sep = "\t",col.names = NA,quote = FALSE)

# Write the correlation table
write.table(my_cor_matrix,"correlation-table.tab",sep = "\t",col.names = NA,quote = FALSE)

# Write the pvalue table
write.table(my_pvl_matrix,"pval-table.tab",sep = "\t",col.names = NA,quote = FALSE)

# Write the number of effective samples table
write.table(my_num_matrix,"support-table.tab",sep = "\t",col.names = NA,quote = FALSE)

# Write the significant correlations
write.table(my_pairs_cutoff,"cutoff-pairs-corr-sign.tab",sep = "\t",col.names = NA,quote = FALSE)

# Write plotted pairs
write.table(corr_pval_cutoff,"plotted-pairs-stat.tab",sep = "\t",col.names = NA,quote = FALSE)


if(!flag) { stop("
    It was not possible to install all required R libraries properly.
                 Please check the installation of all required libraries manually.\n
                 Required libaries:ade4, GUniFrac, phangorn, randomcoloR, Rcpp")
}


#################################################################################
######                           End of Script                             ######
#################################################################################
