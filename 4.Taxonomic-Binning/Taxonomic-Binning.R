#' Version 1.1
#' This script was last modified on 06/10/2026
#' Script: Taxonomic Binning
#'
#' Provides an overview of sample-specific relative abundances for all taxonomic levels
#'
#' Input:
#' 1. Set the path to the directory where the file is stored
#' 2. Write the name of the examined OTU file
#' 3. Number of taxonomic levels in the taxonomy column: 6 or 7 (default 7)
#' 4. Whether the taxonomy goes down to species level: TRUE or FALSE (default FALSE)
#'
#' Supported taxonomy formats (levels separated by ";", a final ";" is optional):
#'   6 levels: Kingdom;Phylum;Class;Order;Family;Genus
#'   7 levels: Domain;Kingdom;Phylum;Class;Order;Family;Genus   (new SILVA)
#'   with species = TRUE one more level is expected at the end: ...;Genus;Species
#' Missing information must still be separated by ";" (e.g. ...;Enterobacteriaceae;;)
#'
#' Output:
#' 1. Relative taxonomic abundance at each taxonomic level for each sample (tab-delimited)
#'    6 levels: 0.Kingdom, 1.Phyla, 2.Classes, 3.Orders, 4.Families, 5.Genera (.all.tab)
#'    7 levels: 0.Domain, 1.Kingdom, 2.Phyla, 3.Classes, 4.Orders, 5.Families, 6.Genera (.all.tab)
#'    species:  6.Species.all.tab (with 6 levels) or 7.Species.all.tab (with 7 levels)
#' 2. Relative taxonomic abundance at all taxonomic levels for each sample (tax.summary.all.tab)
#' 3. Only if the taxonomy column has problems: taxonomy-errors.tab (the script then stops)
#'
#' Graphical Output:
#' 4. Distribution of taxonomic relative abundances across all taxonomic groups for all samples
#'
#' Concept:
#' The taxonomy column is checked: every OTU must have the expected number of levels and
#' no empty level may be followed by a filled level. If any OTU fails, the problems are
#' printed on screen, written to taxonomy-errors.tab, and the script stops.
#' Taxonomic information is split into the different taxonomic levels for each OTU
#' Relative abundance values of OTUs belonging to the same taxonomic group are summed up
#' Individual taxonomic composition for each sample is generated
#'
#' Note:
#' If taxonomic information at a specific level is missing, the entry is replaced by
#' the last available taxonomic level including the prefix "unknown_"
#' Before the check: spaces around ";" are removed, the suffix "--other" is removed, and a
#' level whose whole value is "Incertae Sedis" is set to empty together with all deeper levels
#' (a message on screen reports how many OTUs were affected)
#' Species names are written as "Genus species" (e.g. Enterococcus faecalis)
#' Non-breaking spaces are replaced by normal spaces. Any other character that is not valid
#' UTF-8 (e.g. text saved with Windows encoding) is an error, shown as <xx> in taxonomy-errors.tab


##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/numtax/)
#' Note: the path is denoted by forward slash "/"
setwd("D:/path/to/Rhea/4.Taxonomic-Binning") #<--- CHANGE ACCORDINGLY

#' Please give the file name of the OTU-table containing relative abundances and taxonomic classification 
otu_file<-"OTUs_Table-norm-rel-tax.tab"  #<--- CHANGE ACCORDINGLY

#' Number of taxonomic levels before species: 6 (Kingdom...Genus) or 7 (Domain...Genus, SILVA 144 or newer)
#' For further details read introduction comments
n_levels <- 6  #<--- CHANGE ACCORDINGLY

#' Does the taxonomy go down to species level? (TRUE or FALSE)
species <- FALSE  #<--- CHANGE ACCORDINGLY

######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######

##################################################################################
######                             Main Script                              ###### 
##################################################################################

###################      Function: check taxonomy column      ####################

# ---------------------------------------------------------------------------
# check_taxonomy(): clean the taxonomy column and check that every OTU has the
# expected number of levels and no empty level followed by a filled level.
#   taxonomy   character vector (last column of the OTU table)
#   otu_ids    OTU names, same order
#   n_levels   6 or 7 (levels before species)
#   species    TRUE if a species level follows the genus
# Returns list(levels, fields [matrix OTU x level], errors [data.frame])
# ---------------------------------------------------------------------------
check_taxonomy <- function(taxonomy, otu_ids, n_levels, species) {
  if (!n_levels %in% c(6, 7)) stop("n_levels must be 6 or 7")
  if (!is.logical(species) || length(species) != 1 || is.na(species)) stop("species must be TRUE or FALSE")
  level_names <- if (n_levels == 7) c("Domain","Kingdom","Phylum","Class","Order","Family","Genus")
                 else c("Kingdom","Phylum","Class","Order","Family","Genus")
  if (species) level_names <- c(level_names, "Species")
  n_expected <- length(level_names)

  original <- as.character(taxonomy)
  original[is.na(original)] <- ""
  tax <- original

  # Non-breaking spaces (often copied from web pages or Excel) become normal spaces:
  # as a Windows byte in text that is not valid UTF-8, and as the UTF-8 character
  nbsp_byte <- !validUTF8(tax) & grepl("\xa0", tax, fixed = TRUE, useBytes = TRUE)
  tax[nbsp_byte] <- gsub("\xa0", " ", tax[nbsp_byte], fixed = TRUE, useBytes = TRUE)
  nbsp_utf8 <- validUTF8(tax) & grepl("\xc2\xa0", tax, fixed = TRUE, useBytes = TRUE)
  tax[nbsp_utf8] <- gsub("\xc2\xa0", " ", tax[nbsp_utf8], fixed = TRUE, useBytes = TRUE)
  n_nbsp <- sum(nbsp_byte | nbsp_utf8)

  # Any other character that is not valid UTF-8 is an error; in the error file
  # such characters are shown by their byte code, e.g. <e9>
  invalid <- !validUTF8(tax)
  shown <- original
  shown[invalid] <- iconv(tax[invalid], "latin1", "ASCII", sub = "byte")
  n_found_invalid <- lengths(strsplit(tax[invalid], ";", fixed = TRUE, useBytes = TRUE))
  tax[invalid] <- ""   # these OTUs are not processed further

  # Spaces around ';' and at both ends are removed (spaces inside names are kept)
  cleaned <- gsub("\\s*;\\s*", ";", trimws(tax))
  n_spaces <- sum(cleaned != tax)
  tax <- cleaned
  # The suffix "--other" is removed from any level (e.g. Pseudomonadati--other -> Pseudomonadati)
  no_other <- gsub("--other(?=;|$)", "", tax, perl = TRUE)
  n_other <- sum(no_other != tax)
  tax <- no_other

  # Levels per OTU: strsplit keeps empty levels and drops only the empty piece after
  # a final ';'  ->  "a;b;" = 2 levels, "a;b" = 2 levels, "a;;" = 2 levels
  fields  <- strsplit(tax, ";", fixed = TRUE)
  n_found <- lengths(fields)
  n_found[invalid] <- n_found_invalid

  # A level whose whole value is "Incertae Sedis" becomes empty, together with all deeper levels
  n_incertae <- 0
  fields <- lapply(fields, function(f) {
    hit <- which(tolower(f) == "incertae sedis")
    if (length(hit) > 0) { f[hit[1]:length(f)] <- ""; n_incertae <<- n_incertae + 1 }
    f })

  # Problems of each OTU (a row can have more than one)
  problem <- rep("", length(tax))
  add_problem <- function(idx, txt) problem[idx] <<- ifelse(problem[idx] == "", txt, paste(problem[idx], txt, sep = "; "))
  add_problem(which(invalid), "invalid character (not UTF-8, shown as <xx>)")
  add_problem(which(n_found != n_expected & !invalid), paste0("wrong number of levels (expected ", n_expected, ")"))
  gap <- sapply(fields, function(f) {
    empty <- f == ""
    any(empty) && any(!empty[which(empty)[1]:length(f)]) })
  add_problem(which(gap & !invalid), "empty level followed by a filled level")

  if (n_nbsp > 0)     message("Non-breaking spaces were replaced by normal spaces in ", n_nbsp, " OTU(s)")
  if (n_spaces > 0)   message("Spaces around ';' were removed in ", n_spaces, " OTU(s)")
  if (n_other > 0)    message("'--other' was removed in ", n_other, " OTU(s)")
  if (n_incertae > 0) message("'Incertae Sedis' was set to empty (with all deeper levels) in ", n_incertae, " OTU(s)")

  bad <- problem != ""
  errors <- data.frame(OTU = otu_ids[bad], taxonomy = shown[bad], levels_found = n_found[bad],
                       levels_expected = rep(n_expected, sum(bad)), problem = problem[bad])
  if (any(bad)) return(list(levels = level_names, fields = NULL, errors = errors, n_found = n_found))

  fields <- do.call(rbind, fields)
  dimnames(fields) <- list(otu_ids, level_names)
  # Species written as "Genus species" (a filled species always has a filled genus here)
  if (species) {
    has_sp <- fields[, "Species"] != ""
    fields[has_sp, "Species"] <- paste(fields[has_sp, "Genus"], fields[has_sp, "Species"])
  }
  list(levels = level_names, fields = fields, errors = errors, n_found = n_found)
}

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

###################       Check the taxonomy column          ####################

# Remove an error file left from a previous run
invisible(suppressWarnings(file.remove(file.path("Taxonomic-Binning", "taxonomy-errors.tab"))))

# Taxonomy is the last column
tax_check <- check_taxonomy(otu_table[,dim(otu_table)[2]], row.names(otu_table),
                            n_levels = n_levels, species = species)

# If any OTU has a problem: report it on screen and in taxonomy-errors.tab, then stop
if (nrow(tax_check$errors) > 0) {
  dir.create("Taxonomic-Binning", showWarnings = FALSE)
  write.table(tax_check$errors, file.path("Taxonomic-Binning", "taxonomy-errors.tab"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  counts <- table(tax_check$n_found)
  n_expected <- length(tax_check$levels)
  stop("\n", nrow(tax_check$errors), " of ", length(tax_check$n_found), " OTU(s) have problems in the taxonomy column.\n",
       "Expected ", n_expected, " levels (n_levels = ", n_levels, ", species = ", species, ").\n",
       "Number of levels found: ", paste(paste0(names(counts), " levels in ", counts, " OTU(s)"), collapse = ", "), ".\n",
       if (as.integer(names(counts)[which.max(counts)]) != n_expected)
         "Most OTUs do not have the expected number of levels: please check the settings n_levels and species.\n",
       "The OTUs to correct are listed in Taxonomic-Binning/taxonomy-errors.tab",
       call. = FALSE)
}

# Show the levels used, with the most complete OTU as an example
best <- which.max(rowSums(tax_check$fields != ""))
message("Taxonomic levels used: ", length(tax_check$levels), ". Example OTU ", rownames(tax_check$fields)[best], ":\n",
        paste(sprintf("  %-8s = %s", tax_check$levels, tax_check$fields[best, ]), collapse = "\n"))

# Delete the taxonomy column from the OTU-table
otuFile <- data.matrix(otu_table[, -dim(otu_table)[2], drop = FALSE])

# Names used for each taxonomic level
level_info <- data.frame(
  level  = c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"),
  prefix = c("d__",   "k__",    "p__",   "c__",  "o__",  "f__",   "g__",  "s__"),
  file   = c("Domain","Kingdom","Phyla", "Classes","Orders","Families","Genera","Species"),
  title  = c("Domain","Kingdom","Phyla", "Class","Order","Family","Genus","Species"))
level_info <- level_info[match(tax_check$levels, level_info$level),]
n_tax <- nrow(level_info)

###################     Replace missing levels by unknown_   ####################

taxonomy_new <- tax_check$fields
for (i in 1:dim(taxonomy_new)[1]) {
  splitTax <- taxonomy_new[i,]
  
  # Save the position where the first empty string (sequence of characters) occurs
  value <- which(splitTax == "")[1]
  
  # Save the last known taxa information
  lastTaxa = splitTax[value - 1]
  
  # Replace all empty values by the last taxa information and the prefix "unkown_"
  taxonomy_new[i,] <- replace(splitTax, splitTax == "", paste("unknown_",lastTaxa))
}

# Add level information to all taxonomies (e.g. "k__" for kingdom)
for (m in 1:n_tax) {
  taxonomy_new[,m] <- sub("^", level_info$prefix[m], taxonomy_new[,m])
}

#################################################################################
######          Sum up relative abundances for each taxonomic group        ######
#################################################################################

# One table per taxonomic level: rows = taxa, columns = samples
tax_tables <- list()
for (m in 1:n_tax) {
  # Rows (OTUs) of each taxonomic group, in order of first appearance
  groups <- split(1:dim(otuFile)[1], factor(taxonomy_new[,m], levels = unique(taxonomy_new[,m])))
  
  # Sum up the relative abundances of all OTUs of the same group, for each sample
  sums <- lapply(groups, function(rows) colSums(otuFile[rows, , drop = FALSE]))
  tab <- matrix(unlist(sums), ncol = dim(otuFile)[2], byrow = TRUE,
                dimnames = list(names(groups), colnames(otuFile)))
  
  # Order table according to taxonomic name (all levels except the first one)
  if (m > 1) tab <- tab[order(row.names(tab)), , drop = FALSE]
  tax_tables[[m]] <- tab
}

# Merge all tables
tax_summary <- do.call(rbind.data.frame, tax_tables)

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
for (m in 1:n_tax) {
  write.table(tax_tables[[m]], paste0(m - 1, ".", level_info$file[m], ".all.tab"), sep = "\t", col.names = NA)
}
write.table(tax_summary,"tax.summary.all.tab",sep = "\t",col.names = NA)
suppressWarnings (try(write.table(tax_summary, "../../5.Serial-Group-Comparisons/tax.summary.all.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))


#################################################################################
######                        Write Graphical Output                       ######
#################################################################################

pdf("taxonomic-overview.pdf")
par(xpd=T, mar=par()$mar+c(0,0,0,9))

# One page per taxonomic level
for (m in 1:n_tax) {
  tab <- tax_tables[[m]]
  t_col=rainbow(dim(tab)[1])
  t_col=sample(t_col)
  barplot(tab,col=t_col, cex.names=0.5, ylab="cumulative relative abundance (%)", las=2, main=paste0("Taxonomic binning at ", level_info$title[m], " level"))
  legend(par('usr')[2], par('usr')[4], bty='n',rev(row.names(tab)),cex=0.7,col = rev(t_col),pch = 16,pt.cex = 1.2)
}

dev.off()

#################################################################################
######                           End of Script                             ######
#################################################################################
