#################################################################
## RENAMING SEQUENCES IN A FASTA FILE
## Type: Sequence
## Author: Dennis Maletich Junqueira
## Date: 2025-08-30
##
## DESCRIPTION:
## This script renames sequences in a FASTA file based on a
## provided mapping table in a CSV format.
## 
## USAGE: R
#################################################################
# Setting working directory
setwd("")

# Loading libraries
library(phylotools)

# Loading a .csv file
# First column should contain the header ID
# Second column should contain the new names under the header ID_New
data <- read.csv(file = "", header = TRUE, sep = ",")

# Load the fasta file containing the sequences to be renamed
infile_path <- ""

# Specify the path and name of the new fasta file
outfile_path <- ""

# Renaming
rename.fasta(infile = infile_path,
             ref_table = data,
             outfile = outfile_path)

rm(list=ls())
