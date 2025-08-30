###################################################
## GENBANK SEQUENCE AND METADATA RETRIEVAL
## Type: Sequence Retrieval
## Author: Dennis Maletich Junqueira
## Date: 2025-08-29
##
## Description:
## This script programmatically searches the NCBI Nucleotide database
## for sequences based on a user-provided term. It retrieves full
## sequence data and key metadata (e.g., accession, organism, country,
## and host) and saves the results into a FASTA file and a CSV file.
##
## USAGE: R
###################################################
# Setting working directory
setwd("../Desktop/")

# Loading libraries
library(rentrez)
library(stringr)
library(readr)
library(dplyr)
library(seqinr)

# Function 1: search GenBank and get IDs
#' Search GenBank for a term and return a list of sequence IDs.
#'
#' @param db The database to search (e.g., 'nuccore').
#' @param term The search query (e.g., "coronavirus uruguay").
#' @param retmax The maximum number of IDs to retrieve.
#' @return A character vector of GenBank IDs.
search_and_fetch_ids <- function(db, term, retmax) {
  cat(paste("Searching GenBank for term:", term, "\n"))
  search_results <- entrez_search(db = db, term = term, retmax = retmax)
  cat(paste("Found", search_results$count, "entries.\n"))
  
  if (search_results$count == 0) {
    stop("No sequences found for this search term.")
  }
  
  return(search_results$ids)
}

# Function 2: fetch sequences and save as FASTA
#' Fetch sequences for a list of IDs and save them to a FASTA file.
#'
#' @param db The database (e.g., 'nuccore').
#' @param ids A vector of sequence IDs.
#' @param file_path The path to the output FASTA file.
#' @return A character vector of sequences.
fetch_sequences <- function(db, ids, file_path) {
  cat("Fetching sequences in FASTA format...\n")
  fasta_data <- entrez_fetch(db = db, id = ids, rettype = "fasta")
  
  # A função write_file do pacote readr é ideal para salvar texto
  write_file(fasta_data, file = file_path)
  cat(paste("Sequences saved to:", file_path, "\n"))
  
  return(fasta_data)
}

# Main execution
# Define search criteria
search_term <- "western equine encephalitis virus[organism]"
max_results <- 201 # number of hits to be searched

# Define output folder
output_dir <- "./"
if (!dir.exists(output_dir)) {
  cat(paste("Creating output directory:", output_dir, "\n"))
  dir.create(output_dir)
}

# Define output files
fasta_output_path <- file.path(output_dir, "sequences.fasta")
metadata_output_path <- file.path(output_dir, "metadata.tsv")

# Execute functions
tryCatch({
  # Step 1: Search for IDs
  sequence_ids <- search_and_fetch_ids(db = "nuccore", term = search_term, retmax = max_results)
  
  # Step 2: searching and saving sequences
  fetch_sequences(db = "nuccore", ids = sequence_ids, file_path = fasta_output_path)
  
  # Step 3: searching and saving metadata
  cat("Fetching and parsing metadata...\n")
  
  # Using `entrez_summary` to obtain data in a list structure
  metadata_summary <- entrez_summary(db = "nuccore", id = sequence_ids)
  
  # Creating an empty dataframe
  metadata_df <- data.frame(
    accession = character(),
    title = character(),
    organism = character(),
    stringsAsFactors = FALSE
  )
  
  # Iterate results and extracting metadata
  for (record in metadata_summary) {
    # Extracting accesion and title
    accession <- if (!is.null(record$uid)) record$uid else NA
    title <- if (!is.null(record$title)) record$title else NA
    
    # Extracting organism
    organism <- str_extract(title, "^[^,]+")
    
    # New row
    new_row <- data.frame(
      accession = accession,
      title = title,
      organism = organism,
      stringsAsFactors = FALSE
    )
    
    metadata_df <- bind_rows(metadata_df, new_row)
  }
  
  if (nrow(metadata_df) > 0) {
    write_tsv(metadata_df, metadata_output_path)
    cat(paste("Metadata saved to:", metadata_output_path, "\n"))
  } else {
    cat("No metadata could be extracted.\n")
  }
  
  cat("\nScript finalizado com sucesso.\n")
  
}, error = function(e) {
  cat("\nOcorreu um erro:\n")
  cat(e$message, "\n")
})

# Next steps:
#
# a) Gene identification using a reference sequence: https://github.com/DMJunqueira/LaBEVir/blob/main/002R108_Sequence_GeneIdentificationUsingReference.R
# b) Changing sequence names:

