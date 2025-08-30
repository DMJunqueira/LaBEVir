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

# Main execution
# Define search criteria
search_term <- "western equine encephalitis virus[organism]"
max_results <- 10000 # number of hits to be searched

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
  
  # Step 2: Fetching and saving sequences in FASTA format
  cat("Fetching sequences in FASTA format...\n")
  fasta_data <- entrez_fetch(db = "nuccore", id = sequence_ids, rettype = "fasta")
  
  # Naming sequences with the accession number
  fasta_lines <- str_split(fasta_data, "\n")[[1]]
  cleaned_fasta_lines <- if_else(
    str_detect(fasta_lines, "^>"),
    gsub("^>([^ ]+).*", ">\\1", fasta_lines),
    fasta_lines
  )
  cleaned_fasta_data <- paste(cleaned_fasta_lines, collapse = "\n")
  write_file(cleaned_fasta_data, file = fasta_output_path)
  cat(paste("Sequences saved to:", fasta_output_path, "\n"))
  
  # Step 3: searching and saving metadata
  cat("Fetching and parsing metadata...\n")
  
  # Usando `entrez_summary` para obter dados em uma estrutura de lista
  metadata_summary <- entrez_summary(db = "nuccore", id = sequence_ids)
  
  # Creating a dataframe
  metadata_list <- list()
  
  # Iterate over the results and extract metadata
  for (i in 1:length(metadata_summary)) {
    record <- metadata_summary[[i]]
    
    # Extracting accession, title and organism
    accession <- if (!is.null(record$accession)) record$accession else NA
    title <- if (!is.null(record$title)) record$title else NA
    organism <- str_extract(title, "^[^,]+")
    
    # Adding metadata to a temporary list
    metadata_list[[i]] <- data.frame(
      accession = accession,
      title = title,
      organism = organism,
      stringsAsFactors = FALSE
    )
  }
  
  # Combining all metadata in a unique dataframe
  metadata_df <- bind_rows(metadata_list)
  
  if (nrow(metadata_df) > 0) {
    write_tsv(metadata_df, metadata_output_path)
    cat(paste("Metadata saved to:", metadata_output_path, "\n"))
  } else {
    cat("No metadata could be extracted.\n")
  }
  
  cat("\nScript completed successfully.\n")
  
}, error = function(e) {
  cat("\nError:\n")
  cat(e$message, "\n")
})
