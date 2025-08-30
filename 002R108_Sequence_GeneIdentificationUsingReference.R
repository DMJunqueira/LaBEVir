#################################################################
## GENE IDENTIFICATION USING A REFERENCE SEQUENCE
## Type: Sequence
## Author: Dennis Maletich Junqueira
## Date: 2025-08-27
##
## DESCRIPTION:
## This script aligns user-provided sequences to a reference genome of moderate length,
## identifies the corresponding gene for each sequence, and adds this
## information to a metadata file. It also generates specific
## sub-alignments for each gene found (output: *_alignments.fasta)
## 
## USAGE: terminal
## Rscript annotate_sequences.R --meta my_metadata.csv --aln my_seqs.fasta --ref_seq reference_genome.fasta --annot aanotation.csv --meta_out annotated_metadata.txt
## --meta = metadata file (.csv). It should contain column "Accession".
## --aln = sequence alignment (.fasta). Sequence names should start with "Accession_".
## --ref_seq = reference sequence (.fasta).
## --annot = reference sequence annotation (.csv). It should contain columns: gene_name,start,end
## --meta_out = newly annotated metadata (.txt) with columns gene_aligned and sequence_length
#################################################################

# Loading libraries
rm(list=ls())
library(dplyr)
library(Biostrings)
library(readr)
library(stringr)
library(optparse)

## Parsing arguments
option_list = list(
  make_option(c("-d", "--meta"), type="character", default=NULL,
              help="metadata file path", metavar="character"),
  make_option(c("-a", "--aln"), type="character", default=NULL,
              help="sequences file path (FASTA)", metavar="character"),
  make_option(c("-r", "--ref_seq"), type="character", default=NULL,
              help="reference sequence file path (FASTA)", metavar="character"),
  make_option(c("-o", "--meta_out"), type="character", default=NULL,
              help="metadata output file path", metavar="character"),
  make_option(c("-g", "--annot"), type="character", default=NULL,
              help="gene annotation file path (e.g., TSV or CSV)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Verify required packages are installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("The 'Biostrings' package is required but not installed. Please install it with: install.packages('BiocManager') and then BiocManager::install('Biostrings')", call. = FALSE)
}

# Read files
meta <- read.csv(opt$meta)
user_sequences <- readBStringSet(opt$aln)
ref_sequence <- readBStringSet(opt$ref_seq)
gene_annotations <- read.csv(opt$annot)

# Create new columns for the results
meta$gene_aligned <- NA
meta$sequence_length <- NA

# Create a list to store alignments for each gene
all_gene_alignments <- list()

# Loop through each user sequence
for (i in 1:length(user_sequences)) {
  seq_name <- names(user_sequences)[i]
  
  # Use sequence complete name as ID
  accession_id <- seq_name
  
  # Align the user sequence to the reference
  alignment <- pairwiseAlignment(user_sequences[i], ref_sequence, type = "global")
  aligned_start <- start(subject(alignment))
  aligned_end <- end(subject(alignment))
  
  # Find all genes with overlap
  matching_genes <- c()
  
  for (j in 1:nrow(gene_annotations)) {
    gene_name <- gene_annotations[j, "gene_name"]
    gene_start <- gene_annotations[j, "start"]
    gene_end <- gene_annotations[j, "end"]
    
    # Calculate overlap between the user sequence alignment and the gene annotation
    overlap_start <- max(aligned_start, gene_start)
    overlap_end <- min(aligned_end, gene_end)
    current_overlap <- max(0, overlap_end - overlap_start + 1)
    
    # If there is any overlap, process the gene
    if (current_overlap > 0) {
      matching_genes <- c(matching_genes, gene_name)
      
      # The start and end positions of the gene, relative to the aligned subject sequence
      ref_aln_start <- start(subject(alignment))
      
      # Calculate the start and end positions of the overlapping region within the aligned pattern sequence
      sub_aln_start <- overlap_start - ref_aln_start + 1
      sub_aln_end <- overlap_end - ref_aln_start + 1
      
      # Slice the aligned pattern (user) sequence
      sub_seq_in_gene <- subseq(aligned(pattern(alignment)), start = sub_aln_start, end = sub_aln_end)
      
      # Create a new BStringSet with the aligned user sequence
      sub_aln_output <- BStringSet(toString(sub_seq_in_gene))
      
      # Create proper names for the output sequences
      names(sub_aln_output) <- accession_id
      
      # Add the alignment to the list, grouped by gene name
      if (is.null(all_gene_alignments[[gene_name]])) {
        all_gene_alignments[[gene_name]] <- sub_aln_output
      } else {
        all_gene_alignments[[gene_name]] <- c(all_gene_alignments[[gene_name]], sub_aln_output)
      }
    }
  }
  
  # If no gene was found, set the value to "unknown"
  if (length(matching_genes) == 0) {
    best_match_gene <- "unknown"
  } else {
    best_match_gene <- paste(unique(matching_genes), collapse = ", ")
  }
  
  # Add the results to the metadata table
  meta_row <- which(meta$Accession == accession_id)
  if (length(meta_row) > 0) {
    meta$gene_aligned[meta_row] <- best_match_gene
    meta$sequence_length[meta_row] <- nchar(user_sequences[i])
  }
}

# Write the final processed metadata file
write_delim(meta, opt$meta_out, delim = '\t', quote = 'none')

# Write all gene alignments to separate FASTA files
for (gene_name in names(all_gene_alignments)) {
  
  # Obtain gene posittions from the annotation
  gene_info <- gene_annotations[gene_annotations$gene_name == gene_name, ]
  gene_start <- gene_info$start
  gene_end <- gene_info$end
  
  # Extract the reference sequence e convert it to a BStringSet
  gene_ref_seq <- BStringSet(subseq(ref_sequence[[1]], start = gene_start, end = gene_end))
  names(gene_ref_seq) <- paste0("REF_", gene_name)
  
  # Add reference sequence (only once) to the list of sequences of that gene
  final_alignments <- c(gene_ref_seq, all_gene_alignments[[gene_name]])
  
  # Save fasta file
  output_file_name <- paste0(gene_name, "_alignments.fasta")
  writeXStringSet(final_alignments, file = output_file_name)
}

# Message to inform the user that the script has finished
print("Script finished successfully. All files have been processed.")
