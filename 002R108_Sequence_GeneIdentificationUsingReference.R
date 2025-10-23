#################################################################
## GENE IDENTIFICATION USING A REFERENCE SEQUENCE
## Type: Sequence
## Author: Dennis Maletich Junqueira, Tiago Azevedo Pereira
## AI Disclosure: Gemini generated the first draft
## Date: 2025-08-27
#################################################################
## DESCRIPTION:
## This script aligns user-provided sequences to a reference genome of moderate length,
## identifies the corresponding gene for each sequence, and adds this
## information to a metadata file. It also generates specific
## sub-alignments for each gene found (output: *_alignments.fasta). Additionally, 
## it can eliminate sequences shorter than a specified lenght.
## This script needs Biostrings 3.15 and R 4.2.x.
## 
## RUNNING IN TERMINAL:
## Command: Rscript 002R108_Sequence_GeneIdentificationUsingReference.R --meta my_metadata.csv --aln my_seqs.fasta --ref_seq reference_genome.fasta --annot aanotation.csv --len_min 0 --meta_out annotated_metadata.txt
## --meta = metadata file (.csv). It should contain column "Accession" with the accession ID that is part of the sequence name (example: EPI.ISL.1827381);
## --aln = sequence alignment (.fasta). Sequence names should start with "Accession_" (example: EPI.ISL.1827381_BRAZIL.AC_2023.084). The "Accession" in the metadata file is the extraction of what comes before the "_";
## --ref_seq = reference sequence (.fasta) with all the sequences named following the pattern "Accession_another_information";
## --annot = reference sequence annotation (.tsv). It should contain columns: gene_name,start,end;
## --len_min = minimum sequence length relative to the gene reference (e.75);
## --meta_out = newly annotated metadata (.txt) with columns gene_aligned and sequence_length.
#################################################################
# Loading libraries
rm(list=ls())
library(dplyr)
library(Biostrings)
library(readr)
library(stringr)
library(optparse)

# Start total timer
start_total <- Sys.time()

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
              help="gene annotation file path (e.g., TSV or CSV)", metavar="character"),
  make_option(c("-l", "--len_min"), type="double", default=0,
              help="minimum sequence length relative to the reference gene length", metavar="double")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Verify required packages
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("The 'Biostrings' package is required. Install via BiocManager::install('Biostrings')")
}

# Read files
cat("[INFO]", Sys.time(), "- Reading input files...\n")
meta <- read_delim(opt$meta, delim = ";", quote = "")
user_sequences <- readBStringSet(opt$aln)
ref_sequence <- readBStringSet(opt$ref_seq)
gene_annotations <- read_delim(opt$annot, delim = "\t", quote = "")

# Initialize metadata columns
meta$gene_aligned <- NA
meta$sequence_length <- NA
meta$filtered_status <- "Removed"

# Alignment storage
all_gene_alignments <- list()

cat("[INFO]", Sys.time(), "- Starting sequence processing. Total sequences:", length(user_sequences), "\n")

# Loop through sequences with progress logging
for (i in seq_along(user_sequences)) {
  seq_name <- names(user_sequences)[i]
  
  # Calculate elapsed time and progress
  elapsed <- round(difftime(Sys.time(), start_total, units = "secs"), 2)
  progress_percent <- round((i / length(user_sequences)) * 100, 2)
  cat("[PROGRESS]", Sys.time(), "- Sequence", i, "/", length(user_sequences),
      "(", progress_percent, "% ) processed. Elapsed time:", elapsed, "seconds\n")
  
  # Extract accession
  accession_id <- sub("\\|.*", "", seq_name)  # safe regex
  
  # Align sequence to reference
  alignment <- pairwiseAlignment(user_sequences[i], ref_sequence, type = "global")
  aligned_start <- start(subject(alignment))
  aligned_end <- end(subject(alignment))
  
  matching_genes <- c()
  
  for (j in 1:nrow(gene_annotations)) {
    gene_name <- gene_annotations$gene_name[j]
    gene_start <- gene_annotations$start[j]
    gene_end <- gene_annotations$end[j]
    
    # Overlap calculation
    overlap_start <- max(aligned_start, gene_start)
    overlap_end <- min(aligned_end, gene_end)
    current_overlap <- max(0, overlap_end - overlap_start + 1)
    
    if (current_overlap > 0) {
      matching_genes <- c(matching_genes, gene_name)
      
      sub_aln_start <- overlap_start - aligned_start + 1
      sub_aln_end <- overlap_end - aligned_start + 1
      sub_seq_in_gene <- subseq(aligned(pattern(alignment)), start = sub_aln_start, end = sub_aln_end)
      sub_aln_output <- BStringSet(toString(sub_seq_in_gene))
      names(sub_aln_output) <- accession_id
      
      if (is.null(all_gene_alignments[[gene_name]])) {
        all_gene_alignments[[gene_name]] <- sub_aln_output
      } else {
        all_gene_alignments[[gene_name]] <- c(all_gene_alignments[[gene_name]], sub_aln_output)
      }
    }
  }
  
  best_match_gene <- if(length(matching_genes)==0) "unknown" else paste(unique(matching_genes), collapse=", ")
  
  # Update metadata
  meta_row <- which(meta$Accession == accession_id)
  if (length(meta_row) > 0) {
    meta$gene_aligned[meta_row] <- best_match_gene
    meta$sequence_length[meta_row] <- nchar(user_sequences[i])
  }
}

cat("[INFO]", Sys.time(), "- Finished processing all sequences. Elapsed time:", round(difftime(Sys.time(), start_total, units="secs"),2), "seconds\n")

# Write initial metadata
write_delim(meta, opt$meta_out, delim = '\t', quote = 'none')

# Process alignments per gene
cat("[INFO]", Sys.time(), "- Starting per-gene filtering and saving alignments...\n")
for (gene_name in names(all_gene_alignments)) {
  gene_info <- gene_annotations[gene_annotations$gene_name == gene_name, ]
  gene_start <- gene_info$start[1]
  gene_end <- gene_info$end[1]
  gene_length <- gene_end - gene_start + 1
  gene_ref_seq <- BStringSet(subseq(ref_sequence[[1]], start=gene_start, end=gene_end))
  names(gene_ref_seq) <- paste0("REF_", gene_name)
  
  filtered_alignments <- BStringSet()
  for (i in seq_along(all_gene_alignments[[gene_name]])) {
    seq_to_filter <- all_gene_alignments[[gene_name]][i]
    accession_id <- names(seq_to_filter)
    nucleotide_length <- sum(letterFrequency(seq_to_filter, "ATCG"))
    length_ratio <- nucleotide_length / gene_length
    
    if (nucleotide_length > 0 && length_ratio >= opt$len_min) {
      filtered_alignments <- c(filtered_alignments, seq_to_filter)
      meta_row <- which(meta$Accession == accession_id)
      if (length(meta_row) > 0) meta$filtered_status[meta_row] <- "Kept"
    }
  }
  
  final_alignments <- c(gene_ref_seq, filtered_alignments)
  output_file_name <- paste0(gene_name, "_alignments.fasta")
  writeXStringSet(final_alignments, file = output_file_name)
  cat("[INFO]", Sys.time(), "- Saved alignments for gene", gene_name, "(", length(final_alignments), "sequences )\n")
}

# Write final metadata
write_delim(meta, opt$meta_out, delim = '\t', quote = 'none')
cat("[INFO]", Sys.time(), "- Script finished successfully. Total elapsed time:", round(difftime(Sys.time(), start_total, units="secs"),2), "seconds\n")
