###############################################
## COUNTING NUCLEOTIDES, GAPS AND AMBIGUITIES
## Author: Rafaela Bianchin Mozzaquatro, Dennis Maletich Junqueira
## AI Disclosure: final refinements of the manuscript were made using chatgpt
## Date: 2025-09-03
## 
## DESCRIPTION:
## This script identifies the number of the IUPAC nucleotide and ambiguities for each
## given sequence in a fasta file. Additionally, it counts the number of "-".
##
## USAGE: R
###############################################
# Setting working directory
setwd("")

# Loading libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

# Choose the fasta files containing the sequences to be analysed
fasta_file <- "file.fasta"

# Choose the output file
output_file <- "CountingNucleotides.csv"

# Loading the fasta file
fasta <- readBStringSet(fasta_file, format = "fasta", use.names = TRUE)

# Detect all unique letters
all_letters <- uniqueLetters(fasta)
cat("Unique letters detected in FASTA:\n")
print(all_letters)

# Separate gaps from letters
letters_only <- setdiff(all_letters, "-")

# Count all detected letters
base_counts <- letterFrequency(fasta, letters = letters_only)
gap_counts  <- letterFrequency(fasta, "-", as.prob = FALSE)
colnames(gap_counts) <- "gap"

# Prepare results
results <- data.frame(Name = names(fasta),
                      gap_counts,
                      base_counts,
                      stringsAsFactors = FALSE)

# Group uppercase/lowercase dynamically
grouped <- results %>%
  mutate(across(all_of(letters_only), as.numeric))

# Use capital letters
letters_upper <- toupper(letters_only)

# Uniting columns
for (letter in unique(letters_upper)) {
  cols <- c(letter, tolower(letter))
  cols <- cols[cols %in% colnames(grouped)]  # só usa se existir
  grouped[[paste0("n", letter)]] <- rowSums(grouped[, cols, drop = FALSE])
}

# Summing columns
count_cols <- grep("^n[A-Z]$", colnames(grouped), value = TRUE)
grouped$TotalNucleotides <- rowSums(grouped[, count_cols, drop = FALSE])

final_result <- grouped %>%
  select(Name, all_of(count_cols), gap, TotalNucleotides)

# Save output
write.csv(final_result, file = output_file, row.names = FALSE)

cat("✅ Finished! Results saved in:", output_file, "\n")
