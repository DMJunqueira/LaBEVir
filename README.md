

![](/LaBEVirLogo.png)



# This is a repository for LaBEVir R scripts:




### 002R108_Sequence_GeneIdentificationUsingReferenceGenome: 
This script aligns user-provided sequences to a reference genome of moderate length, 
identifies the corresponding gene for each sequence, and adds this information to a metadata file. It also generates specific sub-alignments 
for each gene found (output: *_alignments.fasta)

### 002R109_Sequence_GenbankSequenceMetadataRetrieval:
This script programmatically searches the NCBI Nucleotide database for sequences based on a user-provided term. It retrieves full sequence 
data and key metadata (e.g., accession, organism, country, and host) and saves the results into a FASTA file and a CSV file.
