###################################################
## NCBI VIRUS SEQUENCE AND METADATA RETRIEVAL
## Author: Gabriel Teixeira de Macedo
## AI Disclosure: ChatGPT generated the first draft
## Date: 2025-09-01
##
## Description:
## Fetch NCBI Virus sequences using flexible queries.
## Supports taxid, virus name, serotype, host, country,
## collection dates, molecule type, etc.
##
## USAGE: R
## Search you sequences through the function (example):
## main_fetch(
##	  virus_name = NULL,
##	  taxid      = 11320,
##	  serotype   = "H5N1",
##	  host       = NULL,
##	  country    = "Brazil",
##	  date_start = NULL,
##	  date_end   = NULL,
##	  molecule_type = "genomic",
##	  max_results = 2000,
##	  output_dir = "./H5N1_output",
##	  fasta_name_fields = c("accession","organism","gene_segment"),
##	  output_format = c("xlsx","csv","fasta")
##	)
## Outputs: FASTA per segment, Excel, CSV with customizable headers.
###################################################
# Setting working directory
setwd("../Desktop/")

# Auto-install packages
packages <- c("rentrez","dplyr","stringr","readr","progress","xml2","writexl")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies=TRUE)
    library(pkg, character.only=TRUE)
  }
}

# Utility functions
`%||%` <- function(a,b) if(!is.null(a)) a else b
sanitize_meta <- function(x) { x <- gsub("[\r\n\t]", "_", x); x <- gsub("[,;]", "_", x); gsub("[ ]+", "_", x) }
sanitize <- function(x) gsub("[^A-Za-z0-9]", "_", x)

# Build Entrez query
build_query <- function(virus_name=NULL, taxid=NULL, serotype=NULL, host=NULL,
                        country=NULL, date_start=NULL, date_end=NULL, molecule_type="genomic") {
  query_parts <- c()
  if(!is.null(taxid)) {
    query_parts <- c(query_parts, paste0("txid",taxid,"[Organism:exp]"))
  } else if(!is.null(virus_name)) {
    query_parts <- c(query_parts, paste0(virus_name,"[Organism]"))
  } else stop("Provide either virus_name or taxid!")
  
  if(!is.null(serotype) && serotype!="NA") query_parts <- c(query_parts, paste0(serotype,"[All Fields]"))
  if(!is.null(host)) query_parts <- c(query_parts, paste0(host,"[Host]"))
  if(!is.null(country)) query_parts <- c(query_parts, paste0(country,"[Country]"))
  if(!is.null(date_start) & !is.null(date_end)) query_parts <- c(query_parts,
      paste0(date_start,":",date_end,"[Publication Date]"))
  if(!is.null(molecule_type)) {
    mol_prop <- switch(molecule_type, genomic="biomol_genomic[PROP]",
                       protein="biomol_protein[PROP]", mRNA="biomol_mRNA[PROP]", NULL)
    if(!is.null(mol_prop)) query_parts <- c(query_parts, mol_prop)
  }
  
  search_term <- paste(query_parts, collapse=" AND ")
  cat("NCBI Entrez query:\n", search_term, "\n\n")
  return(search_term)
}

# Search NCBI
search_ncbi <- function(term, db="nuccore", retmax=500) {
  res <- rentrez::entrez_search(db=db, term=term, retmax=retmax, use_history=TRUE)
  if(res$count==0) stop("No sequences found!")
  cat("✅ Found", res$count,"entries.\n")
  return(res)
}

# Fetch metadata
fetch_metadata <- function(web_history, db="nuccore", batch_size=200) {
  total <- web_history$count
  pb <- progress::progress_bar$new(total=ceiling(total/batch_size), format="Downloading metadata [:bar] :percent :eta")
  metadata_list <- list()
  
  for(start in seq(0,total-1,batch_size)){
    pb$tick()
    ids <- rentrez::entrez_fetch(db=db, web_history=web_history$web_history, rettype="gb", retmode="xml", retstart=start, retmax=batch_size)
    xml_doc <- xml2::read_xml(ids)
    records <- xml2::xml_find_all(xml_doc, "//GBSeq")
    for(rec in records){
      accession <- xml2::xml_text(xml2::xml_find_first(rec,"./GBSeq_accession-version")) %||% NA
      title     <- xml2::xml_text(xml2::xml_find_first(rec,"./GBSeq_definition")) %||% NA
      organism  <- xml2::xml_text(xml2::xml_find_first(rec,"./GBSeq_organism")) %||% NA
      length    <- xml2::xml_text(xml2::xml_find_first(rec,"./GBSeq_length")) %||% NA
      molecule  <- xml2::xml_text(xml2::xml_find_first(rec,"./GBSeq_moltype")) %||% NA
      host <- isolation_source <- location <- collection_date <- isolate <- gene_segment <- NA
      qualifiers <- xml2::xml_find_all(rec,".//GBQualifier")
      for(q in qualifiers){
        name <- xml2::xml_text(xml2::xml_find_first(q,"./GBQualifier_name"))
        value <- xml2::xml_text(xml2::xml_find_first(q,"./GBQualifier_value"))
        if(!is.na(value)){
          if(name=="host") host <- value
          if(name=="isolation_source") isolation_source <- value
          if(name=="collection_date") collection_date <- value
          if(name=="segment") gene_segment <- value
          if(name=="gene") gene_segment <- gene_segment %||% value
          if(name=="isolate") isolate <- value
        }
      }
      match <- stringr::str_match(title,"^[^(]+\\(([^)]+)\\)")
      if(!is.na(match[2])){
        parts <- unlist(strsplit(match[2],"/"))
        if(length(parts)>=3) location <- parts[length(parts)-1]
      }
      metadata_list[[length(metadata_list)+1]] <- data.frame(
        accession=accession, title=title, organism=organism, host=host,
        isolation_source=isolation_source, location=location, collection_date=collection_date,
        length=length, molecule=molecule, gene_segment=gene_segment, isolate=isolate,
        stringsAsFactors=FALSE
      )
    }
  }
  dplyr::bind_rows(metadata_list)
}

# Fetch FASTA per segment
fetch_fasta_by_segment <- function(web_history, metadata_df, batch_size=500, out_dir, fasta_name_fields) {
  total <- web_history$count
  pb <- progress::progress_bar$new(total=ceiling(total/batch_size), format="Downloading FASTA [:bar] :percent :eta")
  fasta_all <- ""
  for(start in seq(0,total-1,batch_size)){
    pb$tick()
    fasta <- rentrez::entrez_fetch(db="nuccore", web_history=web_history$web_history,
                                   rettype="fasta", retstart=start, retmax=batch_size)
    fasta_all <- paste0(fasta_all, fasta)
  }
  
  fasta_lines <- strsplit(fasta_all, "\n")[[1]]
  is_header <- grepl("^>", fasta_lines)
  seq_indices <- which(is_header)
  seq_indices <- c(seq_indices, length(fasta_lines)+1)
  
  seq_list <- vector("list", length(seq_indices)-1)
  headers <- fasta_lines[is_header]
  for(i in seq_along(seq_list)){
    start <- seq_indices[i]+1
    end   <- seq_indices[i+1]-1
    seq_list[[i]] <- paste(fasta_lines[start:end], collapse="")
  }
  
  metadata_df <- metadata_df %>% dplyr::mutate(dplyr::across(where(is.character), sanitize_meta))
  
  # Loop over unique segments
  unique_segments <- unique(metadata_df$gene_segment)
  for(seg in unique_segments){
    seg_indices <- which(metadata_df$gene_segment == seg)
    seg_headers <- sapply(seg_indices, function(i){
      row <- metadata_df[i,]
      paste(sapply(fasta_name_fields, function(f) row[[f]] %||% "NA"), collapse="_")
    })
    seg_seqs <- seq_list[seg_indices]
    fasta_lines_seg <- c(rbind(paste0(">", seg_headers), unlist(seg_seqs)))
    
    # Write segment-specific FASTA
    out_path <- file.path(out_dir, paste0("virus_fasta_segment_", sanitize(seg), ".fasta"))
    writeLines(fasta_lines_seg, out_path)
    cat("FASTA for segment", seg, "saved at:", out_path, "\n")
  }
}

# Main function
main_fetch <- function(virus_name=NULL, taxid=NULL, serotype=NULL, host=NULL,
                       country=NULL, date_start=NULL, date_end=NULL, molecule_type="genomic",
                       max_results=500, output_dir="./", fasta_name_fields=c("accession","organism","gene_segment"),
                       output_format=c("xlsx","csv","fasta")) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  search_term <- build_query(virus_name, taxid, serotype, host, country, date_start, date_end, molecule_type)
  
  tryCatch({
    search_results <- search_ncbi(search_term, retmax=max_results)
    
    cat("\nDownloading metadata...\n")
    metadata_df <- fetch_metadata(search_results)
    
    id_for_file <- if(!is.null(taxid)) taxid else virus_name
    meta_xlsx_out <- file.path(output_dir, paste0("virus_metadata_", sanitize(id_for_file), "_",
                                                 sanitize(serotype), "_", sanitize(host), "_", sanitize(country), ".xlsx"))
    meta_csv_out  <- file.path(output_dir, paste0("virus_metadata_", sanitize(id_for_file), "_",
                                                 sanitize(serotype), "_", sanitize(host), "_", sanitize(country), ".csv"))
    
    if("xlsx" %in% output_format) writexl::write_xlsx(metadata_df, meta_xlsx_out)
    if("csv" %in% output_format) readr::write_csv(metadata_df, meta_csv_out)
    
    cat("\nDownloading FASTA per segment...\n")
    if("fasta" %in% output_format) fetch_fasta_by_segment(search_results, metadata_df, out_dir=output_dir, fasta_name_fields=fasta_name_fields)
    
    cat("\n✅ Completed successfully\n")
    if("xlsx" %in% output_format) cat("Metadata XLSX:", meta_xlsx_out,"\n")
    if("csv" %in% output_format)  cat("Metadata CSV:", meta_csv_out,"\n")
    
  }, error=function(e){
    cat("\n❌ Error:\n", e$message,"\n")
  })
}

# Fetching sequences
# 
main_fetch(
	  virus_name = NULL,
	  taxid      = NULL,
	  serotype   = NULL,
	  host       = NULL,
	  country    = NULL,
	  date_start = NULL,
	  date_end   = NULL,
	  molecule_type = NULL,
	  max_results = NULL,
	  output_dir = "./output",
	  fasta_name_fields = c("accession"),
	  output_format = c("xlsx","csv","fasta")
	)