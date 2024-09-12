#!/usr/bin/env Rscript

library(docopt)
library(data.table)
library(pbapply)
library(magrittr)


'Usage:
  rad.R --fastq_file_or_directory_path=<fastq_file_or_directory_path> --read_layout_path=<read_layout_path> --output_directory_path=<output_directory_path> --compress=<compress> [--nthreads=<nthreads>]

Options:
  --fastq_file_or_directory_path=<fastq_file_or_directory_path>  Path to the FASTQ file or directory containing multiple files.
  --read_layout_path=<read_layout_path>                           Path to the read layout csv.
  --output_directory_path=<output_directory_path>                 Path to the output directory--if it does not exist, it will be made.
  --compress=<compress>                                           Whether to compress output files with .gz (default is TRUE).
  --nthreads=<nthreads>                                           Number of threads to use (defaults to 1).
'->doc

cmd_rad<-function(){
opts<-docopt(doc)
rad<-function(
    fastq_file_or_directory_path, 
  read_layout_path, 
  output_directory_path, 
  compress = TRUE, 
  nthreads = 1){
  if(!dir.exists(output_directory_path)){
    dir.create(output_directory_path)
  }
  prep_read_layout(read_layout_form = read_layout_path)
  misalignment_threshold<-misalignment_stream(input_path = fastq_file_or_directory_path, adapters = adapters, nthreads = nthreads)
  misalignment_threshold<-stat_collector(misalignment_threshold, read_layout)
  sigstream(input_path = fastq_file_or_directory_path, 
    nthreads = nthreads, 
    compress = compress, 
    write_fastq = TRUE, 
    min_length = 100, 
    chunkSize = 50000,
    outputPath = paste0(output_directory_path, "/demuxed_reads.fastq"),
    adapters = adapters,
    misalignment_threshold = misalignment_threshold,
    read_layout = read_layout,
    sigstringsFilePath = paste0(output_directory_path, "/sigstrings.txt")
  )
  gc()
  tabulate_sigs(file_path = paste0(output_directory_path, "/sigstrings.txt", 
    ifelse(test = compress, yes = ".gz", no = "")),
    output_prefix = output_directory_path)
  gc()
  dir.create(path = paste0(output_directory_path,"/barcodes/"))
  tabulate_barcodes(
    input_files = paste0(output_directory_path, "/demuxed_reads.fastq", 
      ifelse(test = compress, yes = ".gz", no = "")), 
    output_prefix = paste0(output_directory_path,"barcodes/"), 
    compress = FALSE)
  barcode_files<-list.files(path = paste0(output_directory_path,"/barcodes/"), full.names = TRUE)
  pbapply::pblapply(X = barcode_files, FUN = function(X){
    process_barcode(barcode_path = X, read_layout = read_layout)
  })
}

rad(
  fastq_file_or_directory_path = opts$`--fastq_file_or_directory_path`,
  read_layout_path = opts$`--read_layout_path`,
  output_directory_path = opts$`--output_directory_path`,
  compress = as.logical(opts$`--compress`),
  nthreads = as.integer(opts$`--nthreads`)
)
}

if (!interactive() && identical(environment(), globalenv())) {
  cmd_rad()
}
