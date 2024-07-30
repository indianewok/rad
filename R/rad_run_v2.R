rad_test<-function(
    fastq_file_or_directory_path, 
    read_layout_path, 
    output_directory_path, 
    compress, 
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
  barcode_files<-list.files(path = paste0(output_directory_path,"/barcodes/"), full.names = TRUE, pattern = "barcode")
  print("Testing barcode whitelist functionality...")
  
  data.table::setkey(read_layout, "class_id")
  barcode_path<-barcode_files[[1]]
  dir_path<-paste0(dirname(barcode_path),"/")
  barcode_id<-basename(tools::file_path_sans_ext(barcode_path))
  
  whitelist_path<-read_layout[barcode_id, whitelist]
  print(paste0("The whitelist path is ", whitelist_path))
  whitelist<-NULL
  whitelist_path<-ifelse(whitelist_path == "10x_3v3", 
                         system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"), 
                         ifelse(whitelist_path == "10x_3v1", 
                         system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"), 
                         whitelist_path))
  print(paste0("The whitelist path is ", whitelist_path))
  
  pbapply::pblapply(X = barcode_files, FUN = function(X){
    process_barcode(barcode_path = X, read_layout = read_layout)
    return(print(paste0("Processed this barcode file: ", X,"!")))
  })
}