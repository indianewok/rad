bajrun<-function(path_layout_form, read_layout_form, 
  test_mode = FALSE){
  prep_seq(read_layout_form = read_layout_form, 
    external_path_form = path_layout_form)
  input_path<-path_layout["file","actual_path"]
  output_path<-path_layout["output_dir","actual_path"]
  if(external_sr_bc == TRUE){
    external_bc_path<-path_layout["external_sr_bc","actual_path"]
    external_bcs<-readLines(external_bc_path[[1]])
  }
  print(paste0("The files are coming in from ", input_path, 
    " and the files are going to be written to ", output_path))
  print("Calculating misalignment distributions with first 100K sequences!")
  if(dir.exists(input_path[[1]])){
    print("Input path is a directory!")
    gz_files<-list.files(path = input_path[[1]], pattern = ".gz", full.names = TRUE)
    file_chunks<-split(gz_files, ceiling(seq_along(gz_files) / chunk_divisor))
    lapply(seq_along(file_chunks), function(i){
      print("Testing to find all file chunks...")
      if(file.exists(file_chunks[[i]][1])){
        print(paste0("Chunk ",i," found!"))
      }
    })
  }
  
  df<-lapply(X=file_chunks[[1]][1:20],
    FUN = function(X){read_fastqas(fn =X, type = "fq")}) %>% 
    data.table::rbindlist(.)
  null_distance<-bajalign_stats(adapters[-grep("poly", x = names(adapters))], 
    sequences = df$fastq_files[1:100000],
    nthreads= 1) %>% stat_collector(., read_layout, mode = "stats")
  
  print("Misalignment Thresholds as follows:")
  
  print(null_distance)
  write.csv(null_distance, file = paste0(output_path[[1]],"/misalignment_distances.csv"))
  pbapply::pblapply(seq_along(file_chunks), function(i){
    print(paste0("now on chunk ", i, " of processing!"))
    df<-lapply(X=file_chunks[[i]],FUN = function(X){jaeger::read_fastqas(fn =X, 
      type = "fq")}) %>% data.table::rbindlist(.)
    df<-df[which(str_length(df$fastq_files) >= 150),]
    sigstrings<-vector(length = length(df$fastq_files))
    print(paste0("Generated a vector to hold your sigstrings that's ", 
      length(df$fastq_files)," spaces long!"))
    sigstrings<-bajalign_sigs(adapters, sequences = df$fastq_files,
      null_distance = null_distance, nthreads_sigstrings) 
    print("Done with sigstringing!")
    invisible(bajbatch(null_distance = null_distance, 
      read_layout = read_layout, sigstrings = sigstrings))
    print("Done with bajbatching!")
    df_new<-baj_extract(sigstrings = sigstrings, whitelist_df = whitelist, 
      df = df, verbose = FALSE, barcorrect = barcorrect,
      nthreads = nthreads_sigstract, jaccard_on = jaccard_on)
    print("Done with baj_extracting!")
    barcodes<-df_new$barcode
    print(paste0(length(unique(barcodes)), " barcodes in total!"))
    if(external_sr_bc == TRUE){
      df_true<-df_new[which(df_new$barcode %in% external_bcs),]
      print(paste0(nrow(df_true), " barcodes in the external list!"))
    } else {
      df_true<-df_new$barcode %>%
        barcodes_to_bits %>% 
        {df_new[which(.%in%whitelist$whitelist_bcs),]}
    }
    df_true$id<-DescTools::StrExtractBetween(x = df_true$sig_id, left = "^", right = "\\*")
    df_true<-df_true %>% tidytable::unite("id", c("id", "umi", "barcode"), sep = "_", remove = FALSE) %>% 
      .[,c("id","filtered_read","filtered_qc")] %>% 
      data.table::setnames(c("id","seq","qual"))
    append_fastq<-file.exists(paste0(output_path, "/filtered_reads.fastq.gz"))
    fastqa_writer(df = df_true, fn = paste0(output_path, "/filtered_reads.fastq.gz"), 
      type = "fq", append = append_fastq)
    append_barcode<-file.exists(paste0(output_path,"/all_barcodes.txt"))
    write(x = barcodes, file = paste0(output_path,"/all_barcodes.txt"), 
      append = append_barcode)
    append_sigstrings<-file.exists(paste0(output_path, "sigsummary.txt"))
    write(x = sigstrings, file = paste0(output_path, "/sigsummary.txt"), 
      append = append_sigstrings)
    print("Done with this chunk!")
  })
}