rad_run<-function(read_layout_form, read_dir, output_dir_, whitelist_path, nthreads = 1, SETUP_DRY_RUN = FALSE){
  print("Welcome to RAD! Now testing that everything's in the right place...")
  #DIRECTORY SETUP & MEMORY
  print(paste0("The files are coming in from ", read_dir, " and the files are going to be written to ", output_dir, "..."))
  if(!dir.exists(output_dir)){
    print("Output directory doesn't exist (yet!) Making it now.")
    dir.create(path = output_dir)
  } else {
    print("Directory exists! Checking subdirectories now...")
    subdirs<-c("/barcode_info","/sigstrings","/plotted_data")
    lapply(subdirs, FUN = function(x){
      if(!dir.exists(paste0(output_dir),x)){
        dir.create(path = paste0(output_dir), x)
        print(paste0("Generated .", x, ": ", file.exists(paste0(output_dir), x)))
      } else {
        print("This file already exists!")
        continue
      }
    })
  }
  print("Done with directory setup! Everything looks good on this end. Now checking the files coming in...")
  print(paste0("Directory with your reads in it: ", read_dir))
  ifelse(test = file.exists(read_dir),yes = {
      print("Directory found! Things are going swimmingly...")
      },no = {
      print("Couldn't find the read directory, check if the path is correct...")
      stop()
    })
  files<-list.files(path = read_dir, pattern = ".fastq.gz", full.names = TRUE)
  file_size<-memuse::Sys.filesize(files[1])
  print(paste0("Found ", length(files), " in the directory with the extension .fastq.gz. Each one is about ", file_size, " big."))
  print("Now importing one to see if we're all good...")
  print(paste0("File found: ", file.exists(files[1])))
  if(file.exists(files[1])){
    test_file<-read_fastqas(fn = files[1], type = "fq")
  } else {
    print("Couldn't find the files inside the directory, something's funky!")
    stop()
  }
  if(nrow(test_file) != 0){
    print(paste0("This .fastq has ", nrow(test_file), " reads in it."))
    print(paste0("This .fastq is ", memuse::memuse(test_file), " large."))
    print(paste0("Set-up requires just about a million reads. That'd be about ", ceiling(1e6/nrow(test_file)), " .fastq file imports and will take about ",
      memuse::memuse(test_file)*ceiling(1e6/nrow(test_file))), ".")
  } else {
    print("Looks like there's a file error! Something's up with the .fastq files.")
    stop()
  }
  if(SETUP_DRY_RUN != FALSE){
    print("Stopping here! Things that work on your computer: read_fastqas, file calculations, and making directories!")
    stop()
  } else {
    print("Now importing...")
  }
  setup_chunk_size<-ceiling(1e6/nrow(test_file))
  df<-read_fastqas(fn = files[1:setup_chunk_size], type = "fq")
  print(paste0("Total size of imported chunk in memory: ", memuse::memuse(df)))
  print(paste0("Total free memory left on your computer: ", memuse::Sys.meminfo()[1]))
  print(paste0("Now that import looks good, let's check our read layout and our misalignment thresholds!"))
  #MISALIGNMENT THRESHOLD & READ LAYOUT
  print(paste0("Read layout form is at "), read_layout_form)
  prep_read_layout(read_layout_form)
  print("Read Layout: ")
  print(read_layout)
  print("Static sequences: ")
  print(adapters)
  print("Now calculating misalignment threshold...")
  baseline_filter <- sum(na.omit(read_layout[direction == "forward", expected_length])) + 100
  df <- df[stringr::str_length(seq) >= baseline_filter]
  misalignment_threshold<-sigalign_stats(adapters = adapters, sequences = df[1:100000,], nthreads = 1) %>% stat_collector(., read_layout, mode = "stats")
  print(misalignment_threshold)
  #INSERT MISALIGNMENT THRESHOLD WRITE HERE
  print("Cool! Got that. Now checking if our end-to-end rad_chunk function works on a read subset...")
  if(nthreads > 1){
    print("Also checking to see if our multi-core code works!")
  }
  df<-rad_chunk(df = df, read_layout = read_layout, misalignment_threshold = misalignment_threshold, nthreads = nthreads, output_file = "", verbose = TRUE)
  print("Looks like it works like a charm! Now generating whitelist from the first subset of reads...")
  original_whitelist<-NULL
  if(file.exists(whitelist_path)){
    original_whitelist<-whitelist_importer(whitelist_path = whitelist_path)
  }
  barcode_handler(df = df, original_whitelist = original_whitelist, generate = TRUE, correct = FALSE)
  #IN PROGRESS
}