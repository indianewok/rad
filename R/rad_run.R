rad_run<-function(read_layout_form, read_dir, output_dir, whitelist_path, nthreads = 1,
  GENERATE_WHITELIST_ONLY = FALSE, SETUP_DRY_RUN = FALSE, FULL_RUN = FALSE, 
  VERBOSE = FALSE, KINDA_VERBOSE = FALSE, 
  ORIGINAL_CHUNK_SIZE = 5e5){
  if(VERBOSE||KINDA_VERBOSE){
    print("Welcome to RAD! Now testing that everything's in the right place...")
    print(paste0("The files are coming in from ", read_dir, " and the files are going to be written to ", output_dir, "..."))
  }
  #DIRECTORY SETUP & MEMORY
  if(is.na(output_dir)){
    print("No output directory provided! Please provide one and try again!")
    stop()
  }
  if(!dir.exists(output_dir)){
    if(VERBOSE){
      print("Output directory doesn't exist (yet!) Making it now.")
    }
    dir.create(path = output_dir)
  } else {
    if(VERBOSE||KINDA_VERBOSE){
      print("Directory exists! Checking subdirectories now...")
    }
  }
    subdirs<-c("/barcode_info","/sigstrings","/output")
    for(i in 1:length(subdirs)){
      if(!dir.exists(paste0(output_dir, subdirs[i]))){
        dir.create(path = paste0(output_dir, subdirs[i]))
        if(VERBOSE){
          print(paste0("Generated .", subdirs[i], ": ", dir.exists(paste0(output_dir, subdirs[i]))))
        }
      } else {
        if(VERBOSE){
          print("This file already exists!")
        }
        next
      }
    }
    subdirs<-paste0(output_dir, subdirs)
    if(VERBOSE){
      print("Done with directory setup! Everything looks good on this end. Now checking the files coming in...")
      print(paste0("Directory with your reads in it: ", read_dir))
      print(subdirs)
    }
    if(VERBOSE){
      ifelse(test = file.exists(read_dir),yes = {
        print("Directory found! Things are going swimmingly...")
      },no = {
        print("Couldn't find the read directory, check if the path is correct...")
        stop()
      })
    }
  files<-list.files(path = read_dir, pattern = "\\d.fastq.gz", full.names = TRUE)
  file_size<-memuse::Sys.filesize(files[1])
  if(VERBOSE){
    print(paste0("Found ", length(files), " in the directory with the extension .fastq.gz. Each one is about ", file_size, " large."))
    print("Now importing one to see if we're all good...")
    print(paste0("File found: ", file.exists(files[1])))
  }
  if(file.exists(files[1])){
    test_file<-read_fastqas(fn = files[1], type = "fq")
  } else {
    print("Couldn't find the files inside the directory, something's funky!")
    stop()
  }
  if(nrow(test_file) != 0){
    file_rows<-nrow(test_file)
    
    if(VERBOSE){
      print(paste0("This .fastq has ", nrow(test_file), " reads in it and will be ", memuse::memuse(test_file), " in total."))
      size<-ceiling(ORIGINAL_CHUNK_SIZE/(nrow(test_file)))
      print(paste0("That'd be about ", size,  " .fastq files and will take about ",
        memuse::memuse(test_file)*size, "."))
    }
  } else {
    if(VERBOSE){
      print("Looks like there's a file error! Something's up with the .fastq files.")
    }
    stop()
  }
  if(SETUP_DRY_RUN != FALSE){
    if(VERBOSE||KINDA_VERBOSE){
      print("Stopping here! Things that work on your computer: read_fastqas, file calculations, and making directories!")
    }
    return(list2env(list(test_file = test_file), .GlobalEnv))
  } else {
    print("Now importing...")
  }
  setup_chunk_size<-ceiling(ORIGINAL_CHUNK_SIZE/nrow(test_file))
  df<-read_fastqas(fn = files[1:setup_chunk_size], type = "fq")
  if(VERBOSE||KINDA_VERBOSE){
    print(nrow(df))
    if(VERBOSE){
      print(paste0("Total size of imported chunk in memory: ", memuse::memuse(df)))
      total_mem<-memuse::Sys.meminfo()[1]
      print(total_mem)
      print(paste0("Now that import looks good, let's check our read layout and our misalignment thresholds!"))
      #MISALIGNMENT THRESHOLD & READ LAYOUT
      print(paste0("Read layout form is at ", read_layout_form))
    }
  }
  prep_read_layout(read_layout_form)
  if(VERBOSE){
    print("Read Layout: ")
    print(read_layout)
    print("Static sequences: ")
    print(adapters)
    print("Now calculating misalignment threshold...")
  }
  baseline_filter <- sum(na.omit(read_layout[direction == "forward", expected_length])) + 100
  df <- df[stringr::str_length(seq) >= baseline_filter]
  strl<-ggpubr::ggdensity(stringr::str_length(df$seq), xscale = "log10")
  ggplot2::ggsave(filename = paste0(output_dir,"/read_length.png"), strl)
  rm(strl)
  misalignment_threshold<-sigalign_stats(adapters = adapters[-grep("poly", x = names(adapters))], 
    sequences = df$seq[1:100000], nthreads = 1) %>% stat_collector(., read_layout, mode = "stats")
  data.table::fwrite(misalignment_threshold, file = paste0(output_dir, "/misalignment_threshold.csv"))
  print(misalignment_threshold)
  #INSERT MISALIGNMENT THRESHOLD WRITE HERE
  if(VERBOSE||KINDA_VERBOSE){
    print("Cool! Got the misalignment thresholds. Now checking if our end-to-end rad_chunk function works on a read subset...")
  }
  if(nthreads > 1){
    if(VERBOSE){
      print("Using multiple cores!")
    }
  }
  df<-rad_chunk(df = df, read_layout = read_layout, misalignment_threshold = misalignment_threshold, 
    nthreads = nthreads, output_dir = output_dir, verbose = VERBOSE)
  if(VERBOSE||KINDA_VERBOSE){
    print("Looks like it works like a charm! Now generating whitelist from the first subset of reads...")
  }
  original_whitelist<-NULL
  if(file.exists(whitelist_path)){
    original_whitelist<-whitelist_importer(whitelist_path = whitelist_path)
  }
  whitelist_generator(df = df, original_whitelist = original_whitelist,
     prefiltered_whitelist = NULL, output_dir = output_dir, stringency = "DEFAULT", verbose = VERBOSE)
  if(GENERATE_WHITELIST_ONLY == TRUE){
    if(!exists("generated_whitelist")||nrow(generated_whitelist) < 5){
      print("Something failed with the whitelist generation!")
      troubleshooting_output<-list(issue_df = df)
      return(list2env(troubleshooting_output))
      stop()
    } else {
      return(print("Whitelist generated!"))
    }
  }
  if(VERBOSE){
    print("Whitelist is generated! Now correcting the df in question...")
  }
  data.table::fwrite(generated_whitelist, file = paste0(output_dir,"/barcode_info/generated_whitelist.csv"))
  barcode_corrector(df, gen_whitelist = generated_whitelist, breadth = 1, depth = 2, high_speed = TRUE,
     nthreads = nthreads, maxDistance = 2, verbose = VERBOSE)
  per_section_pass<-file.exists(paste0(output_dir,"/pass_fail_chunk_stats.csv"))
  pass_fail<-table(df_filtered$pass_fail) %>% data.table::as.data.table(.)
  data.table::fwrite(pass_fail, file = paste0(output_dir, "/pass_fail_chunk_stats.csv"), append = per_section_pass)
  df_filtered<-df_filtered[which(df_filtered$pass_fail == TRUE),]
  df_filtered<-post_processor(df_filtered)
  
  append_fasta<-file.exists(paste0(output_dir, "/output/filtered_reads.fasta.gz"))
  write_fastqas(df = df_filtered, fn = paste0(output_dir, "/output/filtered_reads.fasta.gz"), 
    type = "fa", append = append_fasta, nthreads = nthreads)
  #add dots for internal barcode assignment--fix this part
  aggc()
  if(VERBOSE||KINDA_VERBOSE){
    print("All done with a *first batch*!")
  }
  list2env(x = list(misalignment_threshold = misalignment_threshold, original_whitelist = original_whitelist), envir = .GlobalEnv)
  if(FULL_RUN == TRUE){
    files<-files[-c(1:setup_chunk_size)]
    file_chunk_divisor<-20
    # if(file_rows > 20000){
    #   file_chunk_divisor<-10
    # }
    file_chunks<-split(files, ceiling(seq_along(files) / file_chunk_divisor))
    if(VERBOSE){
      length(file_chunks)
    }
    pbapply::pblapply(file_chunks, FUN = function(x){
      print(length(x))
      single_processor(fastqas = x, nthreads = nthreads, output_dir = output_dir)
      aggc()
    })
  }
  on.exit(expr = aggc())
}