.datatable.aware = TRUE

whitelist_importer<-function(whitelist_path, save = NULL, convert = FALSE){
  print("Importing a little bit of sample data...")
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nrows = 10, nThread = 1,
    col.names = "whitelist_bcs")
  if(class(whitelist$whitelist_bcs) == "character"){
    print("Character type detected! Importing the whitelist.")
    whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1, col.names = "whitelist_bcs")
    if(convert != FALSE){
      print("Converting the whitelist into a more efficient format!")
      whitelist$whitelist_bcs<-barcodes_to_bits(whitelist$whitelist_bcs)
    }
    whitelist<-setkey(whitelist, "whitelist_bcs")
    gc()
    if(!is.null(save)){
      if(save == TRUE){
        new_path<-paste0(dirname(whitelist_path),"/",file_path_sans_ext(whitelist_path)
          %>% basename(.),
          "_bitlist.csv.gz")
        print(paste0("Now saving the integer whitelist to ",new_path,"!"))
        data.table::fwrite(x = whitelist, file = new_path, compress = "auto", nThread = 1, showProgress = TRUE, col.names = FALSE)
        gc()
      }
    }
    return(whitelist)
  }
  else{
    if(class(whitelist$whitelist_bcs) == "integer"){
      print("Integer64-type whitelist detected! Now importing...")
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1,
        col.names = "whitelist_bcs", key = "whitelist_bcs", integer64 = "integer64")
      gc()
      return(whitelist)
    }
  }
}
write_fastqas<-function(df, fn, type, append = FALSE, nthreads = 1){
  if(type == "fq"){
    fq_list<-df %>%
      dplyr::mutate(dummy = "+") %>%
      dplyr::select(id, seq, dummy, qual) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fq_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
  }
  if(type == "fa"){
    if(any(grepl("^>", df$id)) == FALSE){
      df$id<-paste0(">",df$id)
      fa_list<-df %>%
        dplyr::select(id, seq) %>%
        as.matrix() %>%
        t() %>%
        as.character() %>%
        tibble::as_tibble()
      data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
    } else {
    fa_list<-df %>%
      dplyr::select(id, seq) %>%
      as.matrix() %>%
      t() %>%
      as.character() %>%
      tibble::as_tibble()
    data.table::fwrite(fa_list, file = fn, compress = "auto", col.names = FALSE, quote = FALSE, append = append, nThread = nthreads)
    }
  }
}
read_fastqas<-function(fn, type, full_id = FALSE, ...){
  rfq_single<-function(fn, type, full_id = FALSE, ...){
    ext<-tools::file_ext(fn)
    sysname<-Sys.info()[["sysname"]]
    if (ext == "gz"){
      if(sysname == "Linux"){
        cat_cmd<-"zcat"
      }else if(sysname =="Darwin"){
        cat_cmd<-"gunzip -c"
      }else{
        stop("Unsupported OS")
      }
    } else {
      cat_cmd<-"cat"
    }
    if (type == "fq"){
      res<-data.table::fread(
        cmd = glue::glue("{cat_cmd} {fn} | paste - - - - | cut -f1,2,4"),
        col.names = c("id", "seq", "qual"),
        sep = "\t", header = FALSE, quote = "", data.table = TRUE, ...
      )
      res<-res[, .(id = data.table::tstrsplit(id, " ", fixed = TRUE)[[1]], qual, seq)]
      if (full_id) {
        res[, full_id := id]
      }
      return(res)
    }
    if(type == "fa"){
      awk_cmd<-'awk \'/^>/{if (NR>1) printf("\\n"); printf("%s\\t",$0); next} {printf("%s",$0);} END {printf("\\n");}\''
      res<-data.table::fread(
        cmd = glue::glue("{cat_cmd} {fn} | {awk_cmd}"),
        col.names = c("id", "seq"),
        sep = "\t", data.table = TRUE, ...
      )
      return(res)
    }
  }
  if (is.list(fn) || length(fn) > 1) {
    out<-pbapply::pblapply(fn, FUN = function(X) {
      rfq_single(fn = X, type = type, full_id = full_id, ...)
    }, ...)
    dt<-data.table::rbindlist(out)
    return(dt)
  } else {
    return(rfq_single(fn = fn, type = type, full_id = full_id, ...))
  }
}

prep_seq<-function(read_layout_form, external_path_form, create_output_dir = TRUE, data_table_nthreads = 1, ...){
  {
    print("Importing read layout and figuring out its order!")
    tidytable::setDTthreads(threads = 1)
    read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
    
    if(any(read_layout$class != "forw_primer") | any(read_layout$class != "rev_primer")){
      which(read_layout$type == "static") %>% {read_layout[min(.),class := "forw_primer"]}
      which(read_layout$type == "static") %>% {read_layout[max(.),class := "rev_primer"]}
    }
    
    read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])
    
    read_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length),
      expected_length := 12]
    read_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
    read_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]
    
    if(any(duplicated(read_layout$id))){
      read_layout$id[which(duplicated(read_layout$id)|duplicated(read_layout$id, fromLast = TRUE))]<-which(duplicated(read_layout$id)|duplicated(read_layout$id,
        fromLast = TRUE)) %>%
        read_layout$id[.] %>% paste0(., "_", seq(length(.)))
    }
    if(any(duplicated(read_layout$class))){
      read_layout$class[which(duplicated(read_layout$class)|duplicated(read_layout$class, fromLast = TRUE))] %>% unique %>%
        sapply(., FUN = function(x){
          out<-which(read_layout$class == x) %>%
            read_layout[., class_id := paste0(class, "_", seq(.))]
          return(out)
        }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
      read_layout$class_id[which(is.na(read_layout$class_id))]<-read_layout$class[which(is.na(read_layout$class_id))]
    } else {
      read_layout$class_id<-read_layout$class
    }
    read_layout<-read_layout %>%
      .[,direction := "forward"] %>%
       data.table::copy(.) %>%
      .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
      .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"),paste0("rc_", id))] %>%
      .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
      .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
      .[, direction := "reverse"] %>%
      .[order(seq(nrow(.)), decreasing = TRUE),] %>%
      #This part does the reversing
      {rbind(read_layout, .)} %>%
      # This code annotates the reverse variable elements with the "rc" in the type column.
      {.[,class_id := ifelse(test = {.$direction == "reverse" &! (.$class_id == "poly_a" | .$class_id == "poly_t")},
        yes = paste0("rc_", .$class_id),
        no = .$class_id)]} %>%
      {.[,"order" := seq(1:nrow(.))]} %>%
      #{.[,"direction" := NULL]} %>%
      {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
        yes = "poly_tail",
        no = .$class)]}
    #adapters<-read_layout$seq[grep(pattern = "adapter", x = read_layout$type)]
    #names(adapters)<-read_layout$id[grep(pattern = "adapter", x = read_layout$type)] #list2env this one
    adapters<-read_layout$seq[-which(is.na(read_layout$seq))]
    names(adapters)<-read_layout$class_id[-which(is.na(read_layout$seq))]
    #super quick fix to make the alignment adapter stuff work bc of hardcoded adapter name assumptions downstream, fix
  }#preparing read layout and whatnot
  {
    print("Checking executable paths!")
    path_layout<-data.table::fread(file = external_path_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1,
      key = "path_type", data.table = TRUE)
    if(dir.exists(path_layout[path_type == "output_dir", actual_path]) == FALSE & create_output_dir == TRUE){
      dir.create(path = path_layout[path_type == "output_dir", actual_path])
    }
    whitelist<-whitelist_importer(whitelist_path = path_layout[path_type == "whitelist", actual_path], ...)
  }#preparing external paths for input
  return(list2env(
    x = list(read_layout = read_layout,
      path_layout = path_layout,
      adapters = adapters,
      whitelist = whitelist), envir = .GlobalEnv))
}
prep_read_layout<-function(read_layout_form){
  print("Importing read layout and figuring out its order!")
  tidytable::setDTthreads(threads = 1)
  read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
  
  if(any(read_layout$class != "forw_primer") | any(read_layout$class != "rev_primer")){
    which(read_layout$type == "static") %>% {read_layout[min(.),class := "forw_primer"]}
    which(read_layout$type == "static") %>% {read_layout[max(.),class := "rev_primer"]}
  }
  
  read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])
  
  read_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length),
    expected_length := 12]
  read_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
  read_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]
  
  if(any(duplicated(read_layout$id))){
    read_layout$id[which(duplicated(read_layout$id)|duplicated(read_layout$id, fromLast = TRUE))]<-which(duplicated(read_layout$id)|duplicated(read_layout$id,
      fromLast = TRUE)) %>%
      read_layout$id[.] %>% paste0(., "_", seq(length(.)))
  }
  if(any(duplicated(read_layout$class))){
    read_layout$class[which(duplicated(read_layout$class)|duplicated(read_layout$class, fromLast = TRUE))] %>% unique %>%
      sapply(., FUN = function(x){
        out<-which(read_layout$class == x) %>%
          read_layout[., class_id := paste0(class, "_", seq(.))]
        return(out)
      }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
    read_layout$class_id[which(is.na(read_layout$class_id))]<-read_layout$class[which(is.na(read_layout$class_id))]
  } else {
    read_layout$class_id<-read_layout$class
  }
  read_layout<-read_layout %>%
    .[,direction := "forward"] %>%
    data.table::copy(.) %>%
    .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
    .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"),paste0("rc_", id))] %>%
    .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
    .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
    .[, direction := "reverse"] %>%
    .[order(seq(nrow(.)), decreasing = TRUE),] %>%
    #This part does the reversing
    {rbind(read_layout, .)} %>%
    # This code annotates the reverse variable elements with the "rc" in the type column.
    {.[,class_id := ifelse(test = {.$direction == "reverse" &! (.$class_id == "poly_a" | .$class_id == "poly_t")},
      yes = paste0("rc_", .$class_id),
      no = .$class_id)]} %>%
    {.[,"order" := seq(1:nrow(.))]} %>%
    #{.[,"direction" := NULL]} %>%
    {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
      yes = "poly_tail",
      no = .$class)]}
  #adapters<-read_layout$seq[grep(pattern = "adapter", x = read_layout$type)]
  #names(adapters)<-read_layout$id[grep(pattern = "adapter", x = read_layout$type)] #list2env this one
  adapters<-read_layout$seq[-which(is.na(read_layout$seq))]
  names(adapters)<-read_layout$class_id[-which(is.na(read_layout$seq))]
  #super quick fix to make the alignment adapter stuff work bc of hardcoded adapter name assumptions downstream, fix
  return(list2env(x = list(read_layout = read_layout, adapters = adapters), envir = .GlobalEnv))
}
prep_external_path<-function(external_path_form){
  print("Checking executable paths!")
  path_layout<-data.table::fread(file = external_path_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1,
    key = "path_type", data.table = TRUE)
  if(dir.exists(path_layout[path_type == "output_dir", actual_path]) == FALSE & create_output_dir == TRUE){
    dir.create(path = path_layout[path_type == "output_dir", actual_path])
  }
  whitelist<-whitelist_importer(whitelist_path = path_layout[path_type == "whitelist", actual_path], ...)
  return(list2env(x = list(path_layout = path_layout, whitelist = whitelist), envir = .GlobalEnv))
}
stat_collector<-function(df, read_layout, mode = "stats"){
  forward_adapters<-read_layout[type %in% c("static") & direction == "forward" & class != "poly_tail", class_id]
  reverse_adapters<-read_layout[type %in% c("static") & direction == "reverse" & class != "poly_tail", class_id]
  df<-tidytable::as_tidytable(df)
  forward_zero_ids<-df %>%
    tidytable::filter(query_id %in% forward_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance == 0)) %>%
    tidytable::pull(id)
  reverse_zero_ids<-df %>%
    tidytable::filter(query_id %in% reverse_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance == 0)) %>%
    tidytable::pull(id)
  if(mode == "stats"){
    stats_reverse<-df %>%
      tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
      tidytable::group_by(query_id) %>%
      tidytable::summarise(
        misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance))
    stats_forward<-df %>%
      tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
      tidytable::group_by(query_id) %>%
      tidytable::summarise(
        misal_threshold = mean(best_edit_distance),
        misal_sd = sd(best_edit_distance))
    misalignment_thresholds<-rbind(stats_forward, stats_reverse)
    return(misalignment_thresholds)
  }
  if(mode == "graphics"){
    stats_reverse<-df %>%
      tidytable::filter(id %in% forward_zero_ids, query_id %in% reverse_adapters) %>%
      tidytable::group_by(query_id)
    stats_forward<-df %>%
      tidytable::filter(id %in% reverse_zero_ids, query_id %in% forward_adapters) %>%
      tidytable::group_by(query_id)
    total_stats<-data.table::rbindlist(list(stats_forward, stats_reverse))
    return(total_stats)
  }
}
read_paf<-function(path, ...){
  df<-data.table::fread(input = path, 
    fill = TRUE, sep = NULL, strip.white = FALSE, header = FALSE )
  df<-df[[1]] %>% stringr::str_split(., pattern = "\t", simplify = TRUE) %>% {data.table::as.data.table(.)}
  df[df == '']<-NA
  ids<-c("id","ql","qs","qe","qr","sn","sl","ss","se","mm","mt","mq")
  data.table::setnames(df, c(ids,
    df[,ncol(df):ncol(df),] %>% {
      which(!is.na(.))} %>% 
      df[.,13:ncol(df),] %>% 
      {substr(unlist(.), start = 1, stop = 2)} %>% 
      unique))
  for (col in (length(ids) + 1):ncol(df)) {
    tag_prefixes<-substr(df[[col]], start = 1, stop = 2)
    expected_tag<-colnames(df)[col]
    mismatch_indices<-which(tag_prefixes != expected_tag & !is.na(tag_prefixes))
    if (length(mismatch_indices) > 0) {
      # Shifting columns to the right for mismatched rows
      df[mismatch_indices, (col + 1):ncol(df) := .SD, .SDcols = col:(ncol(df) - 1)]
      df[mismatch_indices, (col) := NA]  # Set current column to NA for mismatched rows
    }
  }
  for (col in names(df)) {
    numeric_conversion <- suppressWarnings(as.numeric(df[[col]]))
    if (!any(is.na(numeric_conversion)) || all(is.na(numeric_conversion) == is.na(df[[col]]))) {
      df[[col]] <- numeric_conversion
    }
  }
  return(df)
}
rad_chunk<-function(df, read_layout, misalignment_threshold, nthreads, output_dir, verbose = FALSE, write_out = TRUE){
  baseline_filter <- sum(na.omit(read_layout[direction == "forward", expected_length])) + 100
  df <- df[stringr::str_length(seq) >= baseline_filter]
  write_out<-(nchar(output_dir) > 0)
  if(verbose == TRUE){
    print("Check #1--pre-sigalign.")
    print(lobstr::mem_used())
    print(memuse::Sys.procmem(gcFirst = FALSE))
  }
  
  sigstrings<-vector(length = length(df$seq))
  sigstrings<-sigalign(adapters = adapters, sequences = df$seq, ids = df$id, 
    misalignment_threshold = misalignment_threshold, nthreads = nthreads)
  
  if(write_out){
    output_path<-paste0(output_dir, "/sigstrings")
    append_sigstrings<-file.exists(paste0(output_path, "sigsummary_unprocessed.txt"))
    write(x = sigstrings, file = paste0(output_path, "/sigsummary_unprocessed.txt"), 
      append = append_sigstrings)
  }
  
  if(verbose == TRUE){
    print("Check #2--post-alignment.")
    print(lobstr::mem_used())
    print(memuse::Sys.procmem(gcFirst = FALSE))
  }
  
  sigstrings<-sigrun(read_layout, misalignment_threshold, sigstrings = sigstrings, 
    nthreads = nthreads, verbose = FALSE)
  
  if(write_out){
    output_path<-paste0(output_dir, "/sigstrings")
    append_sigstrings<-file.exists(paste0(output_path, "sigsummary_processed.txt"))
    write(x = sigstrings, file = paste0(output_path, "/sigsummary_processed.txt"), 
      append = append_sigstrings)
  }
  
  if(verbose == TRUE){
    print("Check #3--post-processing.")
    print(lobstr::mem_used())
    print(memuse::Sys.procmem(gcFirst = FALSE))
  }
  
  df<-sig_extractor_v3(read_layout, misalignment_threshold = misalignment_threshold, df = df,
     processed_sigstrings = sigstrings[grep(pattern = "undecided", x = sigstrings, invert = TRUE)],
    verbose = FALSE) %>% data.table::as.data.table(.)
  
   if(verbose == TRUE){
     print("Check #4--extracting sigs to generate the df.")
     print(lobstr::mem_used())
     print(memuse::Sys.procmem(gcFirst = FALSE))
   }
  
   aggc()
   rm(sigstrings)
   
   if(verbose == TRUE){
     print("Check #5--clearing processed_sigs.")
     print(lobstr::mem_used())
     print(memuse::Sys.procmem(gcFirst = FALSE))
   }
   
   length_filter<-apply(df, MARGIN = 2, FUN = function(x){which(stringr::str_length(x) <= 2)},
      simplify = TRUE) %>% unlist %>% unique
   if(length(length_filter) > 0){
     df<-df[-length_filter,]
   }
   
  if(verbose == TRUE){
    print("Final check!")
    print(lobstr::mem_used())
    print(memuse::Sys.procmem())
  }
   
   on.exit(expr = aggc(), add = TRUE)
  return(df)
}
subset_module_analysis<-function(df, whitelist){
  df_out$direction<-str_extract(string = df_out$sig_id, pattern = ":.>$")
  df_out$direction<-DescTools::StrExtractBetween(x = df_out$direction, left = ":", right = ">", greedy = FALSE)
  df_out$concatenate<-grepl(pattern = "FR_RF", x= df_out$sig_id)
  barcodes<-table(df_out$barcode, df_out$concatenate) %>% as.data.table
  colnames(barcodes)<-c("bc","concatenate_extracted","freq")
  barcodes$wl_matched<-barcodes$bc %>% sequence_to_bits(.) %>% {.%in%whitelist$whitelist_bcs}
}
synth_data_processor<-function(fn, type, df_out){
  if(type == "gz"){
    cat_cmd<-"zcat"
  } else {
    cat_cmd<-"cat"
  }
  cmd<-glue::glue("{cat_cmd} {fn} | paste - - - - | awk -F'\t' '{{ split($1, a,
     \" \"); print a[1]\",\"a[2] }}' | cut -d',' -f1,2,3,4,5,6,7")
  res<-data.table::fread(cmd = cmd, sep = ",", header = FALSE, quote = "", fill = TRUE)
  df_out$id<-stringr::str_extract(string = df_out$sig_id, pattern = "<.+?>") %>% 
  DescTools::StrExtractBetween(x = ., left = ":", right = ":")
  df_out$concatenate_state<-ifelse(test = grepl(pattern = "+FR_RF",
     x = df_out$sig_id, fixed = TRUE),
    yes = {
      "TRUE"
    },
    no = {
      "FALSE"
    })
  df_out$id<-gsub(pattern = "\\+FR_RF", replacement = "", x =df_out$id)
  df_new<-dplyr::left_join(df_out, res, by = c("id" = "V1"))
  df_new[, c("original_barcode", "original_umi", "original_read") := data.table::tstrsplit(V2, "_", type.convert = TRUE)]
  df_new$original_barcode<-gsub(pattern = "random", replacement = NA, x = df_new$original_barcode)
  df_new$original_umi<-gsub(pattern = "seq", replacement = NA, x = df_new$original_barcode)
  if(length(grep(pattern = "junk", x = df_new$original_barcode)) > 1){
    df_new<-df_new[-grep(pattern = "junk", x = df_new$original_barcode),]
    #the adapter sequences and original read structure mess up jr here
  }
  return(df_new)
}

whitelist_generator<-function(df, original_whitelist = NULL, prefiltered_whitelist = NULL, 
  output_dir = NULL, verbose = FALSE, stringency = "LOW"){
  if(!is.null(original_whitelist)||!is.null(prefiltered_whitelist)){
    barcodes<-table(df$barcode) %>% data.table::as.data.table(.)
  } else {
    barcodes<-table(df$barcode[grep(pattern = "FR_RF", x = df$sig_id, invert = TRUE)]) %>% data.table::as.data.table(.)
  }
  ratio<-nrow(barcodes)/nrow(df)
  if(verbose){
    print(paste0("Number of total barcodes: ", nrow(barcodes)))
    print(whitelist_size())
    print(paste0("Ratio of unique barcodes to unique reads: ", ratio))
  }
  if(ratio > 0.5){
    rescue_repeats<-6
  } else {
    rescue_repeats<-8
  }
  filter_runs<-sapply(X = c("A","C","T","G"), FUN = function(X){
    paste0(replicate(X, n = rescue_repeats), collapse = "") %>% {which(grepl(pattern = .,  x = barcodes$V1) == TRUE)}
  }) %>% unlist %>% unique
  barcodes<-barcodes[-filter_runs,]
  barcodes$pois_dist<-stats::ppois(q = barcodes$N, lambda = mean(barcodes$N))
  barcodes$V1<-sequence_to_bits(barcodes$V1)
  if(verbose){
    print(head(barcodes))
  }
  if(!is.null(original_whitelist)){
    if(verbose){
      print(table(barcodes$V1 %in% original_whitelist$whitelist_bcs))
    }
    chunk_whitelist<-barcodes[barcodes$V1 %in% original_whitelist$whitelist_bcs,]
    if(verbose){
      print(nrow(chunk_whitelist))
    }
  } else if(!is.null(prefiltered_whitelist)){
    if(verbose){
      print(table(barcodes$V1 %in% prefiltered_whitelist$whitelist_bcs))
    }
    chunk_whitelist<-barcodes[barcodes$V1 %in% prefiltered_whitelist$whitelist_bcs,]
  } else
  if(is.null(prefiltered_whitelist) && is.null(original_whitelist)){
    chunk_whitelist<-barcodes[which(barcodes$pois > 0.9),] #first pass
  }
  {
    if(verbose){
      print(head(chunk_whitelist))
    }
    chunk_whitelist<-chunk_whitelist[order(chunk_whitelist$N, decreasing = TRUE),]
    chunk_whitelist$pois_dist<-stats::ppois(q = chunk_whitelist$N, lambda = mean(chunk_whitelist$N))*100
    suppressWarnings(fitcheck<-MASS::fitdistr(chunk_whitelist$N, "negative binomial"))
    chunk_whitelist$ncsum<-pnbinom(size = fitcheck$estimate[1], mu = fitcheck$estimate[2], chunk_whitelist$N)*100
    chunk_whitelist$csum<-(cumsum(chunk_whitelist$N)/sum(chunk_whitelist$N))*100
    if(verbose){
      print(head(chunk_whitelist))
    }
    #hardcoded for the back-two vars
    if(!is.null(output_dir)){
      output<-ggpubr::ggdensity(chunk_whitelist, x = "pois_dist", weight = "N")
      filename<-paste0(output_dir,"/weighted_poisson_distribution_density.png")
      ggplot2::ggsave(filename, output)
    }
    if(stringency != "EXTRA"){
      if(stringency == "LOW"){
        stringency_params<-c(2,2,2,2)
        chunk_whitelist<-whitelist_filterer(df = chunk_whitelist, 
          stringency_params = stringency_params, verbose= TRUE)
        quant_cutoff<-quantile(chunk_whitelist$pois_dist)[1]
        chunk_whitelist<-chunk_whitelist[which(chunk_whitelist$pois_dist >= quant_cutoff),]
      } else
        if(stringency == "DEFAULT"){
          stringency_params<-c(2,2,2,1)
          chunk_whitelist<-whitelist_filterer(df = chunk_whitelist, 
            stringency_params = stringency_params, verbose= TRUE)
          quant_cutoff<-quantile(chunk_whitelist$pois_dist)[2]
          chunk_whitelist<-chunk_whitelist[which(chunk_whitelist$pois_dist >= quant_cutoff),]
        } else 
          if(stringency == "HIGH"){
            stringency_params<-c(2,2,1,1)
            chunk_whitelist<-whitelist_filterer(df = chunk_whitelist, 
              stringency_params = stringency_params, verbose= TRUE)
            quant_cutoff<-quantile(chunk_whitelist$pois_dist)[3]
            chunk_whitelist<-chunk_whitelist[which(chunk_whitelist$pois_dist >= quant_cutoff),]
        }
    }
    if(stringency == "EXTRA"){
      stringency_params<-c(2,2,1,1)
      chunk_whitelist<-whitelist_filterer(df = chunk_whitelist, 
        stringency_params = stringency_params, verbose= TRUE) %>% 
        {whitelist_filterer(df = ., 
        stringency_params = c(2,1,1,1), verbose= TRUE)}
      quant_cutoff<-quantile(chunk_whitelist$pois_dist)[4]
      chunk_whitelist<-chunk_whitelist[which(chunk_whitelist$pois_dist >= quant_cutoff),]
    }
    if(stringency == "NEUROTIC"){
      
    }
  }#stats module
  if(verbose){
    print(head(chunk_whitelist))
    print(nrow(chunk_whitelist))
  }
  output<-list(generated_whitelist = chunk_whitelist, original_whitelist = original_whitelist)
  return(list2env(output, .GlobalEnv))
}
barcode_corrector<-function(df, gen_whitelist, breadth = 1, depth = 2, high_speed = TRUE, nthreads = 1,
   maxDistance = 2, verbose = FALSE){
  populate_whitelist(barcodes = gen_whitelist$V1, poisson_data = gen_whitelist$pois_dist, counts = gen_whitelist$N)
    gen_whitelist<-wl_to_df()
    if(whitelist_size() < 10){
      stop()
    }
    if(!is.null(original_whitelist)){
      barcodes<-table(df$barcode) %>% data.table::as.data.table(.)
    } else {
      barcodes<-table(df$barcode[grep(pattern = "FR_RF", x = df$sig_id, invert = TRUE)]) %>% data.table::as.data.table(.)
    }
    if(verbose){
      print(paste0("Number of total barcodes: ", nrow(barcodes)))
      print(whitelist_size())
    }
    ratio<-nrow(barcodes)/nrow(df)
    if(verbose){
      print(paste0("Ratio of unique barcodes to unique reads: ", ratio))
    }
    if(ratio > 0.5){
      rescue_repeats<-6
      maxDist<-2
    } else {
      rescue_repeats<-8
    }
    filter_runs<-sapply(X = c("A","C","T","G"), FUN = function(X){
      paste0(replicate(X, n = rescue_repeats), collapse = "") %>% {which(grepl(pattern = .,  x = barcodes$V1) == TRUE)}
    }) %>% unlist %>% unique
    barcodes<-barcodes[-filter_runs,]
    if(verbose){
      print(head(gen_whitelist))
    }
    gen_whitelist$sequence<-bit64::as.integer64(gen_whitelist$sequence)
    barcodes$V1<-sequence_to_bits(barcodes$V1)
    query_list<-barcodes[!barcodes$V1 %in% gen_whitelist$sequence,]
    query_list<-query_list[order(query_list$N, decreasing = TRUE),]
    query_list$corrected<-NA
    
    if(verbose){
      print(paste0("The length of the barcodes are: ", nrow(barcodes)))
      print(head(barcodes))
      print(paste0("The length of the gen_whitelist is ", nrow(gen_whitelist)))
      print(head(gen_whitelist))
      print(paste0("the length of the query list is ", nrow(query_list)))
      print(head(query_list))
      print(class(query_list$V1))
    }
    query_list$corrected<-correct_barcodes(
      barcodes = query_list$V1, 
      breadth = breadth, depth = depth,
      high_speed = high_speed, nthreads = nthreads, maxDistance = 4)
    
    df$barcode<-sequence_to_bits(df$barcode)
    df[query_list, barcode := bit64::as.integer64(i.corrected), on = .(barcode = V1)]
    query_list<-query_list[which(!is.na(query_list$corrected)),]
    
    #df<-df[which(df$barcode %in% chunk_whitelist$V1),]
    
    df$pass_fail<-df$barcode %in% gen_whitelist$sequence
    output<-list(df_filtered = df)
    return(list2env(output, .GlobalEnv))
}
post_processor<-function(df){
  df$sig_id<-stringr::str_extract(string = df$sig_id, pattern = "<.+.(?=:)") 
  df$barcode<-df$barcode %>% bit64::as.integer64(.) %>% bits_to_sequence(.)
  df<-df %>% tidytable::unite("id", c("sig_id", "barcode", "umi"), sep = "_", remove = FALSE) %>%
    .[,c("id","read")] %>% data.table::setnames("read","seq")
 return(df)
}
single_processor<-function(fastqas, nthreads, output_dir){
  df<-read_fastqas(fn = fastqas, type = "fq")
  df<-rad_chunk(df, read_layout, misalignment_threshold, nthreads = nthreads, output_dir = output_dir, verbose = FALSE)
  barcode_corrector(df, gen_whitelist = generated_whitelist, breadth = 1, depth = 2, high_speed = TRUE,
    nthreads = nthreads, maxDistance = 2, verbose = FALSE)
  
  pass_fail<-table(df_filtered$pass_fail) %>% data.table::as.data.table(.)
  per_section_pass<-file.exists(paste0(output_dir,"/pass_fail_chunk_stats.csv"))
  
  data.table::fwrite(pass_fail, file = paste0(output_dir, "/pass_fail_chunk_stats.csv"), append = per_section_pass)
  df_filtered<-df_filtered[which(df_filtered$pass_fail == TRUE),]
  df_filtered<-post_processor(df_filtered)
  append_fasta<-file.exists(paste0(output_dir, "/output/filtered_reads.fasta.gz"))
  write_fastqas(df = df_filtered, fn = paste0(output_dir, "/output/filtered_reads.fasta.gz"), 
    type = "fa", append = append_fasta, nthreads = nthreads)
}
whitelist_filterer<-function(df, stringency_params, verbose = FALSE){
  if(verbose){
    print(stringency_params)
  }
  fpoisson<-table(df$pois_dist) %>% data.table::as.data.table(.)
  if(verbose){
    print(head(fpoisson))
  }
  fpoisson<-fpoisson[order(as.numeric(fpoisson$V1), decreasing = FALSE),]
  min_idx<-diff(fpoisson$N, differences = stringency_params[1], lag = stringency_params[2]) %>% 
    sign %>% {diff(.,differences = stringency_params[3], lag = stringency_params[4]) != 0} %>% #toggle lag for param sensitivity?
    which(. == TRUE) %>% 
    min(.)
  df<-fpoisson$V1[min_idx] %>% 
    as.numeric %>% 
    {df[which(df$pois_dist >= .),]}
  if(verbose){
    print(head(df))
  }
    return(df)
}
aggc<-function(){
  gc()
  sysname<-Sys.info()[["sysname"]]
  if(sysname == "Linux"){
    mallinfo::malloc.trim()
  }
}