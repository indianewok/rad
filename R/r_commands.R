.datatable.aware = TRUE

whitelist_importer<-function(whitelist_path, save = NULL, convert = FALSE){
  cat("Importing a little bit of your whitelist to figure out its type...\n")
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nrows = 10, nThread = 1)
  if(ncol(whitelist) == 1){
    colnames(whitelist)<-"whitelist_bcs"
  }
  if(class(whitelist$whitelist_bcs) == "character"){
    cat("Character type detected! Importing the whitelist.\n")
    whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1, col.names = "whitelist_bcs")
    if(convert != FALSE){
      cat("Converting the whitelist into a more efficient format!\n")
      whitelist$whitelist_bcs<-sequence_to_bits(whitelist$whitelist_bcs)
    }
    whitelist<-data.table::setkey(whitelist, "whitelist_bcs")
    gc()
    if(!is.null(save)){
      if(save == TRUE){
        new_path<-paste0(dirname(whitelist_path),"/",file_path_sans_ext(whitelist_path)
          %>% basename(.),
          "_bitlist.csv.gz")
        cat(paste0("Now saving the integer whitelist to ",new_path,"!\n"))
        data.table::fwrite(x = whitelist, file = new_path, compress = "auto", nThread = 1, showProgress = TRUE, col.names = FALSE)
        gc()
      }
    }
    return(whitelist)
  }
  else{
    if(class(whitelist$whitelist_bcs) == "integer"){
      cat("Integer64-type whitelist detected! Now importing...\n")
      whitelist<-data.table::fread(input = whitelist_path, header = FALSE, data.table = TRUE, nThread = 1,
        col.names = "whitelist_bcs", key = "whitelist_bcs", integer64 = "integer64")
      gc()
      cat(paste0("Your whitelist is ", nrow(whitelist), " barcodes long!\n"))
      print(head(whitelist))
      return(whitelist)
    }
  }
}

write_fastqas<-function(df, fn, type = "fq", append = FALSE, nthreads = 1){
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

read_fastqas<-function(fn, type = "fq", full_id = FALSE, ...){
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

prep_read_layout<-function(read_layout_form) {
  cat("Importing read layout and figuring out its order!\n")
  tidytable::setDTthreads(threads = 1)
  read_layout<-data.table::fread(file = read_layout_form, header = TRUE, fill = TRUE, na.strings = "", skip = 1)
  read_layout$expected_length[which(!is.na(read_layout$seq))]<-stringr::str_length(read_layout$seq[which(!is.na(read_layout$seq))])
  
  read_layout[(class %in% c("poly_a", "poly_t")) & is.na(expected_length), expected_length := 12]
  read_layout[class == "poly_a", seq := paste0("A{", expected_length, ",}+")]
  read_layout[class == "poly_t", seq := paste0("T{", expected_length, ",}+")]
  
  new_row_top<-data.table::data.table(id = "seq_start", seq = NA, expected_length = 0,
    type = "static", class = "start", direction = "forward", class_id = "seq_start")
  new_row_bottom<-data.table::data.table(id = "seq_stop", seq = NA, expected_length = 0, 
    type = "static", class = "stop", direction = "forward", class_id = "seq_stop")
  
  read_layout<-data.table::rbindlist(list(new_row_top, read_layout, new_row_bottom), fill = TRUE)
  
  if(any(read_layout$class != "forw_primer")){
    start_index<-which(read_layout$class == "start")
    if(read_layout$type[start_index+1] == "static"|is.na(read_layout$class[start_index+1])){
      read_layout[start_index+1, class := "forw_primer"]
    }
  }
  
  if(any(read_layout$class != "rev_primer")){
    stop_index<-which(read_layout$class == "stop")
    if(read_layout$type[stop_index-1] == "static"|is.na(read_layout$class[stop_index-1])){
      read_layout[stop_index-1, class := "rev_primer"]
    }
  }
  
  if (any(duplicated(read_layout$id))) {
    read_layout$id[which(duplicated(read_layout$id) | duplicated(read_layout$id, fromLast = TRUE))]<-
      which(duplicated(read_layout$id) | duplicated(read_layout$id, fromLast = TRUE)) %>%
      read_layout$id[.] %>% paste0(., "_", seq(length(.)))
  }
  if (any(duplicated(read_layout$class))) {
    read_layout$class[which(duplicated(read_layout$class) | duplicated(read_layout$class, fromLast = TRUE))] %>% unique %>%
      sapply(., FUN = function(x){
        out<-which(read_layout$class == x) %>%
          read_layout[., class_id := paste0(class, "_", seq(.))]
        return(out)
      }, simplify = FALSE, USE.NAMES = FALSE) %>% invisible(.)
    read_layout$class_id[which(is.na(read_layout$class_id))]<-read_layout$class[which(is.na(read_layout$class_id))]
  } else {
    read_layout$class_id <- ifelse(test = (read_layout$class == "start" | read_layout$class == "stop"), 
      yes = read_layout$class_id, 
      no = read_layout$class)
  }
  
  read_layout<-read_layout %>%
    .[, direction := "forward"] %>%
    data.table::copy(.) %>%
    .[!class %in% c("start", "stop")] %>%
    .[, seq := ifelse(class == "poly_t", "A{12,}+", ifelse(class == "poly_a", "T{12,}+", sapply(seq, revcomp)))] %>%
    .[, id := ifelse(class %in% c("poly_a", "poly_t"), ifelse(class == "poly_a", "poly_t", "poly_a"), paste0("rc_", id))] %>%
    .[, class_id := ifelse(class_id == "poly_a", "poly_t", ifelse(class_id == "poly_t", "poly_a", class_id))] %>%
    .[, class := ifelse(class == "poly_a", "poly_t", ifelse(class == "poly_t", "poly_a", class))] %>%
    .[, direction := "reverse"] %>%
    .[order(seq(nrow(.)), decreasing = TRUE), ] %>%
    {data.table::rbindlist(list(data.table::data.table(id = "seq_start", seq = NA, expected_length = 0, type = "static",
      class = "start", direction = "reverse", class_id = "seq_start"),
      ., data.table::data.table(id = "seq_stop", seq = NA, expected_length = 0, type = "static",
        class = "stop", direction = "reverse", class_id = "seq_stop")),
      fill = TRUE)} %>%
    {rbind(read_layout, .)} %>%
    {.[, class_id := ifelse(test = {.$direction == "reverse" &!
        (.$class_id == "poly_a" | .$class_id == "poly_t")},# | .$class_id == "seq_start" | .$class_id == "seq_stop")},
      yes = paste0("rc_", .$class_id),
      no = .$class_id)]} %>%
    {.[, "order" := seq(1:nrow(.))]} %>%
    {.[, class := ifelse(test = {.$class == "poly_a" | .$class == "poly_t"},
      yes = "poly_tail",
      no = .$class)]}
  
  adapters <- read_layout$seq[-which(is.na(read_layout$seq))]
  names(adapters) <- read_layout$class_id[-which(is.na(read_layout$seq))]
  return(list2env(x = list(read_layout = read_layout, adapters = adapters), envir = .GlobalEnv))
}

stat_collector<-function(df, read_layout, mode = "stats"){
  forward_adapters<-read_layout[type %in% c("static") & direction == "forward" & class != "poly_tail" & class != "start" & class != "stop", class_id]
  reverse_adapters<-read_layout[type %in% c("static") & direction == "reverse" & class != "poly_tail" & class != "start" & class != "stop", class_id]
  df<-tidytable::as_tidytable(df)
  #changed best_edit_distance from == 0 to <= 1. 
  #not the best fix, as it leaves room for perfect concatenates with no errors or something weird like that
  #but it also adds for some more sequences that might otherwise be discarded.
  forward_zero_ids<-df %>%
    tidytable::filter(query_id %in% forward_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance <= 1)) %>%
    tidytable::pull(id)
  reverse_zero_ids<-df %>%
    tidytable::filter(query_id %in% reverse_adapters) %>%
    tidytable::group_by(id) %>%
    tidytable::filter(all(best_edit_distance <= 1)) %>%
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

aggc<-function(){
  gc()
  sysname<-Sys.info()[["sysname"]]
  if(sysname == "Linux"){
    mallinfo::malloc.trim()
  }
}

process_barcodes<-function(barcode_path, read_layout) {
  barcodes<-data.table::fread(barcode_path, col.names = c("seq", "count"))
  data.table::setkey(read_layout, "class_id")
  
  dir_path<-paste0(dirname(barcode_path),"/")
  barcode_id<-basename(tools::file_path_sans_ext(barcode_path))
  
  expected_length<-read_layout[barcode_id, expected_length]
  barcodes[, filtered := NA_character_]
  barcodes[(stringr::str_length(seq) != expected_length), filtered := "length_filter"]
  rescue_filter<-ceiling(expected_length * 0.5)
  patterns<-sprintf("%s{%d,}", c("A", "C", "T", "G"), rescue_filter + 1)
  poly_runs<-data.table::rbindlist(lapply(patterns, function(pattern) {
    barcodes[is.na(filtered) & grepl(pattern, seq), .(seq, count, filtered = "poly_run_filtered")]
  }))
  if (nrow(poly_runs) > 0){
    barcodes[poly_runs, on = .(seq), filtered := i.filtered]
  }
  static_seqs<-read_layout[type == "static" & !is.na(expected_length) & expected_length > 0 & class != "poly_tail", seq]
  static_qgrams<-colnames(stringdist::qgrams(static_seqs, q = expected_length, useNames = TRUE))
  static_qgrams_mutations<-static_qgrams %>% 
    sequence_to_bits(.) %>% 
    generate_recursive_quaternary_mutations(., mutation_rounds =  2, sequence_length = 16) %>%
    bits_to_sequence(.)
  static_qgrams<-append(static_qgrams, static_qgrams_mutations)
  
  barcodes[is.na(filtered) & seq %in% static_qgrams, filtered := "adapter_segment_filtered"]
  whitelist_path<-read_layout[barcode_id, whitelist]
  cat(paste0("The whitelist path is ", whitelist_path,"\n"))
  whitelist<-NULL
  whitelist_path<-ifelse(whitelist_path == "10x_3v3", 
                           system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"), 
                           ifelse(whitelist_path == "10x_3v1", 
                                  system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"), 
                                  whitelist_path))
  cat(paste0("The whitelist path is ", whitelist_path, "\n"))
  
  if(!is.na(whitelist_path) && file.exists(whitelist_path)){
    ext<-tools::file_ext(whitelist_path)
    convert<-ext == "csv" || ext == "txt"
    whitelist<-whitelist_importer(whitelist_path, convert = convert)
  }
  
  barcodes[, pois_dist := stats::ppois(q = count, lambda = mean(count))]
  barcodes[is.na(filtered), ':='(
    int64_seq = sequence_to_bits(seq),
    int64_rcseq = sequence_to_bits(revcomp(seq))
  )]
  if (!is.null(whitelist)){
    barcodes[(is.na(filtered) & stringr::str_length(seq) == expected_length), filtered := ifelse(
      int64_seq %in% whitelist$whitelist_bcs | int64_rcseq %in% whitelist$whitelist_bcs,
      "whitelist_barcode", "barcode_to_correct"
    )]
    barcodes[(filtered == "whitelist_barcode" & pois_dist >= 0.95), filtered := "pois_validated_barcode"]
  }
  data.table::fwrite(x = barcodes, file = barcode_path)
}

correct_barcodes<-function(input_file,
                           verbose = FALSE,
                           nthreads = 1,
                           mutation_rounds = 3,
                           max_shift = 3,
                           read_layout
                           ){
  barcode<-data.table::fread(
    input_file,
    col.names = c("seq","count","filtered","pois_dist","int64_seq","int64_rcseq"))
  data.table::setkey(read_layout, "class_id")
  data.table::setkey(barcode, "filtered")
  barcode_id<-basename(tools::file_path_sans_ext(input_file))
  
  dir_path<-paste0(dirname(input_file),"/")
  input_fastq_path<-paste0(dirname(dir_path),"/demuxed_reads.fastq.gz")
  filtered_fastq_path<-paste0(dirname(dir_path),"/filtered_reads.fastq.gz")
  unfiltered_fastq_path<-paste0(dirname(dir_path),"/unfiltered_reads.fastq.gz")
  barcode_counts_output_path<-paste0(dir_path, barcode_id, "_counts.csv")
  verbose_barcode_output_path<-paste0(dir_path, barcode_id, "_correction_counts.csv")
  expected_length<-read_layout[barcode_id, expected_length]
  
  generate_and_filter_mutations_v6(
    barcode_header = paste0(barcode_id,":"),
    true_barcodes =  barcode["pois_validated_barcode"]$int64_seq, 
    invalid_barcodes = barcode["barcode_to_correct"]$int64_seq,
    true_counts = barcode["pois_validated_barcode"]$count,
    invalid_counts = barcode["barcode_to_correct"]$count,
    mutation_rounds = mutation_rounds,
    sequence_length = expected_length,
    max_shift = max_shift,
    nthread = nthreads,
    input_fastq = input_fastq_path,
    filtered_fastq = filtered_fastq_path,
    unfiltered_fastq = unfiltered_fastq_path,
    counts_output_csv = barcode_counts_output_path,
    detailed_output_csv = verbose_barcode_output_path,
    verbose_output = TRUE,
    verbose = FALSE)
  
  file.copy(filtered_fastq_path, input_fastq_path, overwrite = TRUE)
}

process_sig<-function(file_path,
                      output_prefix,
                      chunk_size = 50000,
                      max_sequences = -1,
                      compress = TRUE,
                      output_type = "summary",
                      nthreads = 1){
  tabulate_sigs(file_path = file_path, 
                output_prefix = output_prefix, 
                chunk_size = chunk_size,
                max_sequences = max_sequences, 
                compress = compress, 
                nthreads = nthreads,
                output_type = output_type)
  summary_results<-data.table::fread(file = paste0(output_prefix,"/summary_table.csv",  
                                                   ifelse(test = compress, yes = ".gz", no = "")),
                                     header = TRUE)
  valid_reads<-sum(summary_results$unique_id_count[which(summary_results$direction == "F"|summary_results$direction == "R")])
  total_reads<-sum(summary_results$unique_id_count)
  summary_plot<-ggpubr::ggbarplot(data = summary_results, 
                                  x = "direction", 
                                  y = "unique_id_count",
                                  color = "concatenate",
                                  fill = "concatenate", 
                                  title = "SigString Summary", 
                                  legend = "right", 
                                  legend.title = "Concatenate Status", 
                                  xlab = "Read (or Filter) Type",
                                  ylab = "Number of Reads", 
                                  order = c("F", 
                                            "R", 
                                            "missing_read_or_barcode_undecided", 
                                            "invalid_barcode_length_undecided",
                                            "invalid_read_length_undecided"), 
                                  subtitle = paste0(
                                    valid_reads, 
                                    " valid reads out of a total ", 
                                    total_reads, " total reads."), 
                                  x.text.angle = 45, palette = "npj")
  ggplot2::ggsave(filename = paste0(output_prefix,"summary_barplot.pdf"), 
                  plot = summary_plot, 
                  device = "pdf", 
                  width = 8, 
                  height = 11, 
                  units = "in", 
                  dpi = "retina")
  results<-data.table::data.table(unique_id_count = paste0(valid_reads, " valid reads captured from ", total_reads, " total reads."))
  summary_results<-data.table::rbindlist(l = list(summary_results, results), fill = TRUE)
  data.table::fwrite(summary_results, paste0(output_prefix,"/summary_table.csv", ifelse(test = compress, yes = ".gz", no = "")))
}

rad_run<-function(
    fastq_file_or_directory_path, 
    read_layout_path, 
    output_directory_path, 
    compress = TRUE,
    generate_sigstring_diagnostics = FALSE,
    sigstring_diagnostic_verbose = FALSE,
    misalignment_threshold_path = NULL,
    tabulated_sigstring_count = -1,
    chunk_size = 50000,
    mutation_rounds = 3,
    nthreads = 1){
  start_time = Sys.time()
  if(!dir.exists(output_directory_path)){
    dir.create(output_directory_path)
  }
  prep_read_layout(read_layout_form = read_layout_path)
  if(!is.null(misalignment_threshold_path)){
    misalignment_threshold<-data.table::fread(misalignment_threshold_path, header = FALSE)
  } else {
    misalignment_threshold<-misalignment_stream(
      input_path = fastq_file_or_directory_path,
      adapters = adapters,
      nthreads = nthreads)
    misalignment_threshold<-stat_collector(misalignment_threshold, read_layout)
    data.table::fwrite(misalignment_threshold, file = paste0(output_directory_path,"/misalignment_threshold.csv"))
  }
  
  sigstream(input_path = fastq_file_or_directory_path, 
            nthreads = nthreads, 
            compress = compress, 
            write_fastq = TRUE, 
            min_length = 100,
            chunkSize = chunk_size,
            outputPath = paste0(output_directory_path, "/demuxed_reads.fastq"),
            adapters = adapters,
            misalignment_threshold = misalignment_threshold,
            read_layout = read_layout,
            sigstringsFilePath = paste0(output_directory_path, "/sigstrings.txt")
  )
  gc()
  if(generate_sigstring_diagnostics){
    if(sigstring_diagnostic_verbose){
      output_type<-"diagnostic"
    } else {
      output_type<-"summary"
    }
    process_sig(
      file_path = paste0(output_directory_path, "/sigstrings.txt", ifelse(test = compress, yes = ".gz", no = "")),
      output_prefix = output_directory_path,
      max_sequences = tabulated_sigstring_count,
      chunk_size = chunk_size,
      compress = compress,
      nthreads = nthreads,
      output_type = output_type
    )
    gc()
  } else {
    cat("Generate sigstring diagnostics set to FALSE--skipping sigstring tabulation.\n")
  }
  
  dir.create(path = paste0(output_directory_path,"/variable_seqs/"))
  cat("Tabulating barcodes...\n")
  tabulate_barcodes(
    input_files = paste0(output_directory_path, "/demuxed_reads.fastq", 
                         ifelse(test = compress, yes = ".gz", no = "")), 
    output_prefix = paste0(output_directory_path,"/variable_seqs/"), 
    compress = FALSE)
  barcode_files<-list.files(path = paste0(output_directory_path,"/variable_seqs/"), full.names = TRUE, pattern = "barcode")
  out<-lapply(X = barcode_files, FUN = function(X){
    process_barcodes(barcode_path = X, read_layout = read_layout)
    correct_barcodes(input_file = X, read_layout = read_layout, nthreads = nthreads)
  })
  end_time = Sys.time()
  cat(paste0("Finished running! Total run time: ", difftime(end_time, start_time, units='mins'), "\n"))
}

generate_synthetic_data<-function(read_layout_path){
  prep_read_layout(read_layout_path)
  data.table::setkey(read_layout, "class_id")
  
  
  
  cat(paste0("The whitelist path is ", whitelist_path,"\n"))
  whitelist<-NULL
  whitelist_path<-ifelse(whitelist_path == "10x_3v3", 
    system.file(package = "rad", "extdata", "3M-february-2018-3v3.txt_bitlist.csv.gz"), 
    ifelse(whitelist_path == "10x_3v1", 
      system.file(package = "rad", "extdata", "737K-august-2016_bitlist.csv.gz"), 
      whitelist_path))
  cat(paste0("The whitelist path is ", whitelist_path, "\n"))
  if(!is.na(whitelist_path) && file.exists(whitelist_path)){
    ext<-tools::file_ext(whitelist_path)
    convert<-ext == "csv" || ext == "txt"
    whitelist<-whitelist_importer(whitelist_path, convert = convert)
  }
  
}


