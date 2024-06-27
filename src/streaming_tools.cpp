#include <Rcpp.h>
#include <zlib.h>
#include "kseq.h"
#include <omp.h>
#include <fstream>
#include "rad.h"

KSEQ_INIT(gzFile, gzread)
  
  struct SearchResult {
    std::string line;
    std::string read_id;
    std::string qual_scores;
  };

std::vector<std::string> sigrun_cpp(
    const std::vector<std::string>& sigstrings, 
    const ReadLayout& read_layout, 
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads = 1, bool verbose = false) {
  VarScan(container, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(container, verbose);
  return process_sigstrings_cpp(container, sigstrings, positionFuncMap, verbose, nthreads);
}

std::vector<SearchResult> readFastqChunk(
    kseq_t* seq, gzFile fp, size_t chunkSize, 
    bool& is_fastq, int& read_count, int max_sequences) {
  std::vector<SearchResult> records;
  for (size_t i = 0; i < chunkSize && kseq_read(seq) >= 0; ++i) {
    if (max_sequences != -1 && read_count >= max_sequences) break;
    SearchResult sr;
    sr.read_id = seq->name.s;
    sr.line = seq->seq.s;
    if (seq->qual.l) {
      sr.qual_scores = seq->qual.s;
      is_fastq = true;
    } else {
      is_fastq = false;
    }
    records.push_back(sr);
    ++read_count;
  }
  return records;
}

Rcpp::CharacterVector sigalign(
    const std::vector<std::string>& queries,
    const std::vector<std::string>& query_names,
    const std::vector<std::string>& sequences,
    const std::vector<std::string>& ids,
    const std::map<std::string, std::pair<double, double>>& null_dist_map,
    int nthreads) {
  
  std::map<int, std::string> signature_map;
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query = queries[j];
      const auto& query_name = query_names[j];
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTail(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        std::set<UniqueAlignment> uniqueAlignments;
        
        char* cquery = const_cast<char*>(query.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csequence, sequence.size(), config);
        std::vector<int> start_positions(cresult.startLocations, cresult.startLocations + cresult.numLocations);
        std::vector<int> end_positions(cresult.endLocations, cresult.endLocations + cresult.numLocations);
        for (int& pos : start_positions) { ++pos; }
        for (int& pos : end_positions) { ++pos; }
        std::sort(start_positions.begin(), start_positions.end());
        std::sort(end_positions.begin(), end_positions.end());
        auto last = std::unique(start_positions.begin(), start_positions.end());
        start_positions.erase(last, start_positions.end());
        std::vector<int> unique_end_positions;
        for (const auto& start : start_positions) {
          auto pos = std::find(start_positions.begin(), start_positions.end(), start) - start_positions.begin();
          unique_end_positions.push_back(end_positions[pos]);
        }
        int edit_distance = cresult.editDistance;
        for (size_t i = 0; i < start_positions.size(); ++i) {
          UniqueAlignment ua = {edit_distance, start_positions[i], unique_end_positions[i]};
          uniqueAlignments.insert(ua);
        }
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1}; ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
#pragma omp critical
{
  std::string best_signature_part = query_name + ":" + std::to_string(best.edit_distance) + ":" + std::to_string(best.start_position) + ":" + std::to_string(best.end_position);
  temp_signature_parts.push_back(best_signature_part);
  if (!skip_second_best) {
    std::string second_best_signature_part = query_name + ":" + std::to_string(secondBest.edit_distance) + ":" + std::to_string(secondBest.start_position) + ":" + std::to_string(secondBest.end_position);
    temp_signature_parts.push_back(second_best_signature_part);
  }
}
        edlibFreeAlignResult(cresult);
      }
    }
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string &a, const std::string &b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      }
    );
    std::string final_signature = std::accumulate(temp_signature_parts.begin(), temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> std::string {
        return a + (a.length() > 0 ? "|" : "") + b;
      });
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

std::pair<std::vector<std::string>, std::vector<std::string>> parseAdapters(
    const Rcpp::CharacterVector& adapters) {
  std::vector<std::string> queries(adapters.size());
  std::vector<std::string> query_names(adapters.size());
  Rcpp::CharacterVector names = adapters.names();
  for (int i = 0; i < adapters.size(); ++i) {
    queries[i] = Rcpp::as<std::string>(adapters[i]);
    query_names[i] = Rcpp::as<std::string>(names[i]);
  }
  return std::make_pair(queries, query_names);
}

std::map<std::string, std::pair<double, double>> parseThresholds(
    const Rcpp::DataFrame& misalignment_threshold) {
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for (int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = std::make_pair((double) misal_threshold[i], (double) misal_sd[i]);
  }
  return null_dist_map;
}

void processChunk(
    const std::vector<SearchResult>& chunk,
    const std::vector<std::string>& queries,
    const std::vector<std::string>& query_names,
    const std::map<std::string, std::pair<double, double>>& null_dist_map,
    const Rcpp::DataFrame& read_layout,
    const Rcpp::DataFrame& misalignment_threshold,
    int nthreads,
    std::vector<std::string>& processed_signatures,
    std::vector<int>& original_indices) {
  
  std::vector<std::string> sequences;
  std::vector<std::string> ids;
  
  for (const auto& sr : chunk) {
    sequences.push_back(sr.line);
    ids.push_back(sr.read_id);
  }
  
  Rcpp::CharacterVector sigalign_results = sigalign(queries, query_names, sequences, ids, null_dist_map, nthreads);
  std::vector<std::string> sigalign_vector = Rcpp::as<std::vector<std::string>>(sigalign_results);
  
  std::vector<std::string> sigrun_results = sigrun_cpp(sigalign_vector, container, misalignment_threshold, nthreads, false);
  
  for (size_t i = 0; i < sigrun_results.size(); ++i) {
    processed_signatures.push_back(sigrun_results[i]);
    original_indices.push_back(i);
  }
}

void writeFasta(const std::string& fileName, const std::vector<std::string>& ids, const std::vector<std::string>& sequences) {
  std::ofstream outputFile(fileName);
  if (!outputFile.is_open()) {
    Rcpp::stop("Failed to open output file.");
  }
  for (size_t i = 0; i < ids.size(); ++i) {
    if(sequences[i].size() > 0){
      outputFile << ">" << ids[i] << "\n" << sequences[i] << "\n";
    } else {
      continue;
    }
    //Rcpp::Rcout << "Writing ID: " << ids[i] << " with sequence length: " << sequences[i].size() << std::endl;
  }
  outputFile.close();
}

void writeSigstrings(const std::string& sigstringsFilePath, const std::vector<std::string>& sigstrings) {
  std::ofstream sigstringsFile(sigstringsFilePath);
  if (!sigstringsFile.is_open()) {
    Rcpp::stop("Failed to open sigstrings file.");
  }
  for (const auto& sig : sigstrings) {
    sigstringsFile << sig << "\n";
  }
  sigstringsFile.close();
}

void sig_extraction_to_fasta(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold,
  const std::vector<std::string>& chunk_ids, const std::vector<std::string>& chunk_sequences,
  const Rcpp::CharacterVector& processed_sigstrings, const std::string& outputFilePath, bool verbose,
  const std::vector<int>& original_indices) {
  // Convert R DataFrame to C++ ReadLayout structure
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  // Process variable scanning and create position function map
  VarScan(readLayout, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(readLayout, verbose);
  // Convert processed_sigstrings to std::vector<std::string>
  std::vector<std::string> sigstrings = Rcpp::as<std::vector<std::string>>(processed_sigstrings);
  
  // Prepare the result vectors for writing to FASTA
  std::vector<std::string> result_ids;
  std::vector<std::string> result_sequences;
  
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < sigstrings.size(); ++i) {
    const std::string& signature = sigstrings[i];
    int original_index = original_indices[i];
    const std::string& sequence = chunk_sequences[original_index];
    const std::string& id = chunk_ids[original_index];
    
    // Parse the signature
    auto read_info_pos = signature.find_last_of('<');
    if (read_info_pos == std::string::npos) {
      Rcpp::Rcout << "Invalid signature format: " << signature << std::endl;
      continue;
    }
    
    std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2);
    std::stringstream lengthTypeStream(lengthTypeStr);
    std::string lengthStr, read_id, type;
    std::getline(lengthTypeStream, lengthStr, ':');
    std::getline(lengthTypeStream, read_id, ':');
    std::getline(lengthTypeStream, type);
    
    if (verbose) {
      Rcpp::Rcout << "Signature: " << signature << std::endl;
      Rcpp::Rcout << "Length: " << lengthStr << ", Read ID: " << read_id << ", Type: " << type << std::endl;
    }
    
    auto plus_pos = read_id.find('+');
    if (plus_pos != std::string::npos) {
      read_id = read_id.substr(0, plus_pos);
    }
    
    if (type == "undecided") {
      continue;
    }
    
    // Process the signature string
    std::stringstream ss(signature.substr(0, read_info_pos));
    std::string token;
    std::string new_id = read_id;
    std::string extracted_sequence;
    
    while (std::getline(ss, token, '|')) {
      std::string token_id;
      int editDistance, startPos, endPos;
      std::stringstream tokenStream(token);
      
      std::getline(tokenStream, token_id, ':');
      tokenStream >> editDistance;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> startPos;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> endPos;
      
      auto it = readLayout.get<id_tag>().find(token_id);
      if (it != readLayout.get<id_tag>().end() && it->type == "variable" && it->global_class != "read") {
        std::string extracted_seq;
        if (startPos > 0 && endPos >= startPos && (sequence.length() >= endPos)) {
          extracted_seq = sequence.substr(startPos - 1, endPos - startPos + 1); // Corrected substring extraction
          if (it->direction != "forward") {
            extracted_seq = revcomp_cpp(extracted_seq);
          }
        } else {
          continue;
        }
        std::string classId = token_id;
        if (classId.substr(0, 3) == "rc_") {
          classId = classId.substr(3);
        }
        new_id += "|" + classId + ":" + extracted_seq; // Append to new_id
      } else if (it != readLayout.get<id_tag>().end() && it->global_class == "read") {
        if (startPos > 0 && endPos >= startPos && (sequence.length() >= endPos)) {
          extracted_sequence = sequence.substr(startPos - 1, endPos - startPos + 1); // Corrected substring extraction
          if (it->direction != "forward") {
            extracted_sequence = revcomp_cpp(extracted_sequence);
          }
        }
      }
    }
#pragma omp critical
{
  //Rcpp::Rcout << "Appending ID: " << new_id << " with sequence length: " << extracted_sequence.size() << std::endl;
  result_ids.push_back(new_id);
  result_sequences.push_back(extracted_sequence);
}
  }
  // Write the results to a FASTA file
  writeFasta(outputFilePath, result_ids, result_sequences);
}

// [[Rcpp::export]]
void sigstream(const std::vector<std::string>& inputFilePaths,
  const std::string& sigstringsFilePath, Rcpp::CharacterVector adapters, 
  const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, 
  int chunkSize, int nthreads, const std::string& fastaOutputPath, 
  int max_sequences = -1, bool verbose = false) {
  
  auto parsed_adapters = parseAdapters(adapters);
  std::vector<std::string> queries = parsed_adapters.first;
  std::vector<std::string> query_names = parsed_adapters.second;
  std::map<std::string, std::pair<double, double>> null_dist_map = parseThresholds(misalignment_threshold);
  
  if (inputFilePaths.size() == 1) {
    std::string inputFilePath = inputFilePaths[0];
    
    gzFile fp = gzopen(inputFilePath.c_str(), "r");
    if (!fp) {
      Rcpp::stop("Failed to open input file.");
    }
    
    kseq_t* seq = kseq_init(fp);
    
    std::ofstream outputFile(sigstringsFilePath);
    if (!outputFile.is_open()) {
      Rcpp::stop("Failed to open output file.");
    }
    
    container = prep_read_layout_cpp(read_layout, misalignment_threshold, false);
    
    std::vector<std::string> chunk_ids;
    std::vector<std::string> chunk_sequences;
    std::vector<std::string> processed_signatures;
    std::vector<int> original_indices;
    
    bool is_fastq = false;
    int read_count = 0;
    while (true) {
      std::vector<SearchResult> chunk = readFastqChunk(seq, fp, chunkSize, is_fastq, read_count, max_sequences);
      if (chunk.empty()) break;
      processChunk(chunk, queries, query_names, null_dist_map, 
        read_layout, misalignment_threshold, nthreads, processed_signatures, original_indices);
      for (const auto& sr : chunk) {
        chunk_ids.push_back(sr.read_id);
        chunk_sequences.push_back(sr.line);
      }
      if (max_sequences != -1 && read_count >= max_sequences) break;
    }
    
    kseq_destroy(seq);
    gzclose(fp);
    outputFile.close();
    
    // Write all signatures to a file
    writeSigstrings(sigstringsFilePath, processed_signatures);
    
    // Call the sig_extraction_to_fasta function
    sig_extraction_to_fasta(read_layout, misalignment_threshold, chunk_ids, chunk_sequences, 
      Rcpp::wrap(processed_signatures), fastaOutputPath, verbose, original_indices);
    
  } else {
    for (const auto& inputFilePath : inputFilePaths) {
      gzFile fp = gzopen(inputFilePath.c_str(), "r");
      if (!fp) {
        Rcpp::stop("Failed to open input file: " + inputFilePath);
      }
      
      kseq_t* seq = kseq_init(fp);
      
      std::ofstream outputFile(sigstringsFilePath, std::ios::app);
      if (!outputFile.is_open()) {
        Rcpp::stop("Failed to open output file.");
      }
      
      container = prep_read_layout_cpp(read_layout, misalignment_threshold, false);
      
      std::vector<std::string> chunk_ids;
      std::vector<std::string> chunk_sequences;
      std::vector<std::string> processed_signatures;
      std::vector<int> original_indices;
      
      bool is_fastq = false;
      int read_count = 0;
      while (true) {
        std::vector<SearchResult> chunk = readFastqChunk(seq, fp, chunkSize, is_fastq, read_count, max_sequences);
        if (chunk.empty()) break;
        processChunk(chunk, queries, query_names, null_dist_map, 
          read_layout, misalignment_threshold, nthreads, processed_signatures, original_indices);
        for (const auto& sr : chunk) {
          chunk_ids.push_back(sr.read_id);
          chunk_sequences.push_back(sr.line);
        }
        if (max_sequences != -1 && read_count >= max_sequences) break;
      }
      
      kseq_destroy(seq);
      gzclose(fp);
      outputFile.close();
      
      // Write all signatures to a file
      writeSigstrings(sigstringsFilePath, processed_signatures);
      
      // Call the sig_extraction_to_fasta function
      sig_extraction_to_fasta(read_layout, misalignment_threshold, chunk_ids, chunk_sequences, 
        Rcpp::wrap(processed_signatures), fastaOutputPath, verbose, original_indices);
    }
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame misalignment_stream(const std::string& fastq_path, Rcpp::CharacterVector adapters, 
  int nthreads = 1, int max_sequences = -1) {
  
  gzFile fastq_file = gzopen(fastq_path.c_str(), "rb");
  if (!fastq_file) {
    Rcpp::stop("Failed to open the FASTQ file.");
  }
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  // Filter out queries with "poly" in their names
  std::vector<std::string> filtered_queries;
  Rcpp::CharacterVector filtered_query_names;
  for (int i = 0; i < query_names.size(); ++i) {
    std::string name = Rcpp::as<std::string>(query_names[i]);
    if (name.find("poly") == std::string::npos) {
      filtered_queries.push_back(queries[i]);
      filtered_query_names.push_back(query_names[i]);
    }
  }
  
  std::map<int, std::vector<AlignmentInfo>> all_alignments;
  size_t initial_buffer_size = 98304;
  std::vector<char> buffer(initial_buffer_size);
  
  std::string line, id, sequence, plus, quality;
  int sequence_counter = 1;
  int processed_sequences = 0;
  auto read_line = [&](gzFile file, std::vector<char>& buffer) -> std::string {
    std::string result;
    while (gzgets(file, buffer.data(), buffer.size()) != Z_NULL) {
      result += std::string(buffer.data());
      if (result.back() == '\n') break;
      buffer.resize(buffer.size() * 2); // Double the buffer size if not complete
    }
    result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
    return result;
  };
  
  omp_set_num_threads(nthreads);
  
  while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL &&
    (max_sequences == -1 || processed_sequences < max_sequences)) {
    line = std::string(buffer.data());
    if (line[0] != '@') {
      Rcpp::stop("Invalid FASTQ format: Missing '@' at the beginning of the ID line.");
    }
    
    id = line.substr(1, line.find(' ') - 1); // Extract ID between '@' and first space
    sequence = read_line(fastq_file, buffer);
    if (sequence.empty()) break;
    
    plus = read_line(fastq_file, buffer);
    if (plus.empty() || plus[0] != '+') {
      Rcpp::Rcout << "Read too big or FASTQ formatting error. Skipping read.\n";
      // Skip the rest of the current sequence
      while (gzgets(fastq_file, buffer.data(), buffer.size()) != Z_NULL && buffer[0] != '@');
      continue;
    }
    
    quality = read_line(fastq_file, buffer);
    if (quality.empty()) break;
    
    // Process the accumulated sequence
    std::vector<AlignmentInfo> alignments(filtered_queries.size());
    
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < filtered_queries.size(); ++j) {
      const auto& query = filtered_queries[j];
      
      EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
      std::set<UniqueAlignment> uniqueAlignments;
      char* cquery = const_cast<char*>(query.c_str());
      char* csequence = const_cast<char*>(sequence.c_str());
      EdlibAlignResult cresult = edlibAlign(cquery, query.size(), csequence, sequence.size(), config);
      std::vector<int> start_positions(cresult.startLocations, cresult.startLocations + cresult.numLocations);
      std::vector<int> end_positions(cresult.endLocations, cresult.endLocations + cresult.numLocations);
      
      for (int& pos : start_positions) {
        ++pos;
      }
      for (int& pos : end_positions) {
        ++pos;
      }
      
      std::sort(start_positions.begin(), start_positions.end());
      std::sort(end_positions.begin(), end_positions.end());
      
      auto last = std::unique(start_positions.begin(), start_positions.end());
      start_positions.erase(last, start_positions.end());
      
      std::vector<int> unique_end_positions;
      for (const auto& start : start_positions) {
        auto pos = std::find(start_positions.begin(), start_positions.end(), start) - start_positions.begin();
        unique_end_positions.push_back(end_positions[pos]);
      }
      
      int edit_distance = cresult.editDistance;
      
      for (size_t k = 0; k < start_positions.size(); ++k) {
        UniqueAlignment ua = {edit_distance, start_positions[k], unique_end_positions[k]};
        uniqueAlignments.insert(ua);
      }
      auto it = uniqueAlignments.begin();
      UniqueAlignment best = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
      ++it;
      UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? *it : UniqueAlignment{-1, -1, -1};
      
      alignments[j] = {
        best.edit_distance,
        best.start_position,
        best.end_position,
        secondBest.edit_distance
      };
      edlibFreeAlignResult(cresult);
    }
    
    all_alignments[sequence_counter] = alignments;
    ++sequence_counter;
    ++processed_sequences;
  }
  gzclose(fastq_file);
  std::vector<int> ids;
  std::vector<std::string> query_ids;  // to store the query ID for each alignment
  std::vector<int> best_edit_distances;
  std::vector<int> best_start_positions;
  std::vector<int> best_stop_positions;
  std::vector<int> second_edit_distances;
  
  for (const auto& pair : all_alignments) {
    int sequence_id = pair.first;
    const auto& alignments = pair.second;
    for (size_t i = 0; i < alignments.size(); ++i) {
      ids.push_back(sequence_id);
      query_ids.push_back(Rcpp::as<std::string>(filtered_query_names[i]));  // get the filtered query names correctly
      best_edit_distances.push_back(alignments[i].best_edit_distance);
      best_start_positions.push_back(alignments[i].best_start_pos);
      best_stop_positions.push_back(alignments[i].best_stop_pos);
      second_edit_distances.push_back(alignments[i].second_edit_distance);
    }
  }
  Rcpp::DataFrame result = Rcpp::DataFrame::create(
    Rcpp::Named("id") = ids,
    Rcpp::Named("query_id") = query_ids,
    Rcpp::Named("best_edit_distance") = best_edit_distances,
    Rcpp::Named("best_start_pos") = best_start_positions,
    Rcpp::Named("best_stop_pos") = best_stop_positions,
    Rcpp::Named("second_edit_distance") = second_edit_distances
  );
  return result;
}
