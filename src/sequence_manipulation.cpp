#include "rad.h"
using namespace std;
using namespace Rcpp;

std::string revcomp_cpp(const std::string& sequence) {
  std::unordered_map<char, char> complement {
    {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
    {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'},
    {'N', 'N'}, {'n', 'n'}
  };
  std::string reversed(sequence.rbegin(), sequence.rend());
  for (char& nucleotide : reversed) {
    nucleotide = complement[nucleotide];
  }
  return reversed;
}

// [[Rcpp::export]]
Rcpp::StringVector revcomp(Rcpp::StringVector sequences) {
  int n = sequences.size();
  Rcpp::StringVector results(n);
  int j = 0; // Counter for results
  for (int i = 0; i < n; ++i) {
    if (sequences[i] == NA_STRING) {
      results[j] = NA_STRING;
      continue;
    }
    results[j] = revcomp_cpp(Rcpp::as<std::string>(sequences[i]));
    ++j;
  }
  return results;
}

std::vector<int64_t> sequence_to_bits_cpp(const std::string& sequence) {
  const int chunk_size = 32; // Number of bases per chunk
  std::vector<int64_t> results;
  for (size_t i = 0; i < sequence.size(); i += chunk_size) {
    int64_t result = 0;
    size_t end = std::min(i + chunk_size, sequence.size());
    for (size_t j = i; j < end; ++j) {
      char c = sequence[j];
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break;
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      }
    }
    results.push_back(result);
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector sequence_to_bits(Rcpp::StringVector sequences) {
  int n = sequences.size();
  std::vector<int64_t> results;
  for(int i = 0; i < n; ++i) {
    std::string sequence = Rcpp::as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = sequence_to_bits_cpp(sequence);
    results.insert(results.end(), chunk_results.begin(), chunk_results.end());
  }
  return Rcpp::toInteger64(results);
}

// [[Rcpp::export]]
List sequence_to_bitlist(Rcpp::StringVector sequences) {
  int n = sequences.size();
  List results(n);
  for(int i = 0; i < n; ++i) {
    std::string sequence = as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = sequence_to_bits_cpp(sequence);
    results[i] = wrap(chunk_results);
  }
  return results;
}

std::vector<int64_t> qual_to_bits_cpp(const std::string& qual) {
  std::vector<int64_t> result;
  const size_t charsPerInt64 = 8; // Number of characters to pack into an int64_t
  size_t totalInt64s = (qual.size() + charsPerInt64 - 1) / charsPerInt64;
  for (size_t i = 0; i < totalInt64s; ++i) {
    int64_t packedValue = 0;
    for (size_t j = 0; j < charsPerInt64; ++j) {
      size_t charIndex = i * charsPerInt64 + j;
      if (charIndex < qual.size()) {
        // Shift the packedValue left by 8 bits and add the new character
        packedValue = (packedValue << 8) | static_cast<unsigned char>(qual[charIndex]);
      }
    }
    result.push_back(packedValue);
  }
  
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector qual_to_bits(Rcpp::StringVector sequences) {
  int n = sequences.size();
  std::vector<int64_t> results;
  for(int i = 0; i < n; ++i) {
    std::string sequence = Rcpp::as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = qual_to_bits_cpp(sequence);
    results.insert(results.end(), chunk_results.begin(), chunk_results.end());
  }
  return Rcpp::toInteger64(results);
}

// [[Rcpp::export]]
List qual_to_bitlist(StringVector sequences) {
  int n = sequences.size();
  List results(n);
  for(int i = 0; i < n; ++i) {
    std::string sequence = as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = qual_to_bits_cpp(sequence);
    results[i] = wrap(chunk_results);
  }
  return results;
}

std::string bits_to_sequence_cpp(const std::vector<int64_t>& input, int sequence_length) {
  int n = input.size();
  std::string final_sequence;
  int chunk_size = 32;
  for (int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    std::string result;
    int bases_to_decode = (i < n - 1) ? chunk_size : sequence_length % chunk_size;
    bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
    for (int j = 0; j < bases_to_decode; ++j) {
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    final_sequence += result;
  }
  return final_sequence;
}

std::vector<unsigned int> bits_to_uint_cpp(const std::vector<int64_t>& input, int sequence_length) {
  int n = input.size();
  std::vector<unsigned int> final_sequence;
  int chunk_size = 32; // Assuming 32 bases per int64_t value
  for (int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    int bases_to_decode = (i < n - 1) ? chunk_size : sequence_length % chunk_size;
    bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
    for (int j = 0; j < bases_to_decode; ++j) {
      unsigned int base = int64_code & 3; // Extract the last two bits
      final_sequence.push_back(base); // Add the base as an unsigned int to the sequence
      int64_code >>= 2; // Move to the next base
    }
  }
  return final_sequence;
}

// [[Rcpp::export]]
Rcpp::StringVector bits_to_sequence(Rcpp::NumericVector input, Rcpp::IntegerVector sequence_length = 16) {
  // Convert the entire NumericVector to std::vector<int64_t>
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(input);
  int n = cpp_input_list.size(); // Number of sequences to process
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  Rcpp::StringVector final_sequences(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int64_t> bits_input = {cpp_input_list[i]}; // Single element vector
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    // Process the single sequence
    std::string sequence = bits_to_sequence_cpp(bits_input, seqlength);
    final_sequences[i] = sequence;
  }
  return final_sequences;
}

std::string bits_to_hex_cpp(const std::vector<int64_t>& input, int sequence_length) {
  std::string final_sequence;
  const int chunk_size = 32;
  for (size_t i = 0; i < input.size(); ++i) {
    int64_t int64_code = input[i];
    std::string result;
    
    int k_mers_to_decode = (i < input.size() - 1) ? chunk_size : sequence_length % chunk_size;
    k_mers_to_decode = (k_mers_to_decode == 0) ? chunk_size : k_mers_to_decode;
    
    for (int j = 0; j < k_mers_to_decode; ++j) {
      int k_mer = int64_code & 0xF;
      if (k_mer < 10) {
        result.insert(result.begin(), '0' + k_mer);
      } else {
        result.insert(result.begin(), 'A' + (k_mer - 10));
      }
      int64_code >>= 4;
    }
    final_sequence += result;
  }
  return final_sequence;
}

// [[Rcpp::export]]
StringVector bits_to_hex(Rcpp::NumericVector input, int sequence_length = 16) {
  Rcpp::StringVector hex_sequences(input.size());
  std::vector<int64_t> cpp_input = Rcpp::as<std::vector<int64_t>>(input);
  int slen = sequence_length/2; //for consistency and compression's sake, we convert from a rep of 2 bits->one char to 4 bits->one char
  for(size_t i = 0; i < cpp_input.size(); ++i) {
    std::vector<int64_t> single_input = {cpp_input[i]};
    std::string hex_sequence = bits_to_hex_cpp(single_input, slen);
    hex_sequences[i] = hex_sequence;
  }
  return hex_sequences;
}

std::vector<int64_t> revcomp_bits_cpp(const std::vector<int64_t>& input, int sequence_length = 16) {
  std::vector<int64_t> reversed(input.rbegin(), input.rend()); // Reverse the vector
  for (int64_t& chunk : reversed) {
    int64_t result = 0;
    for (int i = 0; i < sequence_length; i++) {
      result <<= 2; // Shift result to make room for the new bits
      int64_t base = chunk & 3; // Extract the two rightmost bits
      switch (base) {
      case 0: result |= 2; break; // A (00) becomes T (10)
      case 1: result |= 3; break; // C (01) becomes G (11)
      case 2: result |= 0; break; // T (10) becomes A (00)
      case 3: result |= 1; break; // G (11) becomes C (01)
      }
      chunk >>= 2; // Move to the next base in the chunk
    }
    // Assign the reversed and complemented chunk back
    chunk = result;
  }
  return reversed;
}

// [[Rcpp::export]]
Rcpp::NumericVector revcomp_bits(Rcpp::NumericVector input, Rcpp::IntegerVector sequence_length = 16) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(input);
  int n = cpp_input_list.size(); // Number of sequences to process
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  std::vector<int64_t> results;
  for (int i = 0; i < n; ++i) {
    std::vector<int64_t> bits_input = {cpp_input_list[i]}; // Single element vector
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    std::vector<int64_t> subresult = revcomp_bits_cpp(bits_input, seqlength);
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::toInteger64(results);
}

int hamming_bits_cpp(int64_t a, int64_t b, int sequence_length) {
  int64_t xor_result = a ^ b; // XOR to find differing bits
  int count = 0;
  for (int i = 0; i < sequence_length; i++) {
    if ((xor_result & 3) != 0) {
      count++; 
    }
    xor_result >>= 2;
  }
  return count;
}

// [[Rcpp::export]]
IntegerVector hamming_bits(NumericVector int64_a, NumericVector int64_b, IntegerVector sequence_length = 16) {
  int n_a = int64_a.size();
  int n_b = int64_b.size();
  int n_max = std::max(n_a, n_b); // Determine the maximum length to set as the loop limit
  IntegerVector distances(n_max);
  bool use_single_length = sequence_length.size() == 1;
  int seq_length = use_single_length ? sequence_length[0] : 0; // Use this if only one length is provided
  std::vector<int64_t> int64_avec = Rcpp::as<std::vector<int64_t>>(int64_a);
  std::vector<int64_t> int64_bvec = Rcpp::as<std::vector<int64_t>>(int64_b);
  for (int i = 0; i < n_max; ++i) {
    // Handle the scenario where one vector has a single value and the other has multiple values
    int64_t a = (n_a == 1) ? int64_avec[0] : int64_avec[i]; // Use the single value or the current value
    int64_t b = (n_b == 1) ? int64_bvec[0] : int64_bvec[i]; // Use the single value or the current value
    int cur_seq = use_single_length ? seq_length : sequence_length[std::min(i, static_cast<int>(sequence_length.size()) - 1)];
    // Use the single length or individual lengths, ensuring we don't go out of bounds
    distances[i] = hamming_bits_cpp(a, b, cur_seq);
  } 
  return distances;
}

std::vector<int64_t> generate_bit_mutations_cpp(int64_t sequence, int sequence_length) {
  std::vector<int64_t> mutations;
  for (int position = 0; position < sequence_length; ++position) { // 16 nucleotides in the sequence
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11; // Extract the original nucleotide
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) { // Ensure we're generating a mutation
        int64_t mutation_mask = new_nucleotide << (2 * position); // Set the new nucleotide
        int64_t clear_mask = ~(0b11LL << (2 * position)); // Clear the original nucleotide position
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask; // Apply mutation
        mutations.push_back(mutated_sequence);
      }
    }
  }
  return mutations;
}

// [[Rcpp::export]]
Rcpp::NumericVector generate_bit_mutations(Rcpp::NumericVector sequences, IntegerVector sequence_length = 16) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  int n = cpp_input_list.size(); // Number of sequences to process
  std::vector<int64_t> results;
  for (int i = 0; i < n; ++i) {
    int64_t bits_input = cpp_input_list[i]; // Directly use the int64_t value
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    std::vector<int64_t> subresult = generate_bit_mutations_cpp(bits_input, seqlength); // Corrected call
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::toInteger64(results);
}

std::vector<int64_t> generate_recursive_mutations_cpp(int64_t sequence, int mutation_rounds, int sequence_length = 16) {
  std::unordered_set<int64_t> mutations{sequence};
  std::unordered_set<int64_t> new_mutations;
  for (int round = 0; round < mutation_rounds; ++round) {
    for (const auto& seq : mutations) {
      std::vector<int64_t> current_mutations = generate_bit_mutations_cpp(seq, sequence_length);
      new_mutations.insert(current_mutations.begin(), current_mutations.end());
    }
    mutations.insert(new_mutations.begin(), new_mutations.end());
  }
  mutations.erase(sequence);
  return std::vector<int64_t>(mutations.begin(), mutations.end());
}
// [[Rcpp::export]]
Rcpp::NumericVector generate_recursive_mutations(Rcpp::NumericVector sequences, int mutation_rounds = 2, int sequence_length = 16) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  std::vector<int64_t> results;
  for (int64_t bits_input : cpp_input_list) {
    std::vector<int64_t> subresult = generate_recursive_mutations_cpp(bits_input, mutation_rounds, sequence_length); 
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}

std::vector<int64_t> extract_subsequence_cpp(const std::vector<int64_t>& input, int start_pos, int stop_pos) {
  std::vector<int64_t> result;
  // Calculate start and end elements in the vector
  int start_element = start_pos / 32;
  int end_element = stop_pos / 32;
  int start_bit = (start_pos % 32) * 2;
  int end_bit = ((stop_pos % 32) + 1) * 2;
  if (start_element == end_element) {
    int64_t mask = ((1LL << (end_bit - start_bit)) - 1) << (64 - end_bit);
    result.push_back((input[start_element] & mask) >> (64 - end_bit));
    return result;
  }
  for (int i = start_element; i <= end_element; ++i) {
    int64_t mask;
    if (i == start_element) {
      mask = (1LL << (64 - start_bit)) - 1;
      result.push_back((input[i] & mask) >> start_bit);
    } else if (i == end_element) {
      mask = ((1LL << end_bit) - 1) << (64 - end_bit);
      result.push_back((input[i] & mask) >> (64 - end_bit));
    } else {
      result.push_back(input[i]);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector bit_substring(Rcpp::NumericVector bit_sequences, int start_pos, int stop_pos) {
  start_pos = start_pos - 1;
  stop_pos = stop_pos -1;
  //fixing this index for r to c++
  std::vector<int64_t> sequence = Rcpp::as<std::vector<int64_t>>(bit_sequences);
  std::vector<int64_t> extracted = extract_subsequence_cpp(sequence, start_pos, stop_pos);
  return Rcpp::toInteger64(extracted);
}

std::string to_binary_string(int64_t value, int total_length) {
  std::string binary = "";
  for(int i = total_length - 1; i >= 0; --i) {
    binary += ((value >> i) & 1) ? "1" : "0";
  }
  return binary;
}

std::vector<int64_t> kmer_circ_cpp(int64_t sequence, int base_length, int kmerLength, bool verbose = false) {
  if(base_length > 32) {
    throw std::runtime_error("Base length exceeds 32 base pairs limit for int64_t representation.");
  }
  int64_t circularizedSequence = sequence | (sequence << (base_length * 2));
  if(verbose) {
    Rcpp::Rcout << "Original sequence in binary: " << to_binary_string(sequence, base_length * 2) << "\n";
    Rcpp::Rcout << "Circularized sequence in binary: " << to_binary_string(circularizedSequence, base_length * 4) << "\n";
  }
  std::vector<int64_t> kmers;
  int totalBases = base_length * 2;
  for (int start_pos = 1; start_pos <= totalBases - kmerLength - 1; ++start_pos) {
    int stop_pos = start_pos + kmerLength - 1;
    std::vector<int64_t> subsequence = extract_subsequence_cpp({circularizedSequence}, start_pos, stop_pos);
    if (!subsequence.empty()) {
      kmers.push_back(subsequence[0]);
      if(verbose) {
        Rcpp::Rcout << "Kmer " << start_pos + 1 << " in binary: " << to_binary_string(subsequence[0], kmerLength * 2) << "\n";
      }
    }
  }
  return kmers;
}

std::vector<int64_t> kmer_circ_cpp_v2(int64_t sequence, int base_length, int kmerLength, bool verbose = false) {
  if (base_length > 32) {
    throw std::runtime_error("Base length exceeds 32 base pairs limit for int64_t representation.");
  }
  int64_t circularizedSequence = sequence | (sequence << (base_length * 2));
  if (verbose) {
    Rcpp::Rcout << "Original sequence in binary: " << to_binary_string(sequence, base_length * 2) << "\n";
    Rcpp::Rcout << "Circularized sequence in binary: " << to_binary_string(circularizedSequence, base_length * 4) << "\n";
  }
  std::vector<int64_t> kmers;
  int totalBases = base_length * 2;
  int endRange = totalBases - kmerLength - 1;
  // Extract the first two kmers
  for (int start_pos = 1; start_pos <= 2; ++start_pos) {
    int stop_pos = start_pos + kmerLength - 1;
    std::vector<int64_t> subsequence = extract_subsequence_cpp({circularizedSequence}, start_pos, stop_pos);
    if (!subsequence.empty()) {
      kmers.push_back(subsequence[0]);
      if (verbose) {
        Rcpp::Rcout << "Kmer " << start_pos << " in binary: " << to_binary_string(subsequence[0], kmerLength * 2) << "\n";
      }
    }
  }
  // Extract the last two kmers
  for (int start_pos = endRange - 1; start_pos <= endRange; ++start_pos) {
    int stop_pos = start_pos + kmerLength - 1;
    std::vector<int64_t> subsequence = extract_subsequence_cpp({circularizedSequence}, start_pos, stop_pos);
    if (!subsequence.empty()) {
      kmers.push_back(subsequence[0]);
      if (verbose) {
        Rcpp::Rcout << "Kmer " << start_pos << " in binary: " << to_binary_string(subsequence[0], kmerLength * 2) << "\n";
      }
    }
  }
  return kmers;
}

// [[Rcpp::export]]
Rcpp::NumericVector kmer_circ(SEXP sequence, int base_length = 16, int kmerLength = 16, bool verbose = false) {
  std::vector<int64_t> sequenceVec = Rcpp::as<std::vector<int64_t>>(Rcpp::NumericVector(sequence));
  if (sequenceVec.size() != 1) {
    stop("Expected a single integer64 value.");
  }
  std::vector<int64_t> kmers = kmer_circ_cpp_v2(sequenceVec[0], base_length, kmerLength, verbose);
  return Rcpp::wrap(kmers);
}

int hamming_distance(int64_t a, int64_t b) {
  return std::bitset<64>(a ^ b).count();
}
// Helper function to convert a DNA sequence to int64
int64_t seq_to_int64(const std::string& seq) {
  int64_t result = 0;
  for (char nucleotide : seq) {
    result <<= 2;
    switch (nucleotide) {
    case 'A':
      result |= 0;
      break;
    case 'C':
      result |= 1;
      break;
    case 'G':
      result |= 2;
      break;
    case 'T':
      result |= 3;
      break;
    default:
      Rcpp::stop("Invalid nucleotide in sequence");
    }
  }
  return result;
}

// [[Rcpp::export]]
void bc_correct_module(std::string input_file,
                       std::string output_file,
                       bool verbose = false,
                       int nthreads = 1,
                       int maxDistance = 2,
                       int sequence_length = 16,
                       int depth = 2,
                       int breadth = 1,
                       bool high_speed = true) {
  // Read input CSV file
  std::ifstream infile(input_file);
  if (!infile.is_open()) {
    Rcpp::stop("Could not open input file.");
  }
  
  std::ofstream outfile(output_file);
  if (!outfile.is_open()) {
    Rcpp::stop("Could not open output file.");
  }
  
  // Write the header to the output file
  outfile << "seq,count,filtered,corrected_seq\n";
  
  std::vector<std::string> seq;
  std::vector<int> count;
  std::vector<std::string> filtered;
  std::vector<int64_t> int64_seq;
  std::vector<int64_t> rc_int64seq;
  
  std::string line;
  
  // Skip the header row
  if (!std::getline(infile, line)) {
    Rcpp::stop("Input file is empty.");
  }
  
  while (std::getline(infile, line)) {
    std::istringstream ss(line);
    std::string tmp;
    
    std::getline(ss, tmp, ',');
    seq.push_back(tmp.empty() ? "N" : tmp);  // Fill empty entries with "N"
    
    std::getline(ss, tmp, ',');
    try {
      count.push_back(tmp.empty() ? 0 : std::stoi(tmp));  // Fill empty entries with 0
    } catch (const std::invalid_argument& e) {
      Rcpp::stop("Error parsing count: " + tmp);
    }
    
    std::getline(ss, tmp, ',');
    filtered.push_back(tmp.empty() ? "unknown" : tmp);  // Fill empty entries with "unknown"
    
    // Skip pois_dist column
    std::getline(ss, tmp, ',');
    
    std::getline(ss, tmp, ',');
    try {
      int64_seq.push_back(tmp.empty() ? 0 : std::stoll(tmp));  // Fill empty entries with 0
    } catch (const std::invalid_argument& e) {
      Rcpp::stop("Error parsing int64_seq: " + tmp);
    }
    
    std::getline(ss, tmp, ',');
    try {
      rc_int64seq.push_back(tmp.empty() ? 0 : std::stoll(tmp));  // Fill empty entries with 0
    } catch (const std::invalid_argument& e) {
      Rcpp::stop("Error parsing rc_int64seq: " + tmp);
    }
  }
  infile.close();
  
  // Debugging: Print the parsed input
  if (verbose) {
    Rcpp::Rcout << "Parsed " << seq.size() << " entries from the input file.\n";
  }
  
  // Generate the whitelist from the pois_validated_barcode
  std::unordered_map<int64_t, std::pair<std::string, double>> wl;
  
  for (size_t i = 0; i < seq.size(); ++i) {
    if (filtered[i] == "pois_validated_barcode") {
      wl[int64_seq[i]] = std::make_pair(seq[i], double(count[i]));
    }
  }
  
  // Collapse false positives in the whitelist
  for (size_t i = 0; i < seq.size(); ++i) {
    if (filtered[i] == "whitelist_barcode") {
      int64_t wl_barcode = int64_seq[i];
      for (const auto& validated : wl) {
        if (hamming_distance(wl_barcode, validated.first) <= maxDistance) {
          wl[validated.first].second = std::max(wl[validated.first].second, double(count[i]));
          break;
        }
      }
    }
  }
  
  int corrected_count = 0; // Variable to count the number of corrections
  int scanned_count = 0;   // Variable to count the number of scanned barcodes
  
  omp_set_num_threads(nthreads);
  
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < seq.size(); ++i) {
    if (filtered[i] == "barcode_to_correct") {
      int64_t barcode = int64_seq[i];
      double highest_poisson_score = -1.0;
      int best_distance = std::numeric_limits<int>::max();
      int64_t best_match = 0;
      
      // First check against pois_validated_barcodes
      for (const auto& validated : wl) {
        int distance = hamming_distance(barcode, validated.first);
        if (distance <= maxDistance) {
          double poisson_score = validated.second.second;
          if (poisson_score > highest_poisson_score ||
              (poisson_score == highest_poisson_score && distance < best_distance)) {
            best_match = validated.first;
            highest_poisson_score = poisson_score;
            best_distance = distance;
          }
        }
      }
      
      // If no match found, check against whitelist_barcodes
      if (best_match == 0) {
        for (size_t j = 0; j < seq.size(); ++j) {
          if (filtered[j] == "whitelist_barcode") {
            int64_t wl_barcode = int64_seq[j];
            int distance = hamming_distance(barcode, wl_barcode);
            if (distance <= maxDistance) {
              double poisson_score = double(count[j]);
              if (poisson_score > highest_poisson_score ||
                  (poisson_score == highest_poisson_score && distance < best_distance)) {
                best_match = wl_barcode;
                highest_poisson_score = poisson_score;
                best_distance = distance;
              }
            }
          }
        }
      }
      
      // If no match found, generate mutations and circular sequences
      if (best_match == 0) {
        auto mutations = generate_recursive_mutations_cpp(barcode, depth, sequence_length);
        for (auto& mutation : mutations) {
          if (wl.find(mutation) != wl.end()) {
            int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
            if (dl_distance <= maxDistance) {
              double poisson_score = wl[mutation].second;
              if (poisson_score > highest_poisson_score ||
                  (poisson_score == highest_poisson_score && dl_distance < best_distance)) {
                best_match = mutation;
                highest_poisson_score = poisson_score;
                best_distance = dl_distance;
              }
            }
          }
        }
        if (best_match == 0) {
          auto circularizedSequences = kmer_circ_cpp(barcode, sequence_length, sequence_length, false);
          for (auto& circSeq : circularizedSequences) {
            auto circMutations = generate_recursive_mutations_cpp(circSeq, breadth, sequence_length);
            for (auto& mutation : circMutations) {
              if (wl.find(mutation) != wl.end()) {
                int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
                if (dl_distance <= maxDistance) {
                  double poisson_score = wl[mutation].second;
                  if (poisson_score > highest_poisson_score ||
                      (poisson_score == highest_poisson_score && dl_distance < best_distance)) {
                    best_match = mutation;
                    highest_poisson_score = poisson_score;
                    best_distance = dl_distance;
                  }
                }
              }
            }
          }
        }
      }
      
      // Update the output file within the critical section
#pragma omp critical
{
  if (best_match != 0) {
    wl[best_match].second++;
    filtered[i] = "corrected_barcode";
    outfile << seq[i] << ',' << count[i] << ",corrected_barcode," << wl[best_match].first << '\n';
    corrected_count++;
  } else {
    outfile << seq[i] << ',' << count[i] << ',' << filtered[i] << ',' << "" << '\n';
  }
  
  scanned_count++;
  if (scanned_count % 10000 == 0) {
    Rcpp::Rcout << scanned_count << " barcodes scanned, " << corrected_count << " barcodes corrected.\n";
  }
}
    }
  }
  
  Rcpp::Rcout << "Total barcodes corrected: " << corrected_count << "\n";
  outfile.close();
}