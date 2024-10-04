#include "rad.h"
#include "stringdist_api.h"

using namespace std;
using namespace Rcpp;

KSEQ_INIT(gzFile, gzread)

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
  int slen = sequence_length/2; 
  //for consistency and compression's sake, we convert from a rep of 2 bits->one char to 4 bits->one char
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

// Function to calculate the number of differing bits for int64_t sequences
int hamming_distance(int64_t a, int64_t b) {
  int bit_diff_count = std::bitset<64>(a ^ b).count();
  return std::ceil(bit_diff_count / 2.0);  // Base Hamming distance
}

// [[Rcpp::export]]
int hamming_distance_rcpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  // Convert Rcpp::NumericVector (integer64) to std::vector<int64_t>
  std::vector<int64_t> a_int = Rcpp::fromInteger64(a);
  std::vector<int64_t> b_int = Rcpp::fromInteger64(b);
  
  if (a_int.size() != 1 || b_int.size() != 1) {
    Rcpp::stop("Expecting single integer64 values.");
  }
  
  // Return the Hamming distance between the first (and only) elements
  return hamming_distance(a_int[0], b_int[0]);
}

// Function to generate mutations for a single sequence
std::vector<int64_t> generate_quaternary_mutations_cpp_heap(int64_t sequence, int sequence_length) {
  // Allocate mutations vector on the heap
  auto mutations = std::make_unique<std::vector<int64_t>>();
  for (int position = 0; position < sequence_length; ++position) {
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11;
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) {
        int64_t mutation_mask = new_nucleotide << (2 * position);
        int64_t clear_mask = ~(0b11LL << (2 * position));
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask;
        mutations->push_back(mutated_sequence);
      }
    }
  }
  // Return the vector stored on the heap
  return *mutations;
}

// Function to generate mutations for a single sequence
std::vector<int64_t> generate_quaternary_mutations_cpp(int64_t sequence, int sequence_length) {
  std::vector<int64_t> mutations;
  for (int position = 0; position < sequence_length; ++position) {
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11;
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) {
        int64_t mutation_mask = new_nucleotide << (2 * position);
        int64_t clear_mask = ~(0b11LL << (2 * position));
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask;
        mutations.push_back(mutated_sequence);
      }
    }
  }
  return mutations;
}

// [[Rcpp::export]]
Rcpp::NumericVector generate_quaternary_mutations(Rcpp::NumericVector sequences, Rcpp::IntegerVector sequence_length = 32) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0;
  std::vector<int64_t> results;
  for (size_t i = 0; i < cpp_input_list.size(); ++i) {
    int64_t quaternary_input = cpp_input_list[i];
    int seqlength = use_sl ? seq_length : sequence_length[i];
    std::vector<int64_t> subresult = generate_quaternary_mutations_cpp(quaternary_input, seqlength);
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}

// Function to generate recursive mutations
std::vector<int64_t> generate_recursive_quaternary_mutations_cpp(
    int64_t sequence, 
  int mutation_rounds, 
  int sequence_length = 32) {
  std::unordered_set<int64_t> all_mutations{sequence};
  std::vector<int64_t> current_round{sequence};
  for (int round = 0; round < mutation_rounds; ++round) {
    std::vector<int64_t> next_round;
    for (const auto& seq : current_round) {
      auto mutations = generate_quaternary_mutations_cpp(seq, sequence_length);
      for (const auto& mutation : mutations) {
        if (all_mutations.insert(mutation).second) {
          next_round.push_back(mutation);
        }
      }
    }
    if (next_round.empty()) break;
    current_round = std::move(next_round);
  }
  all_mutations.erase(sequence);
  return std::vector<int64_t>(all_mutations.begin(), all_mutations.end());
}

std::vector<int64_t> generate_recursive_quaternary_mutations_cpp_v2(
    int64_t sequence,
    int mutation_rounds,
    int sequence_length = 32) {
  std::unordered_set<int64_t> all_mutations{sequence};
  std::queue<std::pair<int64_t, int>> sequences_to_mutate;
  sequences_to_mutate.push({sequence, 0});
  
  while (!sequences_to_mutate.empty()) {
    auto current = sequences_to_mutate.front();
    int64_t current_seq = current.first;
    int current_round = current.second;
    sequences_to_mutate.pop();
    
    if (current_round < mutation_rounds) {
      std::vector<int64_t> new_mutations = generate_quaternary_mutations_cpp(current_seq, sequence_length);
      for (const auto& mutation : new_mutations) {
        if (all_mutations.insert(mutation).second) {
          sequences_to_mutate.push(std::make_pair(mutation, current_round + 1));
        }
      }
    }
  }
  return std::vector<int64_t>(all_mutations.begin(), all_mutations.end());
}

std::vector<int64_t> generate_recursive_quaternary_mutations_cpp_v3(
    int64_t sequence,
    int mutation_rounds,
    int sequence_length = 32) {
  // Use unique_ptr to allocate containers on the heap
  auto all_mutations = std::make_unique<std::unordered_set<int64_t>>();
  all_mutations->insert(sequence);
  
  auto sequences_to_mutate = std::make_unique<std::queue<std::pair<int64_t, int>>>();
  sequences_to_mutate->push({sequence, 0});
  
  while (!sequences_to_mutate->empty()) {
    auto current = sequences_to_mutate->front();
    int64_t current_seq = current.first;
    int current_round = current.second;
    sequences_to_mutate->pop();
    
    if (current_round < mutation_rounds) {
      // Store new mutations on the heap
      auto new_mutations = std::make_unique<std::vector<int64_t>>(generate_quaternary_mutations_cpp(current_seq, sequence_length));
      for (const auto& mutation : *new_mutations) {
        if (all_mutations->insert(mutation).second) {
          sequences_to_mutate->push(std::make_pair(mutation, current_round + 1));
        }
      }
    }
  }
  
  // Return vector created from the heap-allocated unordered_set
  return std::vector<int64_t>(all_mutations->begin(), all_mutations->end());
}

// [[Rcpp::export]]
Rcpp::NumericVector generate_recursive_quaternary_mutations(
    Rcpp::NumericVector sequences,
    int mutation_rounds = 2, 
  Rcpp::IntegerVector sequence_length = 32) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  std::vector<int64_t> results;
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0;
  
  for (size_t i = 0; i < cpp_input_list.size(); ++i) {
    int64_t quaternary_input = cpp_input_list[i];
    int seqlength = use_sl ? seq_length : sequence_length[i];
    std::vector<int64_t> subresult = generate_recursive_quaternary_mutations_cpp(quaternary_input, mutation_rounds, seqlength); 
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}

// Helper function to print binary representation of int64_t
std::string int64_to_binary(int64_t value, int length) {
  std::string binary = std::bitset<64>(value).to_string();
  return binary.substr(64 - 2*length);
}

// Helper function to convert int64_t to quaternary string
std::string int64_to_quaternary(int64_t value, int length) {
  std::string result;
  for (int i = length - 1; i >= 0; --i) {
    int base = (value >> (2 * i)) & 0b11;
    switch(base) {
    case 0: result += 'A'; break;
    case 1: result += 'C'; break;
    case 2: result += 'T'; break;
    case 3: result += 'G'; break;
    }
  }
  return result;
}

std::vector<int64_t> generate_shifted_barcodes(int64_t original, int sequence_length, int max_shift = 2) {
  std::unordered_set<int64_t> shifted_barcodes_set;  // Use a set to avoid duplicates
  int64_t mask = (1LL << (2 * sequence_length)) - 1;  // Mask for the original sequence length
  
  // Generate right shifts
  for (int shift = 1; shift <= max_shift; ++shift) {
    for (int padding = 0; padding < 4; ++padding) {
      int64_t shifted = ((original >> (2 * shift)) & (mask >> (2 * shift))) | (padding << (2 * (sequence_length - shift)));
      shifted_barcodes_set.insert(shifted);  // Insert into set
    }
  }
  
  // Generate left shifts
  for (int shift = 1; shift <= max_shift; ++shift) {
    for (int padding = 0; padding < 4; ++padding) {
      int64_t shifted = ((original << (2 * shift)) & mask) | padding;
      shifted_barcodes_set.insert(shifted);  // Insert into set
    }
  }
  
  // Convert set to vector and return
  return std::vector<int64_t>(shifted_barcodes_set.begin(), shifted_barcodes_set.end());
}

// [[Rcpp::export]]
Rcpp::List shift_barcodes(
    SEXP sequences,
    int sequence_length,
    int max_shift = 2,
    bool verbose = false
  ) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  Rcpp::List result(cpp_input_list.size());
  for (size_t seq_index = 0; seq_index < cpp_input_list.size(); ++seq_index) {
    int64_t original = cpp_input_list[seq_index];
    if (verbose) {
      Rcpp::Rcout << "Processing sequence " << seq_index + 1 << ":\n";
      Rcpp::Rcout << "Original sequence: " << int64_to_quaternary(original, sequence_length) << "\n";
      Rcpp::Rcout << "Binary representation: " << int64_to_binary(original, sequence_length) << "\n";
    }
    std::vector<int64_t> shifted = generate_shifted_barcodes(original, sequence_length, max_shift);
    Rcpp::List seq_result(shifted.size());  // +1 to include the original sequence
    //seq_result[0] = Rcpp::wrap(original);  // Include the original sequence
    for (size_t i = 0; i < shifted.size(); ++i) {
      //seq_result[i + 1] = Rcpp::wrap(shifted[i]);
      seq_result[i] = Rcpp::wrap(shifted[i]);
      if (verbose) {
        Rcpp::Rcout << "Shifted " << (i < shifted.size()/2 ? "right" : "left") << ": " 
                    << int64_to_quaternary(shifted[i], sequence_length) << "\n";
        Rcpp::Rcout << "Binary: " << int64_to_binary(shifted[i], sequence_length) << "\n";
      }
    }
    result[seq_index] = seq_result;
  }
  return result;
}

// Helper function to convert int64_t to string representation
std::string int64_to_string(uint64_t value, int length) {
  std::string result;
  for (int i = length - 1; i >= 0; --i) {
    int base = (value >> (2 * i)) & 0b11;
    switch(base) {
    case 0: result += 'A'; break;
    case 1: result += 'C'; break;
    case 2: result += 'T'; break;
    case 3: result += 'G'; break;
    }
  }
  return result;
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List generate_and_filter_mutations_v3(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    std::string input_fastq = "",
    std::string output_fastq = "",
    bool verbose = false,
    bool generate_r_output = true,
    bool update_fastq = true,
    std::string barcode_header = "barcode:") {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  std::vector<int64_t> seq_vec = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  if (seq_vec.empty()) {
    Rcpp::stop("Input sequences vector is empty.");
  }
  
  Rcpp::Rcout << "Number of input correct sequences: " << seq_vec.size() << "\n";
  std::vector<int64_t> invalid_vec = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  Rcpp::Rcout << "Number of barcodes to correct: " << invalid_vec.size() << "\n";
  
  std::unordered_set<int64_t> invalid_set(invalid_vec.begin(), invalid_vec.end());
  if(verbose) {
    Rcpp::Rcout << "Size of invalid_set: " << invalid_set.size() << "\n";
  }
  
  std::unordered_set<int64_t> all_unique_results;
  std::vector<std::pair<int64_t, int64_t>> shifted_results, mutated_results;
  std::unordered_map<int64_t, std::vector<int64_t>> collision_map;
  size_t collision_count = 0;
  
  Rcpp::Rcout << "Generating shifted and mutated barcodes...\n";
  
  size_t total_sequences = seq_vec.size();
  size_t update_frequency = std::max(static_cast<size_t>(1), total_sequences / 100);
  size_t next_update = update_frequency;
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic) reduction(+:collision_count)
  for (size_t i = 0; i < total_sequences; ++i) {
    int64_t current_seq = seq_vec[i];
    // Generate and filter shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(current_seq, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (invalid_set.find(shifted) != invalid_set.end()) {
#pragma omp critical
{
  shifted_results.emplace_back(shifted, current_seq);
  collision_map[shifted].push_back(current_seq);
  if (collision_map[shifted].size() > 1) {
    collision_count++;
  }
}
      }
    }
    // Generate and filter mutations
  std::vector<int64_t> mutations;
  //#pragma omp critical
  //{
    mutations = generate_recursive_quaternary_mutations_cpp_v2(current_seq, mutation_rounds, sequence_length);
  //}
    for (const auto& mutation : mutations) {
      if (invalid_set.find(mutation) != invalid_set.end()) {
#pragma omp critical
{
  mutated_results.emplace_back(mutation, current_seq);
  collision_map[mutation].push_back(current_seq);
  if (collision_map[mutation].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
#pragma omp critical
{
  if (i + 1 >= next_update || i == total_sequences - 1) {
    double progress = (i + 1.0) / total_sequences * 100.0;
    std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% complete" << std::flush;
    next_update = std::min(next_update + update_frequency, total_sequences);
  }
}
  }
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "Resolving collisions...\n";
  
  std::vector<std::pair<int64_t, std::vector<int64_t>>> resolved_collisions;
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  double p = 0.0, bt = 0.0;
  int q = 1, useBytes = 0;
  int method = 2; // Damerau-Levenshtein distance
  
  for (const auto& entry : collision_map) {
    if (entry.second.size() > 1) {
      std::string incorrect_barcode = int64_to_string(entry.first, sequence_length);
      std::vector<std::string> putative_correct_barcodes;
      for (const auto& match : entry.second) {
        putative_correct_barcodes.push_back(int64_to_string(match, sequence_length));
      }
      
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_barcode);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_barcodes);
      
      SEXP dl_dist_results = sd_stringdist(incorrect_sexp, 
        correct_sexp,
        Rcpp::wrap(method), 
        weight,
        Rcpp::wrap(p), 
        Rcpp::wrap(bt), 
        Rcpp::wrap(q), 
        Rcpp::wrap(useBytes),
        Rcpp::wrap(nthread));
      
      std::vector<double> dl_dist = Rcpp::as<std::vector<double>>(dl_dist_results);
      
      double min_distance = *std::min_element(dl_dist.begin(), dl_dist.end());
      std::vector<int> min_indices;
      
      for (int i = 0; i < dl_dist.size(); ++i) {
        if (dl_dist[i] == min_distance) {
          min_indices.push_back(i);
        }
      }
      
      if (min_indices.size() == 1) {
        // Single minimum
        resolved_collisions.emplace_back(entry.first, std::vector<int64_t>{entry.second[min_indices[0]]});
      } else {
        // Tie
        std::vector<int64_t> tied_barcodes;
        for (int index : min_indices) {
          tied_barcodes.push_back(entry.second[index]);
        }
        resolved_collisions.emplace_back(entry.first, tied_barcodes);
      }
    }
  }
  
  Rcpp::Rcout << "Shifted results: " << shifted_results.size() << "\n";
  Rcpp::Rcout << "Mutated results: " << mutated_results.size() << "\n";
  Rcpp::Rcout << "Resolved collisions: " << resolved_collisions.size() << "\n";
  
  // Process FASTQ file if provided
  if (!input_fastq.empty() && update_fastq) {
    std::cout << "Processing FASTQ file...\n";
    // Open input gzipped FASTQ file for reading
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    // Open output gzipped FASTQ file for writing
    gzFile gz_out = gzopen(output_fastq.c_str(), "wb");
    if (gz_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ file.");
    }
    
    kseq_t* seq = kseq_init(fp);
    
    std::unordered_map<int64_t, int64_t> correction_map;
    
    // Prepare correction map
    for (const auto& result : shifted_results) {
      correction_map[result.first] = result.second;
    }
    for (const auto& result : mutated_results) {
      correction_map[result.first] = result.second;
    }
    for (const auto& result : resolved_collisions) {
      if (result.second.size() == 1) {
        correction_map[result.first] = result.second[0];
      }
    }
    
    size_t processed_reads = 0;
    while (kseq_read(seq) >= 0) {
      std::string name(seq->name.s);
      size_t header_pos = name.find(barcode_header);
      if (header_pos != std::string::npos) {
        size_t barcode_start = header_pos + barcode_header.length();
        size_t barcode_end = name.find_first_of(" \t", barcode_start);
        if (barcode_end == std::string::npos) {
          barcode_end = name.length();
        }
        std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
        
        std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
        if (!barcode_bits.empty()) {
          auto it = correction_map.find(barcode_bits[0]);
          if (it != correction_map.end()) {
            // Replace the incorrect barcode with the corrected one
            std::string corrected_barcode = int64_to_string(it->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else {
            gzprintf(gz_out, "@%s\n", name.c_str());
          }
        } else {
          gzprintf(gz_out, "@%s\n", name.c_str());
        }
      } else {
        gzprintf(gz_out, "@%s\n", name.c_str());
      }
      
      gzprintf(gz_out, "%s\n", seq->seq.s);
      gzprintf(gz_out, "+%s\n", seq->comment.l ? seq->comment.s : "");
      gzprintf(gz_out, "%s\n", seq->qual.s);
      
      processed_reads++;
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << processed_reads << " reads..." << std::flush;
      }
    }
    
    std::cout << "\nFASTQ processing complete. Processed " << processed_reads << " reads.\n";
    
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(gz_out);  // Close the output file
    
    std::cout << "FASTQ file updated.\n";
  }
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  std::cout << "Total execution time: " << duration.count() << " seconds\n";
  
  if (generate_r_output) {
    // Prepare results for R
    Rcpp::List r_shifted_results(shifted_results.size());
    Rcpp::List r_mutated_results(mutated_results.size());
    Rcpp::List r_resolved_collisions(resolved_collisions.size());
    
    for (size_t i = 0; i < shifted_results.size(); ++i) {
      r_shifted_results[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(shifted_results[i].first),
        Rcpp::Named("corrected_barcode") = Rcpp::wrap(shifted_results[i].second),
        Rcpp::Named("type") = "shifted"
      );
    }
    
    for (size_t i = 0; i < mutated_results.size(); ++i) {
      r_mutated_results[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(mutated_results[i].first),
        Rcpp::Named("corrected_barcode") = Rcpp::wrap(mutated_results[i].second),
        Rcpp::Named("type") = "mutated"
      );
    }
    
    for (size_t i = 0; i < resolved_collisions.size(); ++i) {
      r_resolved_collisions[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(resolved_collisions[i].first),
        Rcpp::Named("corrected_barcodes") = Rcpp::wrap(resolved_collisions[i].second),
        Rcpp::Named("type") = resolved_collisions[i].second.size() == 1 ? "resolved_single" : "unresolved"
      );
    }
    
    return Rcpp::List::create(
      Rcpp::Named("shifted") = r_shifted_results,
      Rcpp::Named("mutated") = r_mutated_results,
      Rcpp::Named("resolved") = r_resolved_collisions
    );
  } else {
    return Rcpp::List::create(Rcpp::Named("message") = "Processing complete");
  }
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List generate_and_filter_mutations_v4(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    std::string input_fastq = "",
    std::string output_fastq = "",
    bool verbose = false,
    bool generate_r_output = true,
    bool update_fastq = true,
    std::string barcode_header = "barcode:",
    std::string counts_output_csv = "barcode_counts.csv") {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> seq_vec = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  if (seq_vec.empty()) {
    Rcpp::stop("Input sequences vector is empty.");
  }
  
  Rcpp::Rcout << "Number of input correct sequences: " << seq_vec.size() << "\n";
  std::vector<int64_t> invalid_vec = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  Rcpp::Rcout << "Number of barcodes to correct: " << invalid_vec.size() << "\n";
  
  std::unordered_set<int64_t> invalid_set(invalid_vec.begin(), invalid_vec.end());
  if(verbose) {
    Rcpp::Rcout << "Size of invalid_set: " << invalid_set.size() << "\n";
  }
  
  // Step 2: Initialize data structures to store results and track counts
  std::unordered_set<int64_t> all_unique_results;
  std::vector<std::pair<int64_t, int64_t>> shifted_results, mutated_results;
  std::unordered_map<int64_t, std::vector<int64_t>> collision_map;
  std::unordered_map<int64_t, int> original_counts, corrected_counts;  // Track original and corrected counts
  size_t collision_count = 0;
  
  Rcpp::Rcout << "Generating shifted and mutated barcodes...\n";
  
  size_t total_sequences = seq_vec.size();
  size_t update_frequency = std::max(static_cast<size_t>(1), total_sequences / 100);
  size_t next_update = update_frequency;
  
  // Step 3: OpenMP parallel loop for barcode correction
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic) reduction(+:collision_count)
  for (size_t i = 0; i < total_sequences; ++i) {
    int64_t current_seq = seq_vec[i];
    original_counts[current_seq]++;  // Track original counts
    
    // Generate and filter shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(current_seq, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (invalid_set.find(shifted) != invalid_set.end()) {
#pragma omp critical
{
  shifted_results.emplace_back(shifted, current_seq);
  collision_map[shifted].push_back(current_seq);
  if (collision_map[shifted].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
    // Generate and filter mutations
    std::vector<int64_t> mutations = generate_recursive_quaternary_mutations_cpp_v2(current_seq, mutation_rounds, sequence_length);
    for (const auto& mutation : mutations) {
      if (invalid_set.find(mutation) != invalid_set.end()) {
#pragma omp critical
{
  mutated_results.emplace_back(mutation, current_seq);
  collision_map[mutation].push_back(current_seq);
  if (collision_map[mutation].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
#pragma omp critical
{
  if (i + 1 >= next_update || i == total_sequences - 1) {
    double progress = (i + 1.0) / total_sequences * 100.0;
    std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% complete" << std::flush;
    next_update = std::min(next_update + update_frequency, total_sequences);
  }
}
  }
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "Resolving collisions...\n";
  
  // Step 4: Collision resolution using Damerau-Levenshtein distance
  std::vector<std::pair<int64_t, std::vector<int64_t>>> resolved_collisions;
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  double p = 0.0, bt = 0.0;
  int q = 1, useBytes = 0;
  int method = 2;  // Damerau-Levenshtein distance
  
  for (const auto& entry : collision_map) {
    if (entry.second.size() > 1) {
      std::string incorrect_barcode = int64_to_string(entry.first, sequence_length);
      std::vector<std::string> putative_correct_barcodes;
      for (const auto& match : entry.second) {
        putative_correct_barcodes.push_back(int64_to_string(match, sequence_length));
      }
      
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_barcode);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_barcodes);
      
      SEXP dl_dist_results = sd_stringdist(incorrect_sexp, correct_sexp, Rcpp::wrap(method), weight, Rcpp::wrap(p), Rcpp::wrap(bt), Rcpp::wrap(q), Rcpp::wrap(useBytes), Rcpp::wrap(nthread));
      
      std::vector<double> dl_dist = Rcpp::as<std::vector<double>>(dl_dist_results);
      double min_distance = *std::min_element(dl_dist.begin(), dl_dist.end());
      std::vector<int> min_indices;
      
      for (int i = 0; i < dl_dist.size(); ++i) {
        if (dl_dist[i] == min_distance) {
          min_indices.push_back(i);
        }
      }
      
      if (min_indices.size() == 1) {
        // Single minimum
        resolved_collisions.emplace_back(entry.first, std::vector<int64_t>{entry.second[min_indices[0]]});
      } else {
        // Tie
        std::vector<int64_t> tied_barcodes;
        for (int index : min_indices) {
          tied_barcodes.push_back(entry.second[index]);
        }
        resolved_collisions.emplace_back(entry.first, tied_barcodes);
      }
    }
  }
  
  // Step 5: Prepare correction map for fastq updates
  std::unordered_map<int64_t, int64_t> correction_map;
  for (const auto& result : shifted_results) {
    correction_map[result.first] = result.second;
    corrected_counts[result.second]++;  // Update corrected counts
  }
  for (const auto& result : mutated_results) {
    correction_map[result.first] = result.second;
    corrected_counts[result.second]++;
  }
  for (const auto& result : resolved_collisions) {
    if (result.second.size() == 1) {
      correction_map[result.first] = result.second[0];
      corrected_counts[result.second[0]]++;
    }
  }
  
  // Step 6: Process FASTQ file and update barcodes if required
  if (!input_fastq.empty() && update_fastq) {
    std::cout << "Processing FASTQ file...\n";
    // Open input gzipped FASTQ file for reading
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    // Open output gzipped FASTQ file for writing
    gzFile gz_out = gzopen(output_fastq.c_str(), "wb");
    if (gz_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ file.");
    }
    
    kseq_t* seq = kseq_init(fp);
    size_t processed_reads = 0;
    while (kseq_read(seq) >= 0) {
      std::string name(seq->name.s);
      size_t header_pos = name.find(barcode_header);
      if (header_pos != std::string::npos) {
        size_t barcode_start = header_pos + barcode_header.length();
        size_t barcode_end = name.find_first_of(" \t", barcode_start);
        if (barcode_end == std::string::npos) {
          barcode_end = name.length();
        }
        std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
        std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
        if (!barcode_bits.empty()) {
          auto it = correction_map.find(barcode_bits[0]);
          if (it != correction_map.end()) {
            // Replace the incorrect barcode with the corrected one
            std::string corrected_barcode = int64_to_string(it->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else {
            gzprintf(gz_out, "@%s\n", name.c_str());
          }
        } else {
          gzprintf(gz_out, "@%s\n", name.c_str());
        }
      } else {
        gzprintf(gz_out, "@%s\n", name.c_str());
      }
      
      gzprintf(gz_out, "%s\n", seq->seq.s);
      gzprintf(gz_out, "+%s\n", seq->comment.l ? seq->comment.s : "");
      gzprintf(gz_out, "%s\n", seq->qual.s);
      
      processed_reads++;
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << processed_reads << " reads..." << std::flush;
      }
    }
    
    std::cout << "\nFASTQ processing complete. Processed " << processed_reads << " reads.\n";
    
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(gz_out);  // Close the output file
    
    std::cout << "FASTQ file updated.\n";
  }
  
  // Step 7: Write barcode counts to CSV
  std::ofstream csv_file(counts_output_csv);
  csv_file << "barcode,original_counts,corrected_counts,total_counts\n";
  for (const auto& [barcode, orig_count] : original_counts) {
    int corrected_count = corrected_counts[barcode];
    csv_file << int64_to_string(barcode, sequence_length) << ","
             << orig_count << ","
             << corrected_count << ","
             << (orig_count + corrected_count) << "\n";
  }
  csv_file.close();
  
  // Step 8: Generate R output if requested
  if (generate_r_output) {
    // Prepare results for R
    Rcpp::List r_shifted_results(shifted_results.size());
    Rcpp::List r_mutated_results(mutated_results.size());
    Rcpp::List r_resolved_collisions(resolved_collisions.size());
    
    for (size_t i = 0; i < shifted_results.size(); ++i) {
      r_shifted_results[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(shifted_results[i].first),
        Rcpp::Named("corrected_barcode") = Rcpp::wrap(shifted_results[i].second),
        Rcpp::Named("type") = "shifted"
      );
    }
    
    for (size_t i = 0; i < mutated_results.size(); ++i) {
      r_mutated_results[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(mutated_results[i].first),
        Rcpp::Named("corrected_barcode") = Rcpp::wrap(mutated_results[i].second),
        Rcpp::Named("type") = "mutated"
      );
    }
    
    for (size_t i = 0; i < resolved_collisions.size(); ++i) {
      r_resolved_collisions[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(resolved_collisions[i].first),
        Rcpp::Named("corrected_barcodes") = Rcpp::wrap(resolved_collisions[i].second),
        Rcpp::Named("type") = resolved_collisions[i].second.size() == 1 ? "resolved_single" : "unresolved"
      );
    }
    
    return Rcpp::List::create(
      Rcpp::Named("shifted") = r_shifted_results,
      Rcpp::Named("mutated") = r_mutated_results,
      Rcpp::Named("resolved") = r_resolved_collisions
    );
  } else {
    return Rcpp::List::create(Rcpp::Named("message") = "Processing complete");
  }
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List generate_and_filter_mutations_v5(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    std::string input_fastq = "",
    std::string output_fastq = "",
    bool verbose = false,
    bool generate_r_output = true,
    bool update_fastq = true,
    std::string barcode_header = "barcode:",
    std::string counts_output_csv = "barcode_counts.csv") {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> seq_vec = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int> original_counts_input = Rcpp::as<std::vector<int>>(true_counts);
  Rcpp::Rcout << "Number of input correct sequences: " << seq_vec.size() << "\n";
  
  if (seq_vec.empty()) {
    Rcpp::stop("Input sequences vector is empty.");
  }
  
  std::vector<int64_t> invalid_vec = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::vector<int> invalid_counts_input = Rcpp::as<std::vector<int>>(invalid_counts);
  Rcpp::Rcout << "Number of barcodes to correct: " << invalid_vec.size() << "\n";
  
  std::unordered_set<int64_t> invalid_set(invalid_vec.begin(), invalid_vec.end());
  if(verbose) {
    Rcpp::Rcout << "Size of invalid_set: " << invalid_set.size() << "\n";
  }
  
  // Step 2: Initialize data structures to store results and track counts
  std::vector<std::pair<int64_t, int64_t>> shifted_results, mutated_results;
  std::unordered_map<int64_t, std::vector<int64_t>> collision_map;
  std::unordered_map<int64_t, int64_t> correction_map;
  
  std::unordered_map<int64_t, int> original_counts, corrected_counts;
  for (size_t i = 0; i < seq_vec.size(); ++i) {
    original_counts[seq_vec[i]] = original_counts_input[i];
  }
  for (size_t i = 0; i < invalid_vec.size(); ++i) {
    corrected_counts[invalid_vec[i]] = 0;  // Initialize corrected counts for invalid barcodes
  }
  
  size_t collision_count = 0;
  
  Rcpp::Rcout << "Generating shifted and mutated barcodes...\n";
  
  size_t total_sequences = seq_vec.size();
  size_t update_frequency = std::max(static_cast<size_t>(1), total_sequences / 100);
  size_t next_update = update_frequency;
  
  // Step 3: OpenMP parallel loop for barcode correction
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic) reduction(+:collision_count)
  for (size_t i = 0; i < total_sequences; ++i) {
    int64_t current_seq = seq_vec[i];
    
    // Generate and filter shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(current_seq, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (invalid_set.find(shifted) != invalid_set.end()) {
#pragma omp critical
{
  shifted_results.emplace_back(shifted, current_seq);
  collision_map[shifted].push_back(current_seq);
  if (collision_map[shifted].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
    // Generate and filter mutations
    std::vector<int64_t> mutations = generate_recursive_quaternary_mutations_cpp_v2(current_seq, mutation_rounds, sequence_length);
    for (const auto& mutation : mutations) {
      if (invalid_set.find(mutation) != invalid_set.end()) {
#pragma omp critical
{
  mutated_results.emplace_back(mutation, current_seq);
  collision_map[mutation].push_back(current_seq);
  if (collision_map[mutation].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
#pragma omp critical
{
  if (i + 1 >= next_update || i == total_sequences - 1) {
    double progress = (i + 1.0) / total_sequences * 100.0;
    std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% complete" << std::flush;
    next_update = std::min(next_update + update_frequency, total_sequences);
  }
}
  }
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "Resolving collisions...\n";
  
  // Step 4: Collision resolution using Damerau-Levenshtein distance
  std::vector<std::pair<int64_t, std::vector<int64_t>>> resolved_collisions;
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  double p = 0.0, bt = 0.0;
  int q = 1, useBytes = 0;
  int method = 2;  // Damerau-Levenshtein distance
  
  for (const auto& entry : collision_map) {
    if (entry.second.size() > 1) {
      std::string incorrect_barcode = int64_to_string(entry.first, sequence_length);
      std::vector<std::string> putative_correct_barcodes;
      for (const auto& match : entry.second) {
        putative_correct_barcodes.push_back(int64_to_string(match, sequence_length));
      }
      
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_barcode);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_barcodes);
      
      SEXP dl_dist_results = sd_stringdist(incorrect_sexp, correct_sexp, Rcpp::wrap(method), weight, Rcpp::wrap(p), Rcpp::wrap(bt), Rcpp::wrap(q), Rcpp::wrap(useBytes), Rcpp::wrap(nthread));
      
      std::vector<double> dl_dist = Rcpp::as<std::vector<double>>(dl_dist_results);
      double min_distance = *std::min_element(dl_dist.begin(), dl_dist.end());
      std::vector<int> min_indices;
      
      for (int i = 0; i < dl_dist.size(); ++i) {
        if (dl_dist[i] == min_distance) {
          min_indices.push_back(i);
        }
      }
      
      if (min_indices.size() == 1) {
        // Single minimum
        resolved_collisions.emplace_back(entry.first, std::vector<int64_t>{entry.second[min_indices[0]]});
      } else {
        // Tie
        std::vector<int64_t> tied_barcodes;
        for (int index : min_indices) {
          tied_barcodes.push_back(entry.second[index]);
        }
        resolved_collisions.emplace_back(entry.first, tied_barcodes);
      }
    }
  }
  
  // Step 5: Prepare correction map for fastq updates
  for (const auto& result : shifted_results) {
    correction_map[result.first] = result.second;
    corrected_counts[result.second]++;  // Update corrected counts
  }
  for (const auto& result : mutated_results) {
    correction_map[result.first] = result.second;
    corrected_counts[result.second]++;
  }
  for (const auto& result : resolved_collisions) {
    if (result.second.size() == 1) {
      correction_map[result.first] = result.second[0];
      corrected_counts[result.second[0]]++;
    }
  }
  
  // Step 6: Process FASTQ file and update barcodes if required
  if (!input_fastq.empty() && update_fastq) {
    std::cout << "Processing FASTQ file...\n";
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    gzFile gz_out = gzopen(output_fastq.c_str(), "wb");
    if (gz_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ file.");
    }
    
    kseq_t* seq = kseq_init(fp);
    size_t processed_reads = 0;
    while (kseq_read(seq) >= 0) {
      std::string name(seq->name.s);
      size_t header_pos = name.find(barcode_header);
      if (header_pos != std::string::npos) {
        size_t barcode_start = header_pos + barcode_header.length();
        size_t barcode_end = name.find_first_of(" \t", barcode_start);
        if (barcode_end == std::string::npos) {
          barcode_end = name.length();
        }
        std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
        std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
        if (!barcode_bits.empty()) {
          auto it = correction_map.find(barcode_bits[0]);
          if (it != correction_map.end()) {
            std::string corrected_barcode = int64_to_string(it->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else {
            gzprintf(gz_out, "@%s\n", name.c_str());
          }
        } else {
          gzprintf(gz_out, "@%s\n", name.c_str());
        }
      } else {
        gzprintf(gz_out, "@%s\n", name.c_str());
      }
      
      gzprintf(gz_out, "%s\n", seq->seq.s);
      gzprintf(gz_out, "+%s\n", seq->comment.l ? seq->comment.s : "");
      gzprintf(gz_out, "%s\n", seq->qual.s);
      
      processed_reads++;
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << processed_reads << " reads..." << std::flush;
      }
    }
    std::cout << "\nFASTQ processing complete. Processed " << processed_reads << " reads.\n";
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(gz_out);  // Close the output file
    std::cout << "FASTQ file updated.\n";
  }
  // Step 7: Write barcode counts to CSV
  std::ofstream csv_file(counts_output_csv);
  csv_file << "barcode,original_counts,corrected_counts,total_counts\n";
  for (const auto& [barcode, orig_count] : original_counts) {
    int corrected_count = corrected_counts[barcode];
    csv_file << int64_to_string(barcode, sequence_length) << ","
             << orig_count << ","
             << corrected_count << ","
             << (orig_count + corrected_count) << "\n";
  }
  csv_file.close();
  
  // Step 8: Generate R output if requested
  if (generate_r_output) {
    Rcpp::List r_shifted_results(shifted_results.size());
    Rcpp::List r_mutated_results(mutated_results.size());
    Rcpp::List r_resolved_collisions(resolved_collisions.size());
    
    for (size_t i = 0; i < shifted_results.size(); ++i) {
      r_shifted_results[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(shifted_results[i].first),
        Rcpp::Named("corrected_barcode") = Rcpp::wrap(shifted_results[i].second),
        Rcpp::Named("type") = "shifted"
      );
    }
    
    for (size_t i = 0; i < mutated_results.size(); ++i) {
      r_mutated_results[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(mutated_results[i].first),
        Rcpp::Named("corrected_barcode") = Rcpp::wrap(mutated_results[i].second),
        Rcpp::Named("type") = "mutated"
      );
    }
    
    for (size_t i = 0; i < resolved_collisions.size(); ++i) {
      r_resolved_collisions[i] = Rcpp::List::create(
        Rcpp::Named("incorrect_barcode") = Rcpp::wrap(resolved_collisions[i].first),
        Rcpp::Named("corrected_barcodes") = Rcpp::wrap(resolved_collisions[i].second),
        Rcpp::Named("type") = resolved_collisions[i].second.size() == 1 ? "resolved_single" : "unresolved"
      );
    }
    
    return Rcpp::List::create(
      Rcpp::Named("shifted") = r_shifted_results,
      Rcpp::Named("mutated") = r_mutated_results,
      Rcpp::Named("resolved") = r_resolved_collisions
    );
  } else {
    return Rcpp::List::create(Rcpp::Named("message") = "Processing complete");
  }
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List generate_and_filter_mutations_v6(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    std::string input_fastq = "",
    std::string filtered_fastq = "filtered_fastq.gz",
    std::string unfiltered_fastq = "unfiltered_fastq.gz",
    bool verbose = false,
    bool verbose_output = true,
    std::string barcode_header = "barcode:",
    std::string detailed_output_csv = "",
    std::string counts_output_csv = "") {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> seq_vec = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int> original_counts_input = Rcpp::as<std::vector<int>>(true_counts);
  Rcpp::Rcout << "Number of input correct sequences: " << seq_vec.size() << "\n";
  
  if (seq_vec.empty()) {
    Rcpp::stop("Input sequences vector is empty.");
  }
  
  std::vector<int64_t> invalid_vec = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::vector<int> invalid_counts_input = Rcpp::as<std::vector<int>>(invalid_counts);
  Rcpp::Rcout << "Number of barcodes to correct: " << invalid_vec.size() << "\n";
  
  std::unordered_set<int64_t> invalid_set(invalid_vec.begin(), invalid_vec.end());
  if(verbose) {
    Rcpp::Rcout << "Size of invalid_set: " << invalid_set.size() << "\n";
  }
  
  // Step 2: Initialize data structures to store results and track counts
  std::vector<std::pair<int64_t, int64_t>> shifted_results, mutated_results;
  std::unordered_map<int64_t, std::vector<int64_t>> collision_map;
  std::unordered_map<int64_t, int64_t> correction_map;
  
  std::unordered_map<int64_t, int> original_counts, corrected_counts;
  for (size_t i = 0; i < seq_vec.size(); ++i) {
    original_counts[seq_vec[i]] = original_counts_input[i];
  }
  for (size_t i = 0; i < invalid_vec.size(); ++i) {
    corrected_counts[invalid_vec[i]] = 0;  // Initialize corrected counts for invalid barcodes
  }
  
  size_t collision_count = 0;
  Rcpp::Rcout << "Generating shifted and mutated barcodes...\n";
  
  size_t total_sequences = seq_vec.size();
  size_t update_frequency = std::max(static_cast<size_t>(1), total_sequences / 100);
  size_t next_update = update_frequency;
  
  // Step 3: barcode correction
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic) reduction(+:collision_count)
  for (size_t i = 0; i < total_sequences; ++i) {
    int64_t current_seq = seq_vec[i];
    
    // Generate and filter shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(current_seq, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (invalid_set.find(shifted) != invalid_set.end()) {
#pragma omp critical
{
  shifted_results.emplace_back(shifted, current_seq);
  collision_map[shifted].push_back(current_seq);
  if (collision_map[shifted].size() > 1) {
    collision_count++;
  }
}
      }
    }
    // Generate and filter mutations
    std::vector<int64_t> mutations = generate_recursive_quaternary_mutations_cpp_v2(current_seq, mutation_rounds, sequence_length);
    for (const auto& mutation : mutations) {
      if (invalid_set.find(mutation) != invalid_set.end()) {
#pragma omp critical
{
  mutated_results.emplace_back(mutation, current_seq);
  collision_map[mutation].push_back(current_seq);
  if (collision_map[mutation].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
#pragma omp critical
{
  if (i + 1 >= next_update || i == total_sequences - 1) {
    double progress = (i + 1.0) / total_sequences * 100.0;
    std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% complete" << std::flush;
    next_update = std::min(next_update + update_frequency, total_sequences);
  }
}
  }
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "Resolving collisions...\n";
  
  // Step 4: Collision resolution using Damerau-Levenshtein distance
  std::vector<std::pair<int64_t, std::vector<int64_t>>> resolved_collisions;
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  double p = 0.0, bt = 0.0;
  int q = 1, useBytes = 0;
  int method = 2;  // Damerau-Levenshtein distance
  
  for (const auto& entry : collision_map) {
    if (entry.second.size() > 1) {
      std::string incorrect_barcode = int64_to_string(entry.first, sequence_length);
      std::vector<std::string> putative_correct_barcodes;
      for (const auto& match : entry.second) {
        putative_correct_barcodes.push_back(int64_to_string(match, sequence_length));
      }
      
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_barcode);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_barcodes);
      
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(p),
        Rcpp::wrap(bt),
        Rcpp::wrap(q),
        Rcpp::wrap(useBytes),
        Rcpp::wrap(nthread));
      
      std::vector<double> dl_dist = Rcpp::as<std::vector<double>>(dl_dist_results);
      double min_distance = *std::min_element(dl_dist.begin(), dl_dist.end());
      std::vector<int> min_indices;
      
      for (int i = 0; i < dl_dist.size(); ++i) {
        if (dl_dist[i] == min_distance) {
          min_indices.push_back(i);
        }
      }
      
      if (min_indices.size() == 1) {
        // Single minimum
        resolved_collisions.emplace_back(entry.first, std::vector<int64_t>{entry.second[min_indices[0]]});
      } else {
        // Tie
        std::vector<int64_t> tied_barcodes;
        for (int index : min_indices) {
          tied_barcodes.push_back(entry.second[index]);
        }
        resolved_collisions.emplace_back(entry.first, tied_barcodes);
      }
    }
  }
  
  // Step 5: Prepare correction map for fastq updates
  for (const auto& result : shifted_results) {
    correction_map[result.first] = result.second;
    corrected_counts[result.second]++;  // Update corrected counts
  }
  for (const auto& result : mutated_results) {
    correction_map[result.first] = result.second;
    corrected_counts[result.second]++;
  }
  for (const auto& result : resolved_collisions) {
    if (result.second.size() == 1) {
      correction_map[result.first] = result.second[0];
      corrected_counts[result.second[0]]++;
    }
  }
  for (const auto& result : seq_vec) {
    correction_map[result] = result;  // Map correct barcodes to themselves
    corrected_counts[result]++;  // Increment the count for correct barcodes
  }
  
  // Step 6: Stream FASTQ into filtered and unfiltered FASTQ files
  if (!input_fastq.empty()) {
    std::cout << "\nProcessing and splitting FASTQ file...\n";
    // Open input gzipped FASTQ file for reading
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    // Open filtered and unfiltered gzipped FASTQ files
    gzFile gz_filtered_out = gzopen(filtered_fastq.c_str(), "wb");  // Write mode, overwrite filtered output
    gzFile gz_unfiltered_out = gzopen(unfiltered_fastq.c_str(), "ab");  // Append mode for unfiltered output
    
    if (gz_filtered_out == NULL || gz_unfiltered_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ files.");
    }
    
    kseq_t* seq = kseq_init(fp);
    size_t processed_reads = 0;
    
    // Process each read in the input FASTQ
    while (kseq_read(seq) >= 0) {
      std::string name(seq->name.s);
      size_t header_pos = name.find(barcode_header);
      
      gzFile gz_out = gz_unfiltered_out; // Default to unfiltered output
      
      if (header_pos != std::string::npos) {
        size_t barcode_start = header_pos + barcode_header.length();
        size_t barcode_end = name.find_first_of(" \t|", barcode_start);
        if (barcode_end == std::string::npos) {
          barcode_end = name.length();
        }
        std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
        std::string rc_barcode_str = revcomp_cpp(barcode_str);
        std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
        std::vector<int64_t> revcomp_bits = sequence_to_bits_cpp(rc_barcode_str);
        if (!barcode_bits.empty()) {
          auto it = correction_map.find(barcode_bits[0]);
          auto it_revcomp = correction_map.find(revcomp_bits[0]);
          if (it != correction_map.end()) {
            // Found in correction map (barcode corrected)
            std::string corrected_barcode = int64_to_string(it->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gz_out = gz_filtered_out;  // Switch to filtered output
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else if (it_revcomp != correction_map.end()) {
            // Found reverse complement in correction map
            std::string corrected_barcode = int64_to_string(it_revcomp->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gz_out = gz_filtered_out;  // Switch to filtered output
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else {
            // Barcode not found, write to unfiltered
            gzprintf(gz_out, "@%s\n", name.c_str());
          }
        } else {
          gzprintf(gz_out, "@%s\n", name.c_str());
        }
      } else {
        gzprintf(gz_unfiltered_out, "@%s\n", name.c_str());
      }
      gzprintf(gz_out, "%s\n", seq->seq.s);
      gzprintf(gz_out, "+%s\n", seq->comment.l ? seq->comment.s : "");
      gzprintf(gz_out, "%s\n", seq->qual.s);
      processed_reads++;
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << processed_reads << " reads..." << std::flush;
      }
    }
    
    std::cout << "\nFASTQ processing complete. Processed " << processed_reads << " reads...\n";
    // Clean up and close files
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(gz_filtered_out);
    gzclose(gz_unfiltered_out);
    std::cout << "FASTQ files updated.\n";
  }
  
  // Step 7: Write barcode counts to CSV
  std::ofstream csv_file(counts_output_csv);
  csv_file << "barcode,original_counts,corrected_counts,total_counts\n";
  for (const auto& [barcode, orig_count] : original_counts) {
    int corrected_count = corrected_counts[barcode];
    csv_file << int64_to_string(barcode, sequence_length) << ","
             << orig_count << ","
             << corrected_count << ","
             << (orig_count + corrected_count) << "\n";
  }
  csv_file.close();
  
  // Step 8: Write detailed output to CSV if verbose_output is enabled
  if (verbose_output) {
    std::ofstream detailed_csv(detailed_output_csv);
    detailed_csv << "barcode,type,original,corrected\n";
    
    for (const auto& result : shifted_results) {
      detailed_csv << int64_to_string(result.first, sequence_length) << ",shifted,"
                   << int64_to_string(result.first, sequence_length) << ","
                   << int64_to_string(result.second, sequence_length) << "\n";
    }
    
    for (const auto& result : mutated_results) {
      detailed_csv << int64_to_string(result.first, sequence_length) << ",mutated,"
                   << int64_to_string(result.first, sequence_length) << ","
                   << int64_to_string(result.second, sequence_length) << "\n";
    }
    for (const auto& result : resolved_collisions) {
      if (result.second.size() == 1) {
        detailed_csv << int64_to_string(result.first, sequence_length) << ",resolved_single,"
                     << int64_to_string(result.first, sequence_length) << ","
                     << int64_to_string(result.second[0], sequence_length) << "\n";
      } else {
        // Unresolved barcodes: concatenate them with their distances
        std::string unresolved_str;
        for (size_t i = 0; i < result.second.size(); ++i) {
          std::string unresolved_barcode = int64_to_string(result.second[i], sequence_length);
          unresolved_str += unresolved_barcode;
          if (i != result.second.size() - 1) {
            unresolved_str += "|"; // Use '|' as the delimiter
          }
        }
        detailed_csv << int64_to_string(result.first, sequence_length) << ",unresolved,"
                     << int64_to_string(result.first, sequence_length) << ","
                     << unresolved_str << "\n";
      }
    }
    detailed_csv.close();
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  std::cout << "Total execution time: " << duration.count() << " seconds\n";
  return Rcpp::List::create(Rcpp::Named("message") = "Processing complete");
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List generate_and_filter_mutations_v7(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    std::string input_fastq = "",
    std::string filtered_fastq = "filtered_fastq.gz",
    std::string unfiltered_fastq = "unfiltered_fastq.gz",
    bool verbose = false,
    bool verbose_output = true,
    std::string barcode_header = "barcode:",
    std::string detailed_output_csv = "",
    std::string counts_output_csv = "") {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> seq_vec = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int> original_counts_input = Rcpp::as<std::vector<int>>(true_counts);
  Rcpp::Rcout << "Number of input correct sequences: " << seq_vec.size() << "\n";
  
  if (seq_vec.empty()) {
    Rcpp::stop("Input sequences vector is empty.");
  }
  
  std::vector<int64_t> invalid_vec = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::vector<int> invalid_counts_input = Rcpp::as<std::vector<int>>(invalid_counts);
  Rcpp::Rcout << "Number of barcodes to correct: " << invalid_vec.size() << "\n";
  
  std::unordered_set<int64_t> invalid_set(invalid_vec.begin(), invalid_vec.end());
  if(verbose) {
    Rcpp::Rcout << "Size of invalid_set: " << invalid_set.size() << "\n";
  }
  
  // Step 2: Initialize data structures to store results and track counts
  std::vector<std::pair<int64_t, int64_t>> shifted_results, mutated_results;
  std::unordered_map<int64_t, std::vector<int64_t>> collision_map;
  std::unordered_map<int64_t, int64_t> correction_map;
  
  std::unordered_map<int64_t, int> original_counts, corrected_counts;
  for (size_t i = 0; i < seq_vec.size(); ++i) {
    original_counts[seq_vec[i]] = original_counts_input[i];
  }
  
  size_t collision_count = 0;
  Rcpp::Rcout << "Generating shifted and mutated barcodes...\n";
  
  size_t total_sequences = seq_vec.size();
  size_t update_frequency = std::max(static_cast<size_t>(1), total_sequences / 100);
  size_t next_update = update_frequency;
  
  // Step 3: Barcode correction
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic) reduction(+:collision_count)
  for (size_t i = 0; i < total_sequences; ++i) {
    int64_t current_seq = seq_vec[i];
    
    // Generate and filter shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(current_seq, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (invalid_set.find(shifted) != invalid_set.end()) {
#pragma omp critical
{
  shifted_results.emplace_back(shifted, current_seq);
  collision_map[shifted].push_back(current_seq);
  if (collision_map[shifted].size() > 1) {
    collision_count++;
  }
}
      }
    }
    // Generate and filter mutations
    std::vector<int64_t> mutations = generate_recursive_quaternary_mutations_cpp_v2(current_seq, mutation_rounds, sequence_length);
    for (const auto& mutation : mutations) {
      if (invalid_set.find(mutation) != invalid_set.end()) {
#pragma omp critical
{
  mutated_results.emplace_back(mutation, current_seq);
  collision_map[mutation].push_back(current_seq);
  if (collision_map[mutation].size() > 1) {
    collision_count++;
  }
}
      }
    }
    
#pragma omp critical
{
  if (i + 1 >= next_update || i == total_sequences - 1) {
    double progress = (i + 1.0) / total_sequences * 100.0;
    std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << progress << "% complete" << std::flush;
    next_update = std::min(next_update + update_frequency, total_sequences);
  }
}
  }
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "Resolving collisions...\n";
  
  // Step 4: Collision resolution using Damerau-Levenshtein distance
  std::vector<std::pair<int64_t, std::vector<int64_t>>> resolved_collisions;
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  double p = 0.0, bt = 0.0;
  int q = 1, useBytes = 0;
  int method = 2;  // Damerau-Levenshtein distance
  
  for (const auto& entry : collision_map) {
    if (entry.second.size() > 1) {
      std::string incorrect_barcode = int64_to_string(entry.first, sequence_length);
      std::vector<std::string> putative_correct_barcodes;
      for (const auto& match : entry.second) {
        putative_correct_barcodes.push_back(int64_to_string(match, sequence_length));
      }
      
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_barcode);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_barcodes);
      
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(p),
        Rcpp::wrap(bt),
        Rcpp::wrap(q),
        Rcpp::wrap(useBytes),
        Rcpp::wrap(nthread));
      
      std::vector<double> dl_dist = Rcpp::as<std::vector<double>>(dl_dist_results);
      double min_distance = *std::min_element(dl_dist.begin(), dl_dist.end());
      std::vector<int> min_indices;
      
      for (int i = 0; i < dl_dist.size(); ++i) {
        if (dl_dist[i] == min_distance) {
          min_indices.push_back(i);
        }
      }
      
      if (min_indices.size() == 1) {
        // Single minimum
        resolved_collisions.emplace_back(entry.first, std::vector<int64_t>{entry.second[min_indices[0]]});
      } else {
        // Tie
        std::vector<int64_t> tied_barcodes;
        for (int index : min_indices) {
          tied_barcodes.push_back(entry.second[index]);
        }
        resolved_collisions.emplace_back(entry.first, tied_barcodes);
      }
    }
  }
  
  // Step 5: Prepare correction map
  for (const auto& result : shifted_results) {
    correction_map[result.first] = result.second;
  }
  for (const auto& result : mutated_results) {
    correction_map[result.first] = result.second;
  }
  for (const auto& result : resolved_collisions) {
    if (result.second.size() == 1) {
      correction_map[result.first] = result.second[0];
    }
  }
  for (const auto& result : seq_vec) {
    correction_map[result] = result;  // Map correct barcodes to themselves
  }
  // Initialize read-level counts
  size_t total_reads = 0;
  size_t corrected_reads = 0;
  size_t uncorrected_reads = 0;
  std::unordered_map<int64_t, int> read_counts;  // To count the number of reads per corrected barcode
  
  // Step 6: Stream FASTQ into filtered and unfiltered FASTQ files
  if (!input_fastq.empty()) {
    std::cout << "\nProcessing and splitting FASTQ file...\n";
    // Open input gzipped FASTQ file for reading
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    // Open filtered and unfiltered gzipped FASTQ files
    gzFile gz_filtered_out = gzopen(filtered_fastq.c_str(), "wb");  // Write mode, overwrite filtered output
    gzFile gz_unfiltered_out = gzopen(unfiltered_fastq.c_str(), "wb");  // Write mode for unfiltered output
    
    if (gz_filtered_out == NULL || gz_unfiltered_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ files.");
    }
    
    kseq_t* seq = kseq_init(fp);
    size_t processed_reads = 0;
    
    // Process each read in the input FASTQ
    while (kseq_read(seq) >= 0) {
      total_reads++;  // Increment total reads processed
      std::string name(seq->name.s);
      size_t header_pos = name.find(barcode_header);
      
      gzFile gz_out = gz_unfiltered_out; // Default to unfiltered output
      
      if (header_pos != std::string::npos) {
        size_t barcode_start = header_pos + barcode_header.length();
        size_t barcode_end = name.find_first_of(" \t|", barcode_start);
        if (barcode_end == std::string::npos) {
          barcode_end = name.length();
        }
        std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
        std::string rc_barcode_str = revcomp_cpp(barcode_str);
        std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
        std::vector<int64_t> revcomp_bits = sequence_to_bits_cpp(rc_barcode_str);
        if (!barcode_bits.empty()) {
          auto it = correction_map.find(barcode_bits[0]);
          auto it_revcomp = correction_map.find(revcomp_bits[0]);
          if (it != correction_map.end()) {
            // Found in correction map (barcode corrected)
            corrected_reads++;
            read_counts[it->second]++;  // Increment count for this corrected barcode
            std::string corrected_barcode = int64_to_string(it->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gz_out = gz_filtered_out;  // Switch to filtered output
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else if (it_revcomp != correction_map.end()) {
            // Found reverse complement in correction map
            corrected_reads++;
            read_counts[it_revcomp->second]++;  // Increment count for this corrected barcode
            std::string corrected_barcode = int64_to_string(it_revcomp->second, sequence_length);
            std::string corrected_name = name.substr(0, barcode_start) + corrected_barcode + name.substr(barcode_end);
            gz_out = gz_filtered_out;  // Switch to filtered output
            gzprintf(gz_out, "@%s\n", corrected_name.c_str());
          } else {
            // Barcode not found, write to unfiltered
            uncorrected_reads++;
            gzprintf(gz_out, "@%s\n", name.c_str());
          }
        } else {
          // Invalid barcode format
          uncorrected_reads++;
          gzprintf(gz_out, "@%s\n", name.c_str());
        }
      } else {
        // No barcode in header
        uncorrected_reads++;
        gzprintf(gz_out, "@%s\n", name.c_str());
      }
      // Write the sequence and quality lines
      gzprintf(gz_out, "%s\n", seq->seq.s);
      gzprintf(gz_out, "+%s\n", seq->comment.l ? seq->comment.s : "");
      gzprintf(gz_out, "%s\n", seq->qual.s);
      processed_reads++;
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << processed_reads << " reads..." << std::flush;
      }
    }
    
    std::cout << "\nFASTQ processing complete. Processed " << processed_reads << " reads...\n";
    // Clean up and close files
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(gz_filtered_out);
    gzclose(gz_unfiltered_out);
    std::cout << "FASTQ files updated.\n";
  }
  
  // Step 7: Update corrected_counts based on read_counts
  for (const auto& [barcode, count] : read_counts) {
    corrected_counts[barcode] = count;
  }
  
  // Step 8: Write barcode counts to CSV
  std::ofstream csv_file(counts_output_csv);
  csv_file << "barcode,original_counts,corrected_counts,total_counts\n";
  for (const auto& [barcode, orig_count] : original_counts) {
    int corrected_count = corrected_counts[barcode];  // Now reflects actual read counts
    csv_file << int64_to_string(barcode, sequence_length) << ","
             << orig_count << ","
             << corrected_count << ","
             << (orig_count + corrected_count) << "\n";
  }
  csv_file.close();
  
  // Step 9: Write detailed output to CSV if verbose_output is enabled
  if (verbose_output) {
    std::ofstream detailed_csv(detailed_output_csv);
    detailed_csv << "barcode,type,original,corrected\n";
    
    for (const auto& result : shifted_results) {
      detailed_csv << int64_to_string(result.first, sequence_length) << ",shifted,"
                   << int64_to_string(result.first, sequence_length) << ","
                   << int64_to_string(result.second, sequence_length) << "\n";
    }
    
    for (const auto& result : mutated_results) {
      detailed_csv << int64_to_string(result.first, sequence_length) << ",mutated,"
                   << int64_to_string(result.first, sequence_length) << ","
                   << int64_to_string(result.second, sequence_length) << "\n";
    }
    for (const auto& result : resolved_collisions) {
      if (result.second.size() == 1) {
        detailed_csv << int64_to_string(result.first, sequence_length) << ",resolved_single,"
                     << int64_to_string(result.first, sequence_length) << ","
                     << int64_to_string(result.second[0], sequence_length) << "\n";
      } else {
        // Unresolved barcodes: concatenate them with their distances
        std::string unresolved_str;
        for (size_t i = 0; i < result.second.size(); ++i) {
          std::string unresolved_barcode = int64_to_string(result.second[i], sequence_length);
          unresolved_str += unresolved_barcode;
          if (i != result.second.size() - 1) {
            unresolved_str += "|"; // Use '|' as the delimiter
          }
        }
        detailed_csv << int64_to_string(result.first, sequence_length) << ",unresolved,"
                     << int64_to_string(result.first, sequence_length) << ","
                     << unresolved_str << "\n";
      }
    }
    detailed_csv.close();
  }
  
  // Output summary statistics for verification
  std::cout << "Total reads processed: " << total_reads << "\n";
  std::cout << "Corrected reads: " << corrected_reads << "\n";
  std::cout << "Uncorrected reads: " << uncorrected_reads << "\n";
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  std::cout << "Total execution time: " << duration.count() << " seconds\n";
  return Rcpp::List::create(Rcpp::Named("message") = "Processing complete");
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(data.table)]]
// [[Rcpp::export]]
Rcpp::DataFrame generate_and_filter_mutations_v10(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    bool verbose = false) {
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> correct_barcodes = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int> correct_counts = Rcpp::as<std::vector<int>>(true_counts);
  
  std::vector<int64_t> incorrect_barcodes = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::unordered_set<int64_t> incorrect_set(incorrect_barcodes.begin(), incorrect_barcodes.end());
  
  if (correct_barcodes.empty()) {
    Rcpp::stop("No correct barcodes provided.");
  }
  
  if (incorrect_barcodes.empty()) {
    Rcpp::stop("No incorrect barcodes provided.");
  }
  
  // Step 2: Initialize incorrect barcode map (resolved, putative_barcodes, dl_distances)
  std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>> incorrect_map;
  for (const auto& incorrect_barcode : incorrect_barcodes) {
    incorrect_map[incorrect_barcode] = std::make_tuple(false, std::vector<std::pair<int64_t, int>>(), std::vector<double>());
  }
  
  size_t total_sequences = correct_barcodes.size();
  size_t processed_sequences = 0;
  size_t log_frequency = 10000;  // Log progress every 10,000 sequences
  
  // Step 3: Process each correct barcode and generate mutations and shifts
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    int64_t correct_barcode = correct_barcodes[i];
    
    // Generate mutated barcodes
    std::vector<int64_t> mutated_barcodes = generate_recursive_quaternary_mutations_cpp_v2(correct_barcode, mutation_rounds, sequence_length);
    for (const auto& mutated : mutated_barcodes) {
      if (incorrect_set.find(mutated) != incorrect_set.end()) {
#pragma omp critical
        std::get<1>(incorrect_map[mutated]).emplace_back(correct_barcode, 0);  // 0 for mutated
      }
    }
    
    // Generate shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(correct_barcode, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (incorrect_set.find(shifted) != incorrect_set.end()) {
#pragma omp critical
        std::get<1>(incorrect_map[shifted]).emplace_back(correct_barcode, 1);  // 1 for shifted
      }
    }
    
#pragma omp critical
{
  processed_sequences++;
  if (processed_sequences % log_frequency == 0 || processed_sequences == total_sequences) {
    size_t remaining_sequences = total_sequences - processed_sequences;
    std::cout << "Total barcodes to process: " << total_sequences
              << ", Processed: " << processed_sequences
              << ", Remaining: " << remaining_sequences << std::endl;
  }
}
  }
  
  std::cout << "Mutations and shifts generated.\n";
  std::cout << "Populating DL distances...\n";
  
  // Step 4: Populate DL distances (single-threaded, to avoid OpenMP issues)
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  int method = 2; // Damerau-Levenshtein distance method
  
  for (auto& [incorrect_barcode, data] : incorrect_map) {
    auto& [resolved, putative_barcodes, dl_distances] = data;
    
    if (putative_barcodes.size() > 1) {
      // Prepare incorrect barcode string
      std::string incorrect_str = int64_to_string(incorrect_barcode, sequence_length);
      
      // Prepare putative correct barcodes as strings
      std::vector<std::string> putative_correct_strs;
      for (const auto& [correct_barcode, _] : putative_barcodes) {
        putative_correct_strs.push_back(int64_to_string(correct_barcode, sequence_length));
      }
      
      // Compute DL distances using SEXP stringdist implementation (single-threaded)
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_str);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_strs);
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(0.0), Rcpp::wrap(0.0), Rcpp::wrap(1), Rcpp::wrap(0), Rcpp::wrap(1)  // Single thread
      );
      
      dl_distances = Rcpp::as<std::vector<double>>(dl_dist_results);
    }
  }
  
  std::cout << "DL distances populated.\n";
  
  // Step 5: Resolve collisions (parallelized)
  size_t total_incorrect = incorrect_map.size();
  log_frequency = 50000;  // Log progress every 50,000 sequences
  omp_set_num_threads(nthread);
  // Step 1: Create a vector of iterators to avoid using std::next
  std::vector<std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>>::iterator> iterators;
  iterators.reserve(incorrect_map.size());
  for (auto it = incorrect_map.begin(); it != incorrect_map.end(); ++it) {
    iterators.push_back(it);
  }
  // Atomic variable for progress tracking
  std::atomic<size_t> processed_incorrect(0);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < iterators.size(); ++idx) {
    auto iter = iterators[idx];  // Access precomputed iterator
    auto& [incorrect_barcode, data] = *iter;
    auto& [resolved, putative_barcodes, dl_distances] = data;
    // Use atomic increment to avoid critical section for the counter
    processed_incorrect++;
    // 1. Handle the simplest case: empty putative_barcodes
    if (putative_barcodes.empty()) {
      continue;  // Skip this incorrect barcode, no matches to resolve
    }
    // 2. Handle the case where there is exactly one putative barcode
    if (putative_barcodes.size() == 1) {
      resolved = true;  // Only one option, so it's resolved
      continue;  // Move to the next barcode
    }
    // 3. Handle the case where there are multiple putative barcodes (tie-breaking needed)
    if (putative_barcodes.size() > 1) {
      // Find the minimum DL distance
      double min_dl_dist = *std::min_element(dl_distances.begin(), dl_distances.end());
      std::vector<int> min_indices;
      for (int i = 0; i < dl_distances.size(); ++i) {
        if (dl_distances[i] == min_dl_dist) {
          min_indices.push_back(i);
        }
      }
      // 3a. If there's exactly one minimum index, resolve to that one barcode
      if (min_indices.size() == 1) {
        putative_barcodes = {putative_barcodes[min_indices[0]]};
        resolved = true;
      }
      // 3b. Otherwise, reduce to the set of barcodes with the same minimum distance
      else {
        std::vector<std::pair<int64_t, int>> reduced_barcodes;
        for (int index : min_indices) {
          reduced_barcodes.push_back(putative_barcodes[index]);
        }
        putative_barcodes = reduced_barcodes;
        resolved = false;  // Still unresolved due to ties
      }
    }
    // Print progress every 10,000 iterations (wrap print in a critical section)
    if (processed_incorrect % log_frequency == 0 || processed_incorrect == total_incorrect) {
#pragma omp critical
{
  size_t remaining_incorrect = total_incorrect - processed_incorrect.load();
  std::cout << "Total incorrect barcodes to process: " << total_incorrect
            << ", Processed: " << processed_incorrect.load()
            << ", Remaining: " << remaining_incorrect << std::endl;
}
    }
  }
  
  std::cout << "\nCollision resolution completed.\n";
  
  // Step 6: Prepare output using std::vector and parallel processing
  // Collect incorrect barcodes in a vector for direct access
  std::vector<int64_t> incorrect_keys;
  incorrect_keys.reserve(incorrect_map.size());
  for (const auto& pair : incorrect_map) {
    incorrect_keys.push_back(pair.first);
  }
  
  // Prepare the output vectors
  std::vector<std::string> incorrect_barcode_vec(incorrect_map.size());
  std::vector<std::string> corrected_barcode_vec(incorrect_map.size());
  std::vector<bool> resolved_status_vec(incorrect_map.size());
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    // Directly access the data using the pre-collected keys
    int64_t incorrect_barcode = incorrect_keys[idx];
    const auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    // Collect incorrect barcode
    incorrect_barcode_vec[idx] = int64_to_string(incorrect_barcode, sequence_length);
    
    // Collect corrected barcode(s)
    if (resolved) {
      corrected_barcode_vec[idx] = int64_to_string(putative_barcodes[0].first, sequence_length);
    } else {
      // If unresolved, join all putative barcodes with '|'
      std::string unresolved_barcodes;
      for (size_t j = 0; j < putative_barcodes.size(); ++j) {
        unresolved_barcodes += int64_to_string(putative_barcodes[j].first, sequence_length);
        if (j != putative_barcodes.size() - 1) {
          unresolved_barcodes += "|";
        }
      }
      corrected_barcode_vec[idx] = unresolved_barcodes;
    }
    
    // Collect resolution status
    resolved_status_vec[idx] = resolved;
  }
  
  // Step 7: Convert std::vector to Rcpp::StringVector and Rcpp::LogicalVector
  Rcpp::StringVector incorrect_barcode_out = Rcpp::wrap(incorrect_barcode_vec);
  Rcpp::StringVector corrected_barcode_out = Rcpp::wrap(corrected_barcode_vec);
  Rcpp::LogicalVector resolved_status_out = Rcpp::wrap(resolved_status_vec);
  
  // Step 8: Return DataTable as Rcpp::DataFrame
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  std::cout << "Total execution time: " << duration.count() << " seconds\n";
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("incorrect_barcode") = incorrect_barcode_out,
    Rcpp::Named("corrected_barcode") = corrected_barcode_out,
    Rcpp::Named("resolved_status") = resolved_status_out
  );
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(data.table)]]
// [[Rcpp::export]]
Rcpp::DataFrame generate_and_filter_mutations_v11(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    bool verbose = false,
    bool process_fastq = false,          // Flag for FASTQ processing
    std::string input_fastq = "",        // FASTQ input path
    std::string filtered_fastq = "filtered_fastq.gz",   // Output path for filtered FASTQ
    std::string unfiltered_fastq = "unfiltered_fastq.gz", // Output path for unfiltered FASTQ
    std::string barcode_header = "barcode:",  // FASTQ header identifier
    std::string detailed_output_csv = "",     // Output for detailed result CSV
    std::string counts_output_csv = "") {     // Output for counts CSV
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> correct_barcodes = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int> correct_counts = Rcpp::as<std::vector<int>>(true_counts);
  std::vector<int64_t> incorrect_barcodes = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::unordered_set<int64_t> incorrect_set(incorrect_barcodes.begin(), incorrect_barcodes.end());
  
  if (correct_barcodes.empty()) {
    Rcpp::stop("No correct barcodes provided.");
  }
  
  if (incorrect_barcodes.empty()) {
    Rcpp::stop("No incorrect barcodes provided.");
  }
  
  // Step 2: Initialize incorrect barcode map (resolved, putative_barcodes, dl_distances)
  std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>> incorrect_map;
  for (const auto& incorrect_barcode : incorrect_barcodes) {
    incorrect_map[incorrect_barcode] = std::make_tuple(false, std::vector<std::pair<int64_t, int>>(), std::vector<double>());
  }
  
  size_t total_sequences = correct_barcodes.size();
  size_t processed_sequences = 0;
  size_t log_frequency = 10000;  // Log progress every 10,000 sequences
  
  // Step 3: Process each correct barcode and generate mutations and shifts
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    int64_t correct_barcode = correct_barcodes[i];
    
    // Generate mutated barcodes
    std::vector<int64_t> mutated_barcodes = generate_recursive_quaternary_mutations_cpp_v2(correct_barcode, mutation_rounds, sequence_length);
    for (const auto& mutated : mutated_barcodes) {
      if (incorrect_set.find(mutated) != incorrect_set.end()) {
#pragma omp critical
        std::get<1>(incorrect_map[mutated]).emplace_back(correct_barcode, 0);  // 0 for mutated
      }
    }
    
    // Generate shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(correct_barcode, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (incorrect_set.find(shifted) != incorrect_set.end()) {
#pragma omp critical
        std::get<1>(incorrect_map[shifted]).emplace_back(correct_barcode, 1);  // 1 for shifted
      }
    }
    
#pragma omp critical
{
  processed_sequences++;
  if (processed_sequences % log_frequency == 0 || processed_sequences == total_sequences) {
    size_t remaining_sequences = total_sequences - processed_sequences;
    std::cout << "Total barcodes to process: " << total_sequences
              << ", Processed: " << processed_sequences
              << ", Remaining: " << remaining_sequences << std::endl;
  }
}
  }
  
  std::cout << "Mutations and shifts generated.\n";
  std::cout << "Populating DL distances...\n";
  
  // Step 4: Populate DL distances (single-threaded, to avoid OpenMP issues)
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  int method = 2; // Damerau-Levenshtein distance method
  
  for (auto& [incorrect_barcode, data] : incorrect_map) {
    auto& [resolved, putative_barcodes, dl_distances] = data;
    
    if (putative_barcodes.size() > 1) {
      std::string incorrect_str = int64_to_string(incorrect_barcode, sequence_length);
      std::vector<std::string> putative_correct_strs;
      for (const auto& [correct_barcode, _] : putative_barcodes) {
        putative_correct_strs.push_back(int64_to_string(correct_barcode, sequence_length));
      }
      
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_str);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_strs);
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(0.0), Rcpp::wrap(0.0), Rcpp::wrap(1), Rcpp::wrap(0), Rcpp::wrap(1)
      );
      dl_distances = Rcpp::as<std::vector<double>>(dl_dist_results);
    }
  }
  
  std::cout << "DL distances populated.\n";
  
  // Step 5: Resolve collisions (parallelized)
  size_t total_incorrect = incorrect_map.size();
  log_frequency = 50000;  // Log progress every 50,000 sequences
  std::atomic<size_t> processed_incorrect(0);
  std::vector<int64_t> incorrect_keys;
  incorrect_keys.reserve(incorrect_map.size());
  for (const auto& pair : incorrect_map) {
    incorrect_keys.push_back(pair.first);
  }
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    processed_incorrect++;
    
    if (putative_barcodes.empty()) {
      continue;
    }
    
    if (putative_barcodes.size() == 1) {
      resolved = true;
      continue;
    }
    
    if (putative_barcodes.size() > 1) {
      double min_dl_dist = *std::min_element(dl_distances.begin(), dl_distances.end());
      std::vector<int> min_indices;
      for (int i = 0; i < dl_distances.size(); ++i) {
        if (dl_distances[i] == min_dl_dist) {
          min_indices.push_back(i);
        }
      }
      
      if (min_indices.size() == 1) {
        putative_barcodes = {putative_barcodes[min_indices[0]]};
        resolved = true;
      } else {
        std::vector<std::pair<int64_t, int>> reduced_barcodes;
        for (int index : min_indices) {
          reduced_barcodes.push_back(putative_barcodes[index]);
        }
        putative_barcodes = reduced_barcodes;
        resolved = false;
      }
    }
    
    if (processed_incorrect % log_frequency == 0 || processed_incorrect == total_incorrect) {
#pragma omp critical
{
  size_t remaining_incorrect = total_incorrect - processed_incorrect.load();
  std::cout << "Total incorrect barcodes to process: " << total_incorrect
            << ", Processed: " << processed_incorrect.load()
            << ", Remaining: " << remaining_incorrect << std::endl;
}
    }
  }
  
  std::cout << "\nCollision resolution completed.\n";
  
  // Step 6: FASTQ Processing
  if (process_fastq) {
    std::cout << "\nProcessing and splitting FASTQ file...\n";
    
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    
    gzFile gz_filtered_out = gzopen(filtered_fastq.c_str(), "wb");
    gzFile gz_unfiltered_out = gzopen(unfiltered_fastq.c_str(), "wb");
    
    if (gz_filtered_out == NULL || gz_unfiltered_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ files.");
    }
    
    kseq_t* seq = kseq_init(fp);
    size_t processed_reads = 0;
    
    while (kseq_read(seq) >= 0) {
      std::string name(seq->name.s);
      size_t header_pos = name.find(barcode_header);
      gzFile gz_out = gz_unfiltered_out;
      
      if (header_pos != std::string::npos) {
        size_t barcode_start = header_pos + barcode_header.length();
        size_t barcode_end = name.find_first_of(" \t|", barcode_start);
        if (barcode_end == std::string::npos) {
          barcode_end = name.length();
        }
        std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
        std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
        if (!barcode_bits.empty()) {
          auto it = incorrect_map.find(barcode_bits[0]);
          if (it != incorrect_map.end()) {
            gz_out = gz_filtered_out;
          }
        }
      }
      
      gzprintf(gz_out, "@%s\n%s\n+\n%s\n", name.c_str(), seq->seq.s, seq->qual.s);
      processed_reads++;
      
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << processed_reads << " reads..." << std::flush;
      }
    }
    
    std::cout << "\nFASTQ processing complete. Processed " << processed_reads << " reads.\n";
    gzclose(fp);
    gzclose(gz_filtered_out);
    gzclose(gz_unfiltered_out);
  }
  
  // Output processing and return DataFrame
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  std::cout << "Total execution time: " << duration.count() << " seconds\n";
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("incorrect_barcode") = Rcpp::wrap(incorrect_keys)
  );
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(data.table)]]
// [[Rcpp::export]]
Rcpp::List generate_and_filter_mutations_v13(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    bool verbose = false,
    bool process_fastq = false,          // Flag for FASTQ processing
    std::string input_fastq = "",        // FASTQ input path
    std::string filtered_fastq = "filtered_fastq.gz",   // Output path for filtered FASTQ
    std::string unfiltered_fastq = "unfiltered_fastq.gz", // Output path for unfiltered FASTQ
    std::string barcode_header = "barcode:",  // FASTQ header identifier
    std::string detailed_output_csv = "detailed_output.csv", // Output for detailed result CSV
    std::string counts_output_csv = "barcode_counts.csv") {  // Output for barcode counts CSV
  
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> correct_barcodes = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::unordered_map<int64_t, int> original_counts, corrected_counts;  // Track counts for original and corrected barcodes
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    original_counts[correct_barcodes[i]] = 0;
    corrected_counts[correct_barcodes[i]] = 0;
  }
  
  std::vector<int64_t> incorrect_barcodes = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::unordered_set<int64_t> incorrect_set(incorrect_barcodes.begin(), incorrect_barcodes.end());
  
  if (correct_barcodes.empty()) {
    Rcpp::stop("No correct barcodes provided.");
  }
  
  if (incorrect_barcodes.empty()) {
    Rcpp::stop("No incorrect barcodes provided.");
  }
  
  // Step 2: Initialize incorrect barcode map (resolved, [putative_barcodes, mutation_type], dl_distances)
  std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>> incorrect_map;
  for (const auto& incorrect_barcode : incorrect_barcodes) {
    incorrect_map[incorrect_barcode] = std::make_tuple(false, std::vector<std::pair<int64_t, int>>(), std::vector<double>());
  }
  
  size_t total_sequences = correct_barcodes.size();
  size_t processed_sequences = 0;
  size_t log_frequency = 10000;  // Log progress every 10,000 sequences
  
  // Step 3: Process each correct barcode and generate mutations and shifts
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    int64_t correct_barcode = correct_barcodes[i];
    
    // Generate mutated barcodes
    std::vector<int64_t> mutated_barcodes = generate_recursive_quaternary_mutations_cpp_v2(
      correct_barcode,
      mutation_rounds,
      sequence_length);
    for (const auto& mutated : mutated_barcodes) {
      if (incorrect_set.find(mutated) != incorrect_set.end()) {
#pragma omp critical
        std::get<1>(incorrect_map[mutated]).emplace_back(correct_barcode, 0);  // 0 for mutated
      }
    }
    
    // Generate shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(correct_barcode, sequence_length, max_shift);
    for (const auto& shifted : shifted_barcodes) {
      if (incorrect_set.find(shifted) != incorrect_set.end()) {
#pragma omp critical
        std::get<1>(incorrect_map[shifted]).emplace_back(correct_barcode, 1);  // 1 for shifted
      }
    }
    
#pragma omp critical
{
  processed_sequences++;
  if (processed_sequences % log_frequency == 0 || processed_sequences == total_sequences) {
    size_t remaining_sequences = total_sequences - processed_sequences;
    std::cout << "Total barcodes to process: " << total_sequences
              << ", Processed: " << processed_sequences
              << ", Remaining: " << remaining_sequences << std::endl;
  }
}
  }
  
  std::cout << "Mutations and shifts generated.\n";
  std::cout << "Populating DL distances...\n";
  
  // Step 4: Populate DL distances (single-threaded, to avoid OpenMP issues)
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  int method = 2; // Damerau-Levenshtein distance method
  
  for (auto& [incorrect_barcode, data] : incorrect_map) {
    auto& [resolved, putative_barcodes, dl_distances] = data;
    if (putative_barcodes.size() > 1) {
      std::string incorrect_str = int64_to_string(incorrect_barcode, sequence_length);
      std::vector<std::string> putative_correct_strs;
      for (const auto& [correct_barcode, _] : putative_barcodes) {
        putative_correct_strs.push_back(int64_to_string(correct_barcode, sequence_length));
      }
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_str);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_strs);
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(0.0), Rcpp::wrap(0.0), Rcpp::wrap(1), Rcpp::wrap(0), Rcpp::wrap(1)  // Single thread
      );
      dl_distances = Rcpp::as<std::vector<double>>(dl_dist_results);
    }
  }
  
  std::cout << "DL distances populated.\n";
  
  // Step 5: Resolve collisions (parallelized)
  size_t total_incorrect = incorrect_map.size();
  log_frequency = 50000;  // Log progress every 50,000 sequences
  std::atomic<size_t> processed_incorrect(0);
  std::vector<int64_t> incorrect_keys;
  incorrect_keys.reserve(incorrect_map.size());
  for (const auto& pair : incorrect_map) {
    incorrect_keys.push_back(pair.first);
  }
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    processed_incorrect++;
    
    if (putative_barcodes.empty()) {
      resolved = false;
      continue;
    }
    
    if (putative_barcodes.size() == 1) {
      resolved = true;
      continue;
    }
    
    if (putative_barcodes.size() > 1) {
      // Create a map of unique barcodes and their corresponding minimum DL distances
      std::unordered_map<int64_t, double> barcode_to_min_dl_dist;
      
      // Loop through putative_barcodes and dl_distances to build the unique barcode set
      for (size_t i = 0; i < putative_barcodes.size(); ++i) {
        int64_t barcode = putative_barcodes[i].first;
        double dist = dl_distances[i];
        
        // Keep the minimum DL distance for each unique barcode
        if (barcode_to_min_dl_dist.find(barcode) == barcode_to_min_dl_dist.end()) {
          barcode_to_min_dl_dist[barcode] = dist;
        } else {
          barcode_to_min_dl_dist[barcode] = std::min(barcode_to_min_dl_dist[barcode], dist);
        }
      }
      
      // Now find the minimum DL distance from the unique barcode set
      auto min_it = std::min_element(
        barcode_to_min_dl_dist.begin(),
        barcode_to_min_dl_dist.end(),
        [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
          return a.second < b.second;
        }
      );
      
      double min_dl_dist = min_it->second;
      std::vector<int64_t> min_barcodes;
      
      // Collect all barcodes that have this minimum DL distance
      for (const auto& [barcode, dist] : barcode_to_min_dl_dist) {
        if (dist == min_dl_dist) {
          min_barcodes.push_back(barcode);
        }
      }
      
      // If exactly one barcode has the minimum DL distance, resolve to that barcode
      if (min_barcodes.size() == 1) {
        // Find and keep the pair in putative_barcodes that matches this barcode
        for (const auto& pair : putative_barcodes) {
          if (pair.first == min_barcodes[0]) {
            putative_barcodes = {pair};  // Resolve to this barcode
            break;
          }
        }
        resolved = true;  // Only one match, so it's resolved
      }
      // If there are multiple barcodes with the same minimum distance, check for uniqueness
      else {
        std::unordered_set<int64_t> unique_barcodes(min_barcodes.begin(), min_barcodes.end());
        
        // If all barcodes are the same (mutated/shifted), resolve to one of them
        if (unique_barcodes.size() == 1) {
          resolved = true;
          // Find and keep one of the matching pairs in putative_barcodes
          for (const auto& pair : putative_barcodes) {
            if (pair.first == *unique_barcodes.begin()) {
              putative_barcodes = {pair};  // Keep one of the pairs
              break;
            }
          }
        }
        // If not, keep the reduced set of barcodes and mark unresolved
        else {
          std::vector<std::pair<int64_t, int>> reduced_barcodes;
          for (const auto& pair : putative_barcodes) {
            if (unique_barcodes.count(pair.first)) {
              reduced_barcodes.push_back(pair);
            }
          }
          putative_barcodes = reduced_barcodes;
          resolved = false;  // Still unresolved due to multiple candidates
        }
      }
    }
    
    if (processed_incorrect % log_frequency == 0 || processed_incorrect == total_incorrect) {
#pragma omp critical
{
  size_t remaining_incorrect = total_incorrect - processed_incorrect.load();
  std::cout << "Total incorrect barcodes to process: " << total_incorrect
            << ", Processed: " << processed_incorrect.load()
            << ", Remaining: " << remaining_incorrect << std::endl;
}
    }
  }
  
  std::cout << "\nCollision resolution completed.\n";
  
  // Step 6: Prepare output vectors using std::vector and parallel processing
  std::vector<std::string> incorrect_barcode_vec(incorrect_map.size());
  std::vector<std::string> corrected_barcode_vec(incorrect_map.size());
  std::vector<bool> resolved_status_vec(incorrect_map.size());
  std::vector<std::string> mutation_type_vec(incorrect_map.size());
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    incorrect_barcode_vec[idx] = int64_to_string(incorrect_barcode, sequence_length);
    // If resolved, assign the correct barcode and mutation type
    if (resolved) {
      corrected_barcode_vec[idx] = int64_to_string(putative_barcodes[0].first, sequence_length);
      mutation_type_vec[idx] = (putative_barcodes[0].second == 0) ? "mutated" : "shifted";
      resolved_status_vec[idx] = resolved;  // Update resolved status
    }
    // Handle unresolved cases
    else {
      if (putative_barcodes.empty()) {
        mutation_type_vec[idx] = "unresolved_empty";
      } else {
        mutation_type_vec[idx] = "unresolved_multiple";
        // Fix the concatenation of multiple unresolved barcodes
        std::string unresolved_barcodes;
        for (size_t j = 0; j < putative_barcodes.size(); ++j) {
          unresolved_barcodes += int64_to_string(putative_barcodes[j].first, sequence_length);
          if (putative_barcodes.size() > 1 && j != putative_barcodes.size() - 1) {
            unresolved_barcodes += "|";  // Proper delimiter between barcodes
          }
        }
        corrected_barcode_vec[idx] = unresolved_barcodes;  // Store the concatenated unresolved barcodes
      }
      // No mutation type should be assigned for unresolved barcodes
      //mutation_type_vec[idx] = "unresolved";  // Reset mutation type for unresolved cases
    }
    
  }
  // Step 7: FASTQ Processing (Count reads in chunks)
  size_t total_reads = 0, corrected_reads = 0, uncorrected_reads = 0;
  size_t chunk_size = 100000;  // Default chunk size
  
  if (process_fastq) {
    std::cout << "\nProcessing and splitting FASTQ file in chunks...\n";
    
    gzFile fp = gzopen(input_fastq.c_str(), "r");
    if (fp == NULL) {
      Rcpp::stop("Failed to open input FASTQ file.");
    }
    
    gzFile gz_filtered_out = gzopen(filtered_fastq.c_str(), "wb");
    gzFile gz_unfiltered_out = gzopen(unfiltered_fastq.c_str(), "wb");
    
    if (gz_filtered_out == NULL || gz_unfiltered_out == NULL) {
      gzclose(fp);
      Rcpp::stop("Failed to open output FASTQ files.");
    }
    
    kseq_t* seq = kseq_init(fp);
    
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> fastq_chunk;
    fastq_chunk.reserve(chunk_size);
    
    size_t processed_reads = 0;
    
    while (kseq_read(seq) >= 0) {
      std::string name(seq->name.s);
      std::string seq_data(seq->seq.s);
      std::string qual(seq->qual.s);
      std::string comment(seq->comment.l ? seq->comment.s : "+");  // Handle optional comment line
      
      fastq_chunk.emplace_back(name, seq_data, qual, comment);
      total_reads++;
      processed_reads++;
      
      if (fastq_chunk.size() == chunk_size) {
        // Process the current chunk in parallel
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < fastq_chunk.size(); ++i) {
          auto& [name, seq_data, qual, comment] = fastq_chunk[i];
          size_t header_pos = name.find(barcode_header);
          gzFile gz_out = gz_unfiltered_out;
          
          if (header_pos != std::string::npos) {
            size_t barcode_start = header_pos + barcode_header.length();
            size_t barcode_end = name.find_first_of(" \t|", barcode_start);
            if (barcode_end == std::string::npos) {
              barcode_end = name.length();
            }
            std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
            std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
            
            if (!barcode_bits.empty()) {
              auto it = incorrect_map.find(barcode_bits[0]);
              if (it != incorrect_map.end()) {
#pragma omp critical
{
  gz_out = gz_filtered_out;
  corrected_reads++;
  corrected_counts[barcode_bits[0]]++;
}
              } else {
#pragma omp critical
{
  uncorrected_reads++;
  original_counts[barcode_bits[0]]++;
}
              }
            }
          }
          
          // Write to output file
#pragma omp critical
{
  gzprintf(gz_out, "@%s\n%s\n+\n%s\n", name.c_str(), seq_data.c_str(), qual.c_str());
}
        }
        fastq_chunk.clear();  // Clear the chunk buffer
      }
      if (processed_reads % 1000000 == 0) {
        std::cout << "\rProcessed " << total_reads << " reads..." << std::flush;
      }
    }
    
    // Process the remaining chunk if it's not empty
    if (!fastq_chunk.empty()) {
#pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < fastq_chunk.size(); ++i) {
        auto& [name, seq_data, qual, comment] = fastq_chunk[i];
        size_t header_pos = name.find(barcode_header);
        gzFile gz_out = gz_unfiltered_out;
        
        if (header_pos != std::string::npos) {
          size_t barcode_start = header_pos + barcode_header.length();
          size_t barcode_end = name.find_first_of(" \t|", barcode_start);
          if (barcode_end == std::string::npos) {
            barcode_end = name.length();
          }
          std::string barcode_str = name.substr(barcode_start, barcode_end - barcode_start);
          std::vector<int64_t> barcode_bits = sequence_to_bits_cpp(barcode_str);
          
          if (!barcode_bits.empty()) {
            auto it = incorrect_map.find(barcode_bits[0]);
            if (it != incorrect_map.end()) {
#pragma omp critical
{
  gz_out = gz_filtered_out;
  corrected_reads++;
  corrected_counts[barcode_bits[0]]++;
}
            } else {
#pragma omp critical
{
  uncorrected_reads++;
  original_counts[barcode_bits[0]]++;
}
            }
          }
        }
        
#pragma omp critical
{
  gzprintf(gz_out, "@%s\n%s\n+\n%s\n", name.c_str(), seq_data.c_str(), qual.c_str());
}
      }
    }
    
    std::cout << "\nFASTQ processing complete. Total Reads: " << total_reads << ", Corrected Reads: " << corrected_reads << ", Uncorrected Reads: " << uncorrected_reads << "\n";
    gzclose(fp);
    gzclose(gz_filtered_out);
    gzclose(gz_unfiltered_out);
  }
  
  // Step 8: Prepare barcode counts output and DataFrames
  std::vector<std::string> whitelist_barcode_vec(correct_barcodes.size());
  std::vector<int> total_counts_vec(correct_barcodes.size());
  
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    whitelist_barcode_vec[i] = int64_to_string(correct_barcodes[i], sequence_length);
    total_counts_vec[i] = original_counts[correct_barcodes[i]] + corrected_counts[correct_barcodes[i]];
  }
  
  std::vector<std::string> original_barcode_vec, corrected_barcode_output_vec;
  std::vector<int> original_counts_vec, corrected_counts_output_vec, final_total_counts_vec;
  
  for (const auto& pair : original_counts) {
    original_barcode_vec.push_back(int64_to_string(pair.first, sequence_length));
    original_counts_vec.push_back(pair.second);                 
    corrected_counts_output_vec.push_back(corrected_counts[pair.first]); 
    final_total_counts_vec.push_back(pair.second + corrected_counts[pair.first]);
  }
  
  Rcpp::DataFrame barcode_counts = Rcpp::DataFrame::create(
    Rcpp::Named("whitelist_barcode") = Rcpp::wrap(original_barcode_vec),
    Rcpp::Named("original_counts") = Rcpp::wrap(original_counts_vec),
    Rcpp::Named("corrected_counts") = Rcpp::wrap(corrected_counts_output_vec),
    Rcpp::Named("total_counts") = Rcpp::wrap(final_total_counts_vec)
  );
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  std::cout << "Total execution time: " << duration.count() << " seconds\n";
  
  return Rcpp::List::create(
    Rcpp::Named("detailed_output") = Rcpp::DataFrame::create(
      Rcpp::Named("incorrect_barcode") = Rcpp::wrap(incorrect_barcode_vec),
      Rcpp::Named("corrected_barcode") = Rcpp::wrap(corrected_barcode_vec),
      Rcpp::Named("resolved_status") = Rcpp::wrap(resolved_status_vec),
      Rcpp::Named("mutation_type") = Rcpp::wrap(mutation_type_vec)
    ),
    Rcpp::Named("barcode_counts") = barcode_counts
  );
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(data.table)]]
// [[Rcpp::export]]
Rcpp::DataFrame barcode_correction_v3(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    bool verbose = false) {
  //auto start_time = std::chrono::high_resolution_clock::now();
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> correct_barcodes = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int64_t> incorrect_barcodes = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::unordered_map<int64_t, int> original_counts, corrected_counts;  // Track counts for original and corrected barcodes
  
  // Step 2: Deduplicate correct barcodes and consolidate counts
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    original_counts[correct_barcodes[i]] += Rcpp::as<std::vector<int>>(true_counts)[i];  // Consolidate counts
    corrected_counts[correct_barcodes[i]] = 0;  // Initialize corrected counts
  }
  
  // Deduplicate incorrect barcodes
  std::unordered_set<int64_t> incorrect_set(incorrect_barcodes.begin(), incorrect_barcodes.end());
  
  if (correct_barcodes.empty()) {
    Rcpp::stop("No correct barcodes provided.");
  }
  if (incorrect_barcodes.empty()) {
    Rcpp::stop("No incorrect barcodes provided.");
  }
  
  // Step 3: Initialize incorrect barcode map (resolved, [putative_barcodes, mutation_type], dl_distances)
  std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>> incorrect_map;
  for (const auto& incorrect_barcode : incorrect_barcodes) {
    incorrect_map[incorrect_barcode] = std::make_tuple(false, std::vector<std::pair<int64_t, int>>(), std::vector<double>());
  }
  
  size_t total_sequences = correct_barcodes.size();
  size_t processed_sequences = 0;
  size_t log_frequency = 10000;  // Log progress every 10,000 sequences
  
  // Step 4: Generate mutations and shifts, including checking for both mutations and shifts
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    int64_t correct_barcode = correct_barcodes[i];
    
    std::unordered_set<int64_t> mutated_set;
    std::unordered_set<int64_t> shifted_set;
    
    // Generate mutated barcodes
    std::vector<int64_t> mutated_barcodes = generate_recursive_quaternary_mutations_cpp_v2(
      correct_barcode,
      mutation_rounds,
      sequence_length);
    
    mutated_set.insert(mutated_barcodes.begin(), mutated_barcodes.end());
    
    // Generate shifted barcodes
    if(max_shift > 0){
      std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(correct_barcode, sequence_length, max_shift);
      shifted_set.insert(shifted_barcodes.begin(), shifted_barcodes.end());
    }
    // Combine both sets to track "mutated", "shifted", and "mutated_and_shifted"
    std::unordered_set<int64_t> unique_barcodes(mutated_set.begin(), mutated_set.end());
    unique_barcodes.insert(shifted_set.begin(), shifted_set.end());
    
    for (const auto& barcode : unique_barcodes) {
      if (incorrect_set.find(barcode) != incorrect_set.end()) {
        int mutation_type = 0;  // Default: mutated
        if (mutated_set.find(barcode) != mutated_set.end() && shifted_set.find(barcode) != shifted_set.end()) {
          mutation_type = 2;  // mutated and shifted
        } else if (shifted_set.find(barcode) != shifted_set.end()) {
          mutation_type = 1;  // shifted
        }
#pragma omp critical
        std::get<1>(incorrect_map[barcode]).emplace_back(correct_barcode, mutation_type);
      }
    }
    
#pragma omp critical
{
  processed_sequences++;
  if (processed_sequences % log_frequency == 0 || processed_sequences == total_sequences) {
    size_t remaining_sequences = total_sequences - processed_sequences;
    std::cout << "Total barcodes to process: " << total_sequences
              << ", Processed: " << processed_sequences
              << ", Remaining: " << remaining_sequences << std::endl;
  }
}
  }
  
  std::cout << "Mutations and shifts generated.\n";
  std::cout << "Populating DL distances...\n";
  
  // Step 5: Populate DL distances
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  int method = 2; // Damerau-Levenshtein distance method
  
  for (auto& [incorrect_barcode, data] : incorrect_map) {
    auto& [resolved, putative_barcodes, dl_distances] = data;
    if (putative_barcodes.size() > 1) {
      std::string incorrect_str = int64_to_string(incorrect_barcode, sequence_length);
      std::vector<std::string> putative_correct_strs;
      for (const auto& [correct_barcode, _] : putative_barcodes) {
        putative_correct_strs.push_back(int64_to_string(correct_barcode, sequence_length));
      }
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_str);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_strs);
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(0.0), Rcpp::wrap(0.0), Rcpp::wrap(1), Rcpp::wrap(0), Rcpp::wrap(1)  // Single thread
      );
      dl_distances = Rcpp::as<std::vector<double>>(dl_dist_results);
    }
  }
  
  std::cout << "DL distances populated.\n";
  
  // Step 6: Resolve collisions
  //size_t total_incorrect = incorrect_map.size();
  log_frequency = 50000;  // Log progress every 50,000 sequences
  std::atomic<size_t> processed_incorrect(0);
  std::vector<int64_t> incorrect_keys;
  incorrect_keys.reserve(incorrect_map.size());
  for (const auto& pair : incorrect_map) {
    incorrect_keys.push_back(pair.first);
  }
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    processed_incorrect++;
    
    if (putative_barcodes.empty()) {
      resolved = false;
      continue;
    }
    
    if (putative_barcodes.size() == 1) {
      resolved = true;
      continue;
    }
    
    // Find the minimum DL distance and handle uniqueness
    std::unordered_map<int64_t, double> barcode_to_min_dl_dist;
    for (size_t i = 0; i < putative_barcodes.size(); ++i) {
      int64_t barcode = putative_barcodes[i].first;
      double dist = dl_distances[i];
      barcode_to_min_dl_dist[barcode] = std::min(barcode_to_min_dl_dist[barcode], dist);
    }
    
    // Find the barcode with the smallest DL distance
    auto min_it = std::min_element(
      barcode_to_min_dl_dist.begin(),
      barcode_to_min_dl_dist.end(),
      [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
        return a.second < b.second;
      }
    );
    
    double min_dl_dist = min_it->second;
    std::vector<int64_t> min_barcodes;
    
    for (const auto& [barcode, dist] : barcode_to_min_dl_dist) {
      if (dist == min_dl_dist) {
        min_barcodes.push_back(barcode);
      }
    }
    
    if (min_barcodes.size() == 1) {
      for (const auto& pair : putative_barcodes) {
        if (pair.first == min_barcodes[0]) {
          putative_barcodes = {pair};
          break;
        }
      }
      resolved = true;
    } else {
      resolved = false;
    }
  }
  
  std::cout << "Collision resolution completed.\n";
  
  // Step 7: Prepare outputs
  std::vector<std::string> incorrect_barcode_vec(incorrect_map.size());
  std::vector<std::string> corrected_barcode_vec(incorrect_map.size());
  std::vector<bool> resolved_status_vec(incorrect_map.size());
  std::vector<std::string> mutation_type_vec(incorrect_map.size());
  
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    incorrect_barcode_vec[idx] = int64_to_string(incorrect_barcode, sequence_length);
    if (resolved) {
      corrected_barcode_vec[idx] = int64_to_string(putative_barcodes[0].first, sequence_length);
      mutation_type_vec[idx] = (putative_barcodes[0].second == 0) ? "mutated" :
        (putative_barcodes[0].second == 1) ? "shifted" : "mutated_and_shifted";
    } else {
      if (putative_barcodes.empty()) {
        mutation_type_vec[idx] = "unresolved_empty";
      } else {
        mutation_type_vec[idx] = "unresolved_multiple";
        std::string unresolved_barcodes;
        for (size_t j = 0; j < putative_barcodes.size(); ++j) {
          unresolved_barcodes += int64_to_string(putative_barcodes[j].first, sequence_length);
          if (putative_barcodes.size() > 1 && j != putative_barcodes.size() - 1) {
            unresolved_barcodes += "|";
          }
        }
        corrected_barcode_vec[idx] = unresolved_barcodes;
      }
    }
    resolved_status_vec[idx] = resolved;
  }
  // Step 8: Return as DataFrame
  return Rcpp::DataFrame::create(
    Rcpp::Named("incorrect_barcodes") = Rcpp::wrap(incorrect_barcode_vec),
    Rcpp::Named("corrected_barcodes") = Rcpp::wrap(corrected_barcode_vec),
    Rcpp::Named("resolved_status") = Rcpp::wrap(resolved_status_vec),
    Rcpp::Named("mutation_type") = Rcpp::wrap(mutation_type_vec)
  );
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(data.table)]]
// [[Rcpp::export]]
Rcpp::DataFrame barcode_correction_v4(
    SEXP true_barcodes,
    SEXP invalid_barcodes,
    SEXP true_counts,
    SEXP invalid_counts,
    int mutation_rounds = 3,
    int sequence_length = 16,
    int nthread = 1,
    int max_shift = 3,
    bool verbose = false) {
  
  // Step 1: Initialize input vectors and validate
  std::vector<int64_t> correct_barcodes = Rcpp::as<std::vector<int64_t>>(true_barcodes);
  std::vector<int64_t> incorrect_barcodes = Rcpp::as<std::vector<int64_t>>(invalid_barcodes);
  std::unordered_map<int64_t, int> original_counts, corrected_counts;  // Track counts for original and corrected barcodes
  
  // Step 2: Deduplicate correct barcodes and consolidate counts
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    original_counts[correct_barcodes[i]] += Rcpp::as<std::vector<int>>(true_counts)[i];  // Consolidate counts
    corrected_counts[correct_barcodes[i]] = 0;  // Initialize corrected counts
  }
  
  // Deduplicate incorrect barcodes
  std::unordered_set<int64_t> incorrect_set(incorrect_barcodes.begin(), incorrect_barcodes.end());
  
  if (correct_barcodes.empty()) {
    Rcpp::stop("No correct barcodes provided.");
  }
  if (incorrect_barcodes.empty()) {
    Rcpp::stop("No incorrect barcodes provided.");
  }
  
  // Step 3: Initialize incorrect barcode map (resolved, [putative_barcodes, mutation_type], dl_distances)
  std::unordered_map<int64_t, std::tuple<bool, std::vector<std::pair<int64_t, int>>, std::vector<double>>> incorrect_map;
  for (const auto& incorrect_barcode : incorrect_barcodes) {
    incorrect_map[incorrect_barcode] = std::make_tuple(false, std::vector<std::pair<int64_t, int>>(), std::vector<double>());
  }
  
  size_t total_sequences = correct_barcodes.size();
  size_t processed_sequences = 0;
  size_t log_frequency = 10000;  // Log progress every 10,000 sequences
  
  // Step 4: Generate mutations and shifts, including checking for both mutations and shifts
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < correct_barcodes.size(); ++i) {
    int64_t correct_barcode = correct_barcodes[i];
    
    std::unordered_set<int64_t> mutated_set;
    std::unordered_set<int64_t> shifted_set;
    
    // Generate mutated barcodes
    std::vector<int64_t> mutated_barcodes = generate_recursive_quaternary_mutations_cpp_v2(
      correct_barcode,
      mutation_rounds,
      sequence_length);
    
    mutated_set.insert(mutated_barcodes.begin(), mutated_barcodes.end());
    
    // Generate shifted barcodes
    std::vector<int64_t> shifted_barcodes = generate_shifted_barcodes(correct_barcode, sequence_length, max_shift);
    
    shifted_set.insert(shifted_barcodes.begin(), shifted_barcodes.end());
    
    // Combine both sets to track "mutated", "shifted", and "mutated_and_shifted"
    std::unordered_set<int64_t> unique_barcodes(mutated_set.begin(), mutated_set.end());
    unique_barcodes.insert(shifted_set.begin(), shifted_set.end());
    
    for (const auto& barcode : unique_barcodes) {
      if (incorrect_set.find(barcode) != incorrect_set.end()) {
        int mutation_type = 0;  // Default: mutated
        if (mutated_set.find(barcode) != mutated_set.end() && shifted_set.find(barcode) != shifted_set.end()) {
          mutation_type = 2;  // mutated and shifted
        } else if (shifted_set.find(barcode) != shifted_set.end()) {
          mutation_type = 1;  // shifted
        }
#pragma omp critical
        std::get<1>(incorrect_map[barcode]).emplace_back(correct_barcode, mutation_type);
      }
    }
    
#pragma omp critical
{
  processed_sequences++;
  if (processed_sequences % log_frequency == 0 || processed_sequences == total_sequences) {
    size_t remaining_sequences = total_sequences - processed_sequences;
    std::cout << "Total barcodes to process: " << total_sequences
              << ", Processed: " << processed_sequences
              << ", Remaining: " << remaining_sequences << std::endl;
  }
}
  }
  
  std::cout << "Mutations and shifts generated.\n";
  std::cout << "Populating DL distances...\n";
  
  // Step 5: Populate DL distances
  Rcpp::NumericVector weight = Rcpp::NumericVector::create(1.0, 1.0, 1.0, 1.0);
  int method = 2; // Damerau-Levenshtein distance method
  
  for (auto& [incorrect_barcode, data] : incorrect_map) {
    auto& [resolved, putative_barcodes, dl_distances] = data;
    if (putative_barcodes.size() > 1) {
      std::string incorrect_str = int64_to_string(incorrect_barcode, sequence_length);
      std::vector<std::string> putative_correct_strs;
      for (const auto& [correct_barcode, _] : putative_barcodes) {
        putative_correct_strs.push_back(int64_to_string(correct_barcode, sequence_length));
      }
      SEXP incorrect_sexp = Rcpp::wrap(incorrect_str);
      SEXP correct_sexp = Rcpp::wrap(putative_correct_strs);
      SEXP dl_dist_results = sd_stringdist(
        incorrect_sexp,
        correct_sexp,
        Rcpp::wrap(method),
        weight, Rcpp::wrap(0.0), Rcpp::wrap(0.0), Rcpp::wrap(1), Rcpp::wrap(0), Rcpp::wrap(1)  // Single thread
      );
      dl_distances = Rcpp::as<std::vector<double>>(dl_dist_results);
    }
  }
  
  std::cout << "DL distances populated.\n";
  
  // Step 6: Resolve collisions and update corrected counts
  log_frequency = 50000;  // Log progress every 50,000 sequences
  std::atomic<size_t> processed_incorrect(0);
  std::vector<int64_t> incorrect_keys;
  incorrect_keys.reserve(incorrect_map.size());
  for (const auto& pair : incorrect_map) {
    incorrect_keys.push_back(pair.first);
  }
  
  omp_set_num_threads(nthread);
#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    processed_incorrect++;
    
    if (putative_barcodes.empty()) {
      resolved = false;
      continue;
    }
    
    if (putative_barcodes.size() == 1) {
      resolved = true;
      // Increment corrected counts
      int64_t correct_barcode = putative_barcodes[0].first;
#pragma omp critical
{
  corrected_counts[correct_barcode] += original_counts[correct_barcode];
}
continue;
    }
    
    // Find the minimum DL distance and handle uniqueness
    std::unordered_map<int64_t, double> barcode_to_min_dl_dist;
    for (size_t i = 0; i < putative_barcodes.size(); ++i) {
      int64_t barcode = putative_barcodes[i].first;
      double dist = dl_distances[i];
      barcode_to_min_dl_dist[barcode] = std::min(barcode_to_min_dl_dist[barcode], dist);
    }
    
    // Find the barcode with the smallest DL distance
    auto min_it = std::min_element(
      barcode_to_min_dl_dist.begin(),
      barcode_to_min_dl_dist.end(),
      [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
        return a.second < b.second;
      }
    );
    
    double min_dl_dist = min_it->second;
    std::vector<int64_t> min_barcodes;
    
    for (const auto& [barcode, dist] : barcode_to_min_dl_dist) {
      if (dist == min_dl_dist) {
        min_barcodes.push_back(barcode);
      }
    }
    
    if (min_barcodes.size() == 1) {
      for (const auto& pair : putative_barcodes) {
        if (pair.first == min_barcodes[0]) {
          putative_barcodes = {pair};
          break;
        }
      }
      resolved = true;
      // Increment corrected counts
      int64_t correct_barcode = putative_barcodes[0].first;
#pragma omp critical
{
  corrected_counts[correct_barcode] += original_counts[correct_barcode];
}
    } else {
      resolved = false;
    }
  }
  
  std::cout << "Collision resolution completed.\n";
  
  // Step 7: Prepare outputs
  std::vector<std::string> incorrect_barcode_vec(incorrect_map.size());
  std::vector<std::string> corrected_barcode_vec(incorrect_map.size());
  std::vector<bool> resolved_status_vec(incorrect_map.size());
  std::vector<std::string> mutation_type_vec(incorrect_map.size());
  std::vector<int> corrected_count_vec(incorrect_map.size());
  
  for (size_t idx = 0; idx < incorrect_keys.size(); ++idx) {
    int64_t incorrect_barcode = incorrect_keys[idx];
    auto& [resolved, putative_barcodes, dl_distances] = incorrect_map[incorrect_barcode];
    
    incorrect_barcode_vec[idx] = int64_to_string(incorrect_barcode, sequence_length);
    if (resolved) {
      corrected_barcode_vec[idx] = int64_to_string(putative_barcodes[0].first, sequence_length);
      mutation_type_vec[idx] = (putative_barcodes[0].second == 0) ? "mutated" :
        (putative_barcodes[0].second == 1) ? "shifted" : "mutated_and_shifted";
      corrected_count_vec[idx] = corrected_counts[putative_barcodes[0].first];
    } else {
      if (putative_barcodes.empty()) {
        mutation_type_vec[idx] = "unresolved_empty";
        corrected_count_vec[idx] = 0;
      } else {
        mutation_type_vec[idx] = "unresolved_multiple";
        std::string unresolved_barcodes;
        for (size_t j = 0; j < putative_barcodes.size(); ++j) {
          unresolved_barcodes += int64_to_string(putative_barcodes[j].first, sequence_length);
          if (putative_barcodes.size() > 1 && j != putative_barcodes.size() - 1) {
            unresolved_barcodes += "|";
          }
        }
        corrected_barcode_vec[idx] = unresolved_barcodes;
        corrected_count_vec[idx] = 0;
      }
    }
    resolved_status_vec[idx] = resolved;
  }
  
  // Step 8: Return as DataFrame
  return Rcpp::DataFrame::create(
    Rcpp::Named("incorrect_barcodes") = Rcpp::wrap(incorrect_barcode_vec),
    Rcpp::Named("corrected_barcodes") = Rcpp::wrap(corrected_barcode_vec),
    Rcpp::Named("resolved_status") = Rcpp::wrap(resolved_status_vec),
    Rcpp::Named("mutation_type") = Rcpp::wrap(mutation_type_vec),
    Rcpp::Named("corrected_count") = Rcpp::wrap(corrected_count_vec)
  );
}