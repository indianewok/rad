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
  std::vector<int64_t> shifted_barcodes;
  int64_t mask = (1LL << (2 * sequence_length)) - 1;  // Mask for the original sequence length
  // Generate right shifts
  for (int shift = 1; shift <= max_shift; ++shift) {
    for (int padding = 0; padding < 4; ++padding) {
      int64_t shifted = ((original >> (2 * shift)) & (mask >> (2 * shift))) | (padding << (2 * (sequence_length - shift)));
      shifted_barcodes.push_back(shifted);
    }
  }
  // Generate left shifts
  for (int shift = 1; shift <= max_shift; ++shift) {
    for (int padding = 0; padding < 4; ++padding) {
      int64_t shifted = ((original << (2 * shift)) & mask) | padding;
      shifted_barcodes.push_back(shifted);
    }
  }
  return shifted_barcodes;
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
            std::string corrected_barcode = int64_to_string(it->second, sequence_length);
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
