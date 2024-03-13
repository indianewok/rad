#include "rad.h"
using namespace std;
using namespace Rcpp;

whitelist wl; 

#define MIN(a,b) (((a)<(b))?(a):(b))
typedef struct dictionary item;
static __inline item* push(unsigned int key, item* curr){
  item* head;
  head = (item*)malloc(sizeof(item));
  head->key = key;
  head->value = 0;
  head->next = curr;
  return head;
}
static __inline item* find(item* head, unsigned int key){
  item* iterator = head;
  while(iterator){
    if(iterator->key == key){
      return iterator;
    }
    iterator = iterator->next;
  }
  return NULL;
}
static __inline item* uniquePush(item* head, unsigned int key){
  item* iterator = head;
  while(iterator){
    if(iterator->key == key){
      return head;
    }
    iterator = iterator->next;
  }
  return push(key, head);
}
static void dict_free(item* head){
  item* iterator = head;
  while(iterator){
    item* temp = iterator;
    iterator = iterator->next;
    free(temp);
  }
  head = NULL;
}

static int dl_dist_main(int64_t src_input, int64_t tgt_input, unsigned int maxDistance){
  
  const unsigned int x = 16; // Hardcoded sequence length
  const unsigned int y = 16; // Hardcoded sequence length
  
  std::vector<int64_t> src_input_vec = {src_input}; // Wrap src_input in a vector
  std::vector<int64_t> tgt_input_vec = {tgt_input}; // Wrap tgt_input in a vector
  
  std::vector<unsigned int> src_vec = bits_to_uint_cpp(src_input_vec, x);
  std::vector<unsigned int> tgt_vec = bits_to_uint_cpp(tgt_input_vec, y);
  
  unsigned int* src = src_vec.data();  // Equivalent to &src_vec[0]
  unsigned int* tgt = tgt_vec.data();  // Equivalent to &tgt_vec[0]
  
  item *head = NULL;
  unsigned int swapCount, swapScore, targetCharCount, i, j;
  unsigned int *scores = (unsigned int*)malloc( (x + 2) * (y + 2) * sizeof(unsigned int) );
  unsigned int score_ceil = x + y;
  unsigned int curr_score;
  unsigned int diff = x > y ? (x - y) : (y - x);
  if(maxDistance != 0 && diff > maxDistance) {
    free(scores);
    return -1;
  }
  scores[0] = score_ceil;
  scores[1 * (y + 2) + 0] = score_ceil;
  scores[0 * (y + 2) + 1] = score_ceil;
  scores[1 * (y + 2) + 1] = 0;
  head = uniquePush(uniquePush(head, src[0]), tgt[0]);
  for(i = 1; i <= x; i++) {
    if(i < x)
      head = uniquePush(head, src[i]);
    scores[(i + 1) * (y + 2) + 1] = i;
    scores[(i + 1) * (y + 2) + 0] = score_ceil;
    swapCount = 0;
    for(j = 1; j <= y; j++) {
      if(i == 1) {
        if(j < y)
          head = uniquePush(head, tgt[j]);
        scores[1 * (y + 2) + (j + 1)] = j;
        scores[0 * (y + 2) + (j + 1)] = score_ceil;
      }
      curr_score = 0;
      targetCharCount = find(head, tgt[j - 1])->value;
      swapScore = scores[targetCharCount * (y + 2) + swapCount] + i - targetCharCount - 1 + j - swapCount;
      if(src[i - 1] != tgt[j - 1]){
        scores[(i + 1) * (y + 2) + (j + 1)] = MIN(swapScore, (MIN(scores[i * (y + 2) + j], MIN(scores[(i + 1) * (y + 2) + j], scores[i * (y + 2) + (j + 1)])) + 1));
      } else {
        swapCount = j;
        scores[(i + 1) * (y + 2) + (j + 1)] = MIN(scores[i * (y + 2) + j], swapScore);
      }
      curr_score = MIN(curr_score, scores[(i + 1) * (y + 2) + (j + 1)]);
    }
    if(maxDistance != 0 && curr_score > maxDistance) {
      dict_free(head);
      free(scores);
      return -1;
    }
    find(head, src[i - 1])->value = i;
  }
  {
    unsigned int score = scores[(x + 1) * (y + 2) + (y + 1)];
    dict_free(head);
    free(scores);
    return (maxDistance != 0 && maxDistance < score) ? (-1) : score;
  }
}

std::vector<int> dl_dist_cpp(const std::vector<int64_t>& src_input_vec, 
  const std::vector<int64_t>& tgt_input_vec, 
  unsigned int maxDistance = 10) {
  if (src_input_vec.size() != 1) {
    throw std::invalid_argument("src_input_vec must contain exactly one element.");
  }
  std::vector<int> distances;
  distances.reserve(tgt_input_vec.size()); // Reserve space to avoid reallocations
  for (const auto& tgt_input : tgt_input_vec) {
    // Call dl_main for each element in tgt_input_vec with the single src_input
    int distance = dl_dist_main(src_input_vec[0], tgt_input, maxDistance);
    distances.push_back(distance);
  }
  return distances;
}

// [[Rcpp::export]]
Rcpp::IntegerVector dl_dist_rcpp(Rcpp::NumericVector src_input, Rcpp::NumericVector tgt_input, unsigned int maxDistance) {
  // Convert Rcpp::NumericVector to std::vector<int64_t>
  std::vector<int64_t> src_input_vec = Rcpp::as<std::vector<int64_t>>(src_input);
  std::vector<int64_t> tgt_input_vec = Rcpp::as<std::vector<int64_t>>(tgt_input);
  // Call the C++ version of dl_dist
  std::vector<int> distances = dl_dist_cpp(src_input_vec, tgt_input_vec, maxDistance);
  // Convert std::vector<int> back to Rcpp::IntegerVector to return to R
  return Rcpp::wrap(distances);
}

void generate_whitelist(const std::vector<int64_t>& barcodes) {
  wl.clear();
  for (int64_t barcode : barcodes) {
    wl[barcode] = barcode_data{};
  }
}

// [[Rcpp::export]]
void populate_whitelist(NumericVector barcodes, 
  NumericVector poisson_data = NumericVector(), 
  NumericVector counts = NumericVector(), bool verbose = false) {
  std::vector<int64_t> cpp_barcodes = as<std::vector<int64_t>>(barcodes);
  if (wl.empty()) {
    generate_whitelist(cpp_barcodes);
    if(verbose){
      Rcpp::Rcout << "This whitelist is empty!\n"; 
    }
  }
  bool has_poisson_data = poisson_data.size() == barcodes.size();
  bool has_count_data = counts.size() == barcodes.size();
  for (size_t i = 0; i < cpp_barcodes.size(); ++i) {
    int64_t barcode = cpp_barcodes[i];
    auto it = wl.find(barcode);
    if (it != wl.end()) {
      if(verbose){
        Rcpp::Rcout << "Barcode has been found!\n";
        Rcpp::Rcout << barcode << "\n";
      }
      // Existing entry found, update it
      if (has_poisson_data && !Rcpp::NumericVector::is_na(poisson_data[i])) {
        it->second.poisson_score = static_cast<double>(poisson_data[i]); // Ensure correct type
      }
      if (has_count_data && !Rcpp::NumericVector::is_na(counts[i])) {
        it->second.count = static_cast<int>(counts[i]);
      }
    } else {
      // Entry not found, insert a new one
      if(verbose){
        Rcpp::Rcout << "Barcode not found!\n";
      }
      barcode_data data;
      if (has_poisson_data && !Rcpp::NumericVector::is_na(poisson_data[i])) {
        data.poisson_score = static_cast<double>(poisson_data[i]);
      }
      if (has_count_data && !Rcpp::NumericVector::is_na(counts[i])) {
        data.count = static_cast<int>(counts[i]);
      }
      wl[barcode] = data;
    }
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame wl_to_df() {
  Rcpp::NumericVector sequences; // Use NumericVector to store int64_t values
  std::vector<int> counts;
  std::vector<double> poisson_scores;
  for (const auto& item : wl) {
    // Directly push the int64_t barcode value, cast to double for R compatibility
    sequences.push_back(static_cast<double>(item.first));
    counts.push_back(item.second.count);
    poisson_scores.push_back(item.second.poisson_score >= 0.0 ? item.second.poisson_score : NA_REAL);
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("sequence") = sequences,
    Rcpp::Named("count") = counts,
    Rcpp::Named("poisson_score") = poisson_scores,
    Rcpp::_["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
void clear_whitelist() {
  Rcpp::Rcout << "Clearing whitelist, previous size: " << wl.size() << "\n";
  wl.clear();
  Rcpp::Rcout << "Whitelist cleared, new size: " << wl.size() << "\n";
}

// [[Rcpp::export]]
int whitelist_size() {
  return wl.size();
}

// [[Rcpp::export]]
NumericVector correct_barcodes_v1(NumericVector barcodes, bool verbose = false, int nthreads = 1, 
  int depth = 2, int breadth = 1) {
  int n = barcodes.size();
  NumericVector results(n, NA_REAL);
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    std::vector<int64_t> candidates; // Store all possible candidates
    // First, handle original barcode mutations
    std::vector<int64_t> mutations = generate_recursive_mutations_cpp(barcode, depth);
    for (auto& mutation : mutations) {
      if (wl.find(mutation) != wl.end()) {
        candidates.push_back(mutation);
      }
    }
    // Then, handle circularized sequence mutations if no candidates found from original mutations
    if (candidates.empty()) {
      std::vector<int64_t> circularizedSequences = kmer_circ_cpp(barcode, 16, 16, verbose);
      for (auto& circSeq : circularizedSequences) {
        std::vector<int64_t> circMutations = generate_recursive_mutations_cpp(circSeq, breadth);
        for (auto& mutation : circMutations) {
          if (wl.find(mutation) != wl.end()) {
            candidates.push_back(mutation);
          }
        }
      }
    }
    // Rank candidates based on Poisson score and DL distance
    int64_t best_candidate = 0;
    double highest_poisson_score = -1.0;
    int best_dl_distance = INT_MAX;
    for (auto& candidate : candidates) {
      auto& candidate_data = wl[candidate];
      int dl_distance = dl_dist_cpp({barcode}, {candidate}, 10)[0];
      if(verbose){
        Rcpp::Rcout << dl_distance << "\n";
      }
      if ((candidate_data.poisson_score > highest_poisson_score) ||
        (candidate_data.poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
        best_candidate = candidate;
        highest_poisson_score = candidate_data.poisson_score;
        best_dl_distance = dl_distance;
      }
    }
    if (best_candidate != 0) {
      wl[best_candidate].count++; // Increment count for the best match
      results[i] = static_cast<double>(best_candidate); // Store the best candidate as a double
      if (verbose) {
#pragma omp critical
        Rcpp::Rcout << "Best match found for barcode: " << barcode << " with Poisson score: " << 
          highest_poisson_score << ", DL distance: " << best_dl_distance << std::endl;
      }
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes_v2(Rcpp::NumericVector barcodes, bool verbose = false,
  int nthreads = 1, int maxDistance = 3) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL);
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    int64_t best_match = 0;
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    for (const auto& wl_entry : wl) {
      std::vector<int> dl_distances = dl_dist_cpp({barcode}, {wl_entry.first}, maxDistance);
      if (!dl_distances.empty() && dl_distances[0] >= 0 && dl_distances[0] <= maxDistance) {
        double poisson_score = wl_entry.second.poisson_score;
        if (poisson_score > highest_poisson_score || 
          (poisson_score == highest_poisson_score && dl_distances[0] < best_dl_distance)) {
          best_match = wl_entry.first;
          highest_poisson_score = poisson_score;
          best_dl_distance = dl_distances[0];
        }
      }
    }
    if (best_match != 0) {
      wl[best_match].count++;
      results[i] = static_cast<double>(best_match);
      if (verbose) {
#pragma omp critical
        Rcpp::Rcout << "Best match for barcode " << barcode << ": " << best_match 
                    << " with Poisson score: " << highest_poisson_score 
                    << " and DL distance: " << best_dl_distance << std::endl;
      }
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes_v3(Rcpp::NumericVector barcodes, bool verbose = false, 
  int nthreads = 1, int breadth = 2, int depth = 2, int maxDistance = 3) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL);
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    int64_t best_match = 0;
    for (const auto& wl_entry : wl) {
      int dl_distance = dl_dist_cpp({barcode}, {wl_entry.first}, maxDistance)[0];
      if (dl_distance >= 0 && dl_distance <= maxDistance) {
        double poisson_score = wl_entry.second.poisson_score;
        if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
          best_match = wl_entry.first;
          highest_poisson_score = poisson_score;
          best_dl_distance = dl_distance;
        }
      }
    }
    if (best_match == 0) {
      auto mutations = generate_recursive_mutations_cpp(barcode, depth);
      for (auto mutation : mutations) {
        if (wl.find(mutation) != wl.end()) {
          int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
          if (dl_distance <= maxDistance) {
            double poisson_score = wl[mutation].poisson_score;
            if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
              best_match = mutation;
              highest_poisson_score = poisson_score;
              best_dl_distance = dl_distance;
            }
          }
        }
      }
      // Circularized sequence mutations
      if (best_match == 0) {
        auto circularizedSequences = kmer_circ_cpp(barcode, 16, 16, verbose);
        for (auto circSeq : circularizedSequences) {
          auto circMutations = generate_recursive_mutations_cpp(circSeq, breadth);
          for (auto mutation : circMutations) {
            if (wl.find(mutation) != wl.end()) {
              int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
              if (dl_distance <= maxDistance) {
                double poisson_score = wl[mutation].poisson_score;
                if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
                  best_match = mutation;
                  highest_poisson_score = poisson_score;
                  best_dl_distance = dl_distance;
                }
              }
            }
          }
        }
      }
    }
    // Update results
    if (best_match != 0) {
      wl[best_match].count++; // Increment the count for the best match
      results[i] = static_cast<double>(best_match); // Convert the best match to double and assign
      if (verbose) {
#pragma omp critical
        Rcpp::Rcout << "Best match found for barcode: " << barcode << " is: " << best_match 
                    << " with Poisson score: " << highest_poisson_score << " and DL distance: " << best_dl_distance << std::endl;
      }
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes_v4(Rcpp::NumericVector barcodes, bool verbose = false, 
  int nthreads = 1, int depth = 2, int breadth = 2, bool high_speed = false, int maxDistance = 3) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL);
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    int64_t best_match = 0;
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    auto mutations = generate_recursive_mutations_cpp(barcode, depth);
    for (auto mutation : mutations) {
      if (wl.find(mutation) != wl.end()) {
        int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
        if (dl_distance <= maxDistance && dl_distance > -1) {
          double poisson_score = wl[mutation].poisson_score;
          if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
            best_match = mutation;
            highest_poisson_score = poisson_score;
            best_dl_distance = dl_distance;
          }
        }
      }
    }
    if (!high_speed && best_match == 0) {
      auto circularizedSequences = kmer_circ_cpp(barcode, 16, 16, verbose);
      for (auto circSeq : circularizedSequences) {
        auto circMutations = generate_recursive_mutations_cpp(circSeq, breadth);
        for (auto mutation : circMutations) {
          if (wl.find(mutation) != wl.end()) {
            int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
            if (dl_distance <= maxDistance) {
              double poisson_score = wl[mutation].poisson_score;
              if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
                best_match = mutation;
                highest_poisson_score = poisson_score;
                best_dl_distance = dl_distance;
              }
            }
          }
        }
      }
    }
    if (best_match != 0) {
      wl[best_match].count++;
      results[i] = static_cast<double>(best_match);
      if (verbose) {
#pragma omp critical
        Rcpp::Rcout << "Best match found for barcode: " << barcode << " is: " << best_match 
                    << " with Poisson score: " << highest_poisson_score << " and DL distance: " << best_dl_distance << std::endl;
      }
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes_v5(Rcpp::NumericVector barcodes, bool verbose = false, 
  int nthreads = 1, int depth = 2, int breadth = 2, bool high_speed = false, int maxDistance = 3) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL); // Initialize results with NA_REAL to handle barcodes without valid corrections
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads); // Set the desired number of threads
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    int64_t best_match = 0;
    // First, try to find a match through mutations or circularized sequences if high_speed is false
    if (!high_speed) {
      auto mutations = generate_recursive_mutations_cpp(barcode, depth);
      for (auto mutation : mutations) {
        if (wl.find(mutation) != wl.end()) {
          int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
          if (dl_distance <= maxDistance) {
            double poisson_score = wl[mutation].poisson_score;
            if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
              best_match = mutation;
              highest_poisson_score = poisson_score;
              best_dl_distance = dl_distance;
            }
          }
        }
      }
      if (best_match == 0) {
        auto circularizedSequences = kmer_circ_cpp(barcode, 16, 16, verbose);
        for (auto circSeq : circularizedSequences) {
          auto circMutations = generate_recursive_mutations_cpp(circSeq, breadth);
          for (auto mutation : circMutations) {
            if (wl.find(mutation) != wl.end()) {
              int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
              if (dl_distance <= maxDistance) {
                double poisson_score = wl[mutation].poisson_score;
                if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
                  best_match = mutation;
                  highest_poisson_score = poisson_score;
                  best_dl_distance = dl_distance;
                }
              }
            }
          }
        }
      }
    }
    // If no match is found yet and high_speed is FALSE, do DL distance comparison against each member of the whitelist
    if (best_match == 0 && !high_speed) {
      for (const auto& wl_entry : wl) {
        int dl_distance = dl_dist_cpp({barcode}, {wl_entry.first}, maxDistance)[0];
        if (dl_distance <= maxDistance) {
          double poisson_score = wl_entry.second.poisson_score;
          if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
            best_match = wl_entry.first;
            highest_poisson_score = poisson_score;
            best_dl_distance = dl_distance;
          }
        }
      }
    }
    // Update results
    if (best_match != 0) {
      wl[best_match].count++; // Increment the count for the best match
      results[i] = static_cast<double>(best_match); // Convert the best match to double and assign
      if (verbose) {
#pragma omp critical
        Rcpp::Rcout << "Best match found for barcode: " << barcode << " is: " << best_match 
                    << " with Poisson score: " << highest_poisson_score << " and DL distance: " << best_dl_distance << std::endl;
      }
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes_v6(Rcpp::NumericVector barcodes, bool verbose = false, 
  int nthreads = 1, int maxDistance = 3) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL); // Initialize results with NA_REAL to handle barcodes without valid corrections
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads); // Set the desired number of threads
  
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    int64_t best_match = 0;
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    // Check against each member of the whitelist for DL distance
    for (const auto& wl_entry : wl) {
      int dl_distance = dl_dist_cpp({barcode}, {wl_entry.first}, maxDistance)[0];
      if (dl_distance < 0) { // If DL distance is under 0, it's considered a mismatch or error
        continue; // Skip to the next iteration without updating the best match
      }
      if (dl_distance <= maxDistance) {
        double poisson_score = wl_entry.second.poisson_score;
        if (poisson_score > highest_poisson_score || (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
          best_match = wl_entry.first;
          highest_poisson_score = poisson_score;
          best_dl_distance = dl_distance;
        }
      }
    }
    // If a match is found, update the results
    if (best_match != 0) {
      wl[best_match].count++; // Increment the count for the best match
      results[i] = static_cast<double>(best_match); // Convert the best match to double and assign
    } else {
      results[i] = NA_REAL; // Assign NA if no suitable match is found
    }
    if (verbose && best_match != 0) {
#pragma omp critical
      Rcpp::Rcout << "Best match found for barcode: " << barcode << " is: " << best_match 
                  << " with Poisson score: " << highest_poisson_score << " and DL distance: " << best_dl_distance << std::endl;
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes_v7(Rcpp::NumericVector barcodes, bool verbose = false, 
  int nthreads = 1, int depth = 2, int breadth = 2, bool high_speed = false, int maxDistance = 3) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL); // Initialize results with NA_REAL
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads); // Set the number of threads
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    int64_t best_match = 0;
    // First, check against each member of the whitelist for DL distance
    for (const auto& wl_entry : wl) {
      int dl_distance = dl_dist_cpp({barcode}, {wl_entry.first}, maxDistance)[0];
      if (dl_distance >= 0 && dl_distance <= maxDistance) {
        double poisson_score = wl_entry.second.poisson_score;
        if (poisson_score > highest_poisson_score || 
          (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
          best_match = wl_entry.first;
          highest_poisson_score = poisson_score;
          best_dl_distance = dl_distance;
        }
      }
    }
    // If no match is found and high-speed mode is false, proceed with mutations and circularized sequences checks
    if (!high_speed && best_match == 0) {
      // Generate mutations for the original barcode
      auto mutations = generate_recursive_mutations_cpp(barcode, depth);
      for (auto& mutation : mutations) {
        if (wl.find(mutation) != wl.end()) {
          int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
          if (dl_distance <= maxDistance) {
            double poisson_score = wl[mutation].poisson_score;
            if (poisson_score > highest_poisson_score || 
              (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
              best_match = mutation;
              highest_poisson_score = poisson_score;
              best_dl_distance = dl_distance;
            }
          }
        }
      }
      // If still no match, proceed with circularized sequences
      if (best_match == 0) {
        auto circularizedSequences = kmer_circ_cpp(barcode, 16, 16, verbose);
        for (auto& circSeq : circularizedSequences) {
          auto circMutations = generate_recursive_mutations_cpp(circSeq, breadth);
          for (auto& mutation : circMutations) {
            if (wl.find(mutation) != wl.end()) {
              int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
              if (dl_distance <= maxDistance) {
                double poisson_score = wl[mutation].poisson_score;
                if (poisson_score > highest_poisson_score || 
                  (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
                  best_match = mutation;
                  highest_poisson_score = poisson_score;
                  best_dl_distance = dl_distance;
                }
              }
            }
          }
        }
      }
    }
    // Update results based on the best match found
    if (best_match != 0) {
      wl[best_match].count++; // Increment the count for the best match
      results[i] = static_cast<double>(best_match); // Convert the best match to double and assign
    } else {
      results[i] = NA_REAL; // Assign NA if no suitable match is found
    }
    
    if (verbose && best_match != 0) {
#pragma omp critical
      Rcpp::Rcout << "Best match found for barcode: " << barcode << " is: " << best_match 
                  << " with Poisson score: " << highest_poisson_score << " and DL distance: " << best_dl_distance << std::endl;
    }
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector correct_barcodes(Rcpp::NumericVector barcodes, bool verbose = false, 
  int nthreads = 1, int depth = 2, int breadth = 1, bool high_speed = true, int maxDistance = 2) {
  int n = barcodes.size();
  Rcpp::NumericVector results(n, NA_REAL); // Initialize results with NA_REAL
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads); // Set the number of threads
  
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    int64_t barcode = cpp_barcodes[i];
    double highest_poisson_score = -1.0;
    int best_dl_distance = std::numeric_limits<int>::max();
    int64_t best_match = 0;
    // Generate mutations for the original barcode as the first step
    auto mutations = generate_recursive_mutations_cpp(barcode, depth);
    for (auto& mutation : mutations) {
      if (wl.find(mutation) != wl.end()) {
        int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
        if (dl_distance <= maxDistance) {
          double poisson_score = wl[mutation].poisson_score;
          if (poisson_score > highest_poisson_score || 
            (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
            best_match = mutation;
            highest_poisson_score = poisson_score;
            best_dl_distance = dl_distance;
          }
        }
      }
    }
    // Generate circularized sequences if no match from mutations
    if (best_match == 0) {
      auto circularizedSequences = kmer_circ_cpp(barcode, 16, 16, verbose);
      for (auto& circSeq : circularizedSequences) {
        auto circMutations = generate_recursive_mutations_cpp(circSeq, breadth);
        for (auto& mutation : circMutations) {
          if (wl.find(mutation) != wl.end()) {
            int dl_distance = dl_dist_cpp({barcode}, {mutation}, maxDistance)[0];
            if (dl_distance <= maxDistance) {
              double poisson_score = wl[mutation].poisson_score;
              if (poisson_score > highest_poisson_score || 
                (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
                best_match = mutation;
                highest_poisson_score = poisson_score;
                best_dl_distance = dl_distance;
              }
            }
          }
        }
      }
    }
    // If high-speed mode is false and still no match, check against each whitelist entry
    if (!high_speed && best_match == 0) {
      for (const auto& wl_entry : wl) {
        int dl_distance = dl_dist_cpp({barcode}, {wl_entry.first}, maxDistance)[0];
        if (dl_distance >= 0 && dl_distance <= maxDistance) {
          double poisson_score = wl_entry.second.poisson_score;
          if (poisson_score > highest_poisson_score || 
            (poisson_score == highest_poisson_score && dl_distance < best_dl_distance)) {
            best_match = wl_entry.first;
            highest_poisson_score = poisson_score;
            best_dl_distance = dl_distance;
          }
        }
      }
    }
    // Update results based on the best match found
    if (best_match != 0) {
      wl[best_match].count++; // Increment the count for the best match
      results[i] = static_cast<double>(best_match); // Convert the best match to double and assign
    } else {
      results[i] = NA_REAL; // Assign NA if no suitable match is found
    }
    if (verbose && best_match != 0) {
#pragma omp critical
      Rcpp::Rcout << "Best match found for barcode: " << barcode << " is: " << best_match 
                  << " with Poisson score: " << highest_poisson_score << " and DL distance: " << best_dl_distance << std::endl;
    }
  }
  return results;
}
