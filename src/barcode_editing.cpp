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
  for (int64_t barcode : barcodes) {
      wl[barcode] = barcode_data{boost::none, boost::none};
  }
}

// [[Rcpp::export]]
void populate_whitelist(NumericVector barcodes, 
  NumericVector poisson_data = NumericVector(), 
  NumericVector counts = NumericVector()) {
  std::vector<int64_t> cpp_barcodes = as<std::vector<int64_t>>(barcodes);
  if (wl.empty()) {
    generate_whitelist(cpp_barcodes);
  }
  bool has_poisson_data = poisson_data.size() == barcodes.size();
  bool has_count_data = counts.size() == barcodes.size();
  for (size_t i = 0; i < cpp_barcodes.size(); ++i) {
    int64_t barcode = cpp_barcodes[i];
    auto it = wl.find(barcode);
    
    if (it != wl.end()) {
      // Existing entry found, update it carefully
      if (has_poisson_data && !Rcpp::NumericVector::is_na(poisson_data[i])) {
        it->second.poisson_score = static_cast<double>(poisson_data[i]); // Ensure correct type
      }
      if (has_count_data && !Rcpp::NumericVector::is_na(counts[i])) {
        it->second.count = static_cast<int>(counts[i]);
      }
    } else {
      // Entry not found, insert a new one
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
  std::vector<std::string> sequences;
  std::vector<int> counts;
  std::vector<double> poisson_scores;
  for (const auto& item : wl) {
    // Create a temporary vector with a single int64_t value for barcode conversion
    std::vector<int64_t> barcodeVec = {item.first};
    // Assuming a fixed sequence length, e.g., 16
    sequences.push_back(bits_to_sequence_cpp(barcodeVec, 16));
    counts.push_back(item.second.count ? *(item.second.count) : NA_INTEGER);
    poisson_scores.push_back(item.second.poisson_score ? *(item.second.poisson_score) : NA_REAL);
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
  wl.clear();
}

// [[Rcpp::export]]
Rcpp::StringVector correct_barcodes(Rcpp::NumericVector barcodes, bool verbose = false, int nthreads = 1) {
  int n = barcodes.size();
  Rcpp::StringVector results(n);
  std::vector<int64_t> cpp_barcodes = Rcpp::as<std::vector<int64_t>>(barcodes);
  omp_set_num_threads(nthreads); // Set the desired number of threads
  #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
      int64_t barcode = cpp_barcodes[i];
      std::string result;
      auto wl_it = wl.find(barcode);
      if (wl_it != wl.end()) {
        #pragma omp critical
      {
        if (wl_it->second.count) *(wl_it->second.count) += 1;
        else wl_it->second.count = boost::optional<int>(1);
      }
        std::vector<int64_t> barcodeVec = {barcode};
        result = "match:" + bits_to_sequence_cpp(barcodeVec, 16);
        if(verbose) {
            #pragma omp critical
        Rcpp::Rcout << "Direct match found for barcode: " << barcode << std::endl;
      }
    } else {
      std::unordered_set<int64_t> combined;
      std::vector<int64_t> circularized_sequences = kmer_circ_cpp(barcode, 16, 16, false);
      for (auto& circ_seq : circularized_sequences) {
        std::vector<int64_t> mutations = generate_recursive_mutations_cpp(circ_seq, 2);
        combined.insert(mutations.begin(), mutations.end());
      }
      std::vector<int64_t> original_mutations = generate_recursive_mutations_cpp(barcode, 3);
      combined.insert(original_mutations.begin(), original_mutations.end());
      double max_poisson = -1;
      int64_t best_seq_poisson = 0;
      for (const auto& seq : combined) {
        auto seq_it = wl.find(seq);
        if (seq_it != wl.end() && seq_it->second.poisson_score && *seq_it->second.poisson_score > max_poisson) {
          max_poisson = *seq_it->second.poisson_score;
          best_seq_poisson = seq;
        }
      }
      if (best_seq_poisson != 0) {
                #pragma omp critical
                            {
                              auto best_it = wl.find(best_seq_poisson);
                              if (best_it != wl.end()) {
                              if (best_it->second.count) *(best_it->second.count) += 1;
                              else best_it->second.count = boost::optional<int>(1);
                            }
                      }
        std::vector<int64_t> bestSeqVec = {best_seq_poisson};
        result = bits_to_sequence_cpp(bestSeqVec, 16);
        if(verbose) {
            #pragma omp critical
          Rcpp::Rcout << "Poisson score match found for barcode: " << bits_to_sequence_cpp(bestSeqVec, 16) << std::endl;
        }
      } else { 
        result = "no_match";
      }
    }
    results[i] = result;
  }
  return results;
}
  
