#include <Rcpp.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <functional>
#include <unordered_map>

using namespace Rcpp;
using boost::multi_index_container;
using namespace boost::multi_index;
using namespace std;

// Define tags for multi-indexing
struct id_tag {};
struct length_tag {};
struct type_tag {};
struct order_tag {};
struct direction_tag {};
struct global_class_tag {};

struct PositionInfo {
  std::tuple<std::string, std::string, std::string, std::string> position_data; // Tuple of primary and secondary start positions
};

struct ReadElement {
  std::string class_id;
  std::string seq; // Sequence information, there's nothing in here tbh
  std::string global_class;
  boost::optional<int> expected_length; // Optional expected length
  std::string type; // Type of the element (e.g., "static", "variable")
  int order; // Order in the layout
  std::string direction; // Direction (e.g., "forward", "reverse")
  boost::optional<std::tuple<int, int, int>> misalignment_threshold; // Optional misalignment threshold
  boost::optional<PositionInfo> position_data; // Optional position information for variable elements
  
  ReadElement(
    const std::string& class_id,
    const std::string& seq,
    const boost::optional<int>& expected_length,
    const std::string& type,
    int order,
    const std::string& direction,
    const std::string& global_class,
    const boost::optional<std::tuple<int, int, int>>& misalignment_threshold = boost::none,
    const boost::optional<PositionInfo>& position_data = boost::none
  ) : class_id(class_id),
  seq(seq),
  expected_length(expected_length),
  type(type),
  order(order),
  direction(direction),
  global_class(global_class),
  misalignment_threshold(misalignment_threshold),
  position_data(position_data) {}
};

typedef multi_index_container<
  ReadElement,
  indexed_by<
    ordered_unique<tag<id_tag>, member<ReadElement, std::string, &ReadElement::class_id>>,
    ordered_non_unique<tag<length_tag>, member<ReadElement, boost::optional<int>, &ReadElement::expected_length>>,
    ordered_non_unique<tag<type_tag>, member<ReadElement, std::string, &ReadElement::type>>,
    ordered_non_unique<tag<order_tag>, member<ReadElement, int, &ReadElement::order>>,
    ordered_non_unique<tag<direction_tag>, member<ReadElement, std::string, &ReadElement::direction>>,
    ordered_non_unique<tag<global_class_tag>, member<ReadElement, std::string, &ReadElement::global_class>>
  >
> ReadLayout;

ReadLayout prep_read_layout_cpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold) {
  ReadLayout container;
  // Assume 'query_id' matches 'class_id' from 'read_layout'
  std::unordered_map<std::string, std::tuple<int, int, int>> thresholdsMap;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for (int i = 0; i < query_id.size(); ++i) {
    double threshold = misal_threshold[i];
    double sd = std::ceil(misal_sd[i]);
    
    int lower_bound = static_cast<int>(std::floor(threshold - sd));
    int upper_bound = static_cast<int>(std::ceil(threshold + sd));
    
    thresholdsMap[Rcpp::as<std::string>(query_id[i])] = std::make_tuple(
      lower_bound,
      static_cast<int>(threshold),
      upper_bound
    );
  }
  
  Rcpp::StringVector class_id = read_layout["class_id"];
  Rcpp::StringVector class_column = read_layout["class"];
  Rcpp::StringVector seq = read_layout["seq"];
  Rcpp::IntegerVector expected_length = read_layout["expected_length"];
  Rcpp::StringVector type = read_layout["type"];
  Rcpp::IntegerVector order = read_layout["order"];
  Rcpp::StringVector direction = read_layout["direction"];
  
  for (int i = 0; i < class_id.size(); ++i) {
    std::string cid = Rcpp::as<std::string>(class_id[i]);
    std::string class_value = Rcpp::as<std::string>(class_column[i]); // Get the class value
    
    boost::optional<std::tuple<int, int, int>> misalignment_threshold_opt;
    if (thresholdsMap.count(cid)) {
      misalignment_threshold_opt = boost::optional<std::tuple<int, int, int>>(thresholdsMap[cid]);
    }else{
      misalignment_threshold_opt = boost::none;
    }
    
    boost::optional<int> expected_length_opt;
    if (class_value != "read") {
      expected_length_opt = expected_length[i];
    }
    
    container.insert(ReadElement(
        cid,
        Rcpp::as<std::string>(seq[i]),
        expected_length_opt,
        Rcpp::as<std::string>(type[i]),
        order[i],
        Rcpp::as<std::string>(direction[i]),
        Rcpp::as<std::string>(class_column[i]),
        misalignment_threshold_opt
    ));
  }
  return container;
}

struct SigElement {
  std::string class_id;
  std::string global_class;
  boost::optional<int> edit_distance; // Optional edit distance
  std::pair<int, int> position;
  std::string type;
  int order;
  std::string direction;
  boost::optional<bool> element_pass;
  
  SigElement(
    std::string class_id,
    std::string global_class,
    boost::optional<int> edit_distance,
    std::pair<int, int> position,
    std::string type,
    int order,
    std::string direction,
    boost::optional<bool> element_pass = boost::none
  ) : class_id(std::move(class_id)),
  global_class(std::move(global_class)),
  edit_distance(edit_distance),
  position(position),
  type(std::move(type)),
  order(order),
  direction(std::move(direction)),
  element_pass(element_pass) {}
};

struct sig_id_tag {};
struct sig_global_tag {};
struct sig_ed_tag {};
struct sig_dir_tag {};
struct sig_order_tag {};
struct sig_pass_tag {}; //added to be able to index by element_pass

// Define the multi-index container for SigString
typedef multi_index_container<
  SigElement,
  indexed_by<
    ordered_non_unique<tag<sig_id_tag>, member<SigElement, std::string, &SigElement::class_id>>,
    ordered_non_unique<tag<sig_global_tag>, member<SigElement, std::string, &SigElement::global_class>>,
    ordered_non_unique<tag<sig_ed_tag>, member<SigElement, boost::optional<int>, &SigElement::edit_distance>>,
    ordered_non_unique<tag<sig_order_tag>, member<SigElement, int, &SigElement::order>>,
    ordered_non_unique<tag<sig_dir_tag>, member<SigElement, std::string, &SigElement::direction>>,
    ordered_non_unique<tag<sig_pass_tag>, member<SigElement, boost::optional<bool>, &SigElement::element_pass>>
  >
> SigString;

struct StringInfo {
  int str_length;
  std::string str_type;
  std::string str_id;  // Make sure this is added
  
  StringInfo(int length = 0, std::string id = "NA", std::string type = "undefined")
    : str_length(length), str_type(std::move(type)), str_id(std::move(id)) {}
};

struct ReadData {
  std::unordered_map<std::string, StringInfo> string_info;
  SigString sigstring;
  
  void addStringInfo(const std::string& key, int length, const std::string& type, const std::string& id) {
    string_info[key] = StringInfo(length, type, id);
  }
};

struct stat_elem {
  int forw_pass = 0;
  int rev_pass = 0;
  int unique_forw_pass = 0; 
  int unique_rev_pass = 0;
  int term_forw_elements = 0;
  int term_rev_elements = 0;
};

struct PositionCalcFunc {
  std::function<int(const ReadData&)> primaryStartFunc;
  std::function<int(const ReadData&)> secondaryStartFunc;
  std::function<int(const ReadData&)> primaryStopFunc;
  std::function<int(const ReadData&)> secondaryStopFunc;
};

// Parse static elements to the left of a read
void parse_static_left(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions, bool static_mode) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto& direction_index = readLayout.get<direction_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::max_element(segment_orders.begin(), segment_orders.end());
  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    if(static_mode){
      if (it_var->type != "variable" || !it_var->expected_length || it_var->direction != direction) {
        continue;
      }
      if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
        continue;
      }
    } else {
      if(it_var -> global_class != "read"){
        continue;
      }
    }
    if (it_var->order <= read_order) {
      // Find previous element
      if (verbose) {
        Rcpp::Rcout << "Processing varElem: " << it_var->class_id << "\n";
      }
      if (it_var->order > 1) {
        auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
        std::string prevType = prevElem->type;
        std::string primaryStart;
        std::string primaryStop;
        std::string secondaryStart;
        std::string secondaryStop;
        
        if (prevElem != order_index.end()) {
          std::string primaryStart;
          primaryStart = prevElem->class_id + "|stop+1";
          
          if(static_mode){
            primaryStop = prevElem->class_id + "|stop+" + std::to_string(*it_var->expected_length);
          } else {
            auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
            primaryStop = nextElem->class_id + "|start-1";
          }
          // Check for chained variables
          if(static_mode){
            if (prevType == "variable" && assignedPositions.find(prevElem->class_id) != assignedPositions.end()) {
              auto& prevAssigned = assignedPositions[prevElem->class_id];
              primaryStart = primaryStart + "|chained_start";
              primaryStop = primaryStop + "|chained_stop";
             int read_order = *std::max_element(segment_orders.begin(), segment_orders.end());
              auto NextElem = readLayout.get<order_tag>().find(it_var->order+1);
                if(NextElem->type == "static"){
                  secondaryStart = NextElem->class_id + "|start-" + std::to_string(*NextElem->expected_length);
                  secondaryStop = NextElem->class_id + "|start-1";
                }
              // Assuming it_var is initialized correctly before this loop
              // while (it_var->order < read_order + 1) {
              //    auto NextElem = readLayout.get<order_tag>().find(it_var->order+1);
              //   if (verbose) {
              //     Rcpp::Rcout << "In the while loop!\n";
              //     Rcpp::Rcout << "Current class_id is " << it_var->class_id << "\n";
              //     Rcpp::Rcout << "Current order is " << it_var->order << "\n";
              //     Rcpp::Rcout << "Read order is " << read_order + 1 << "\n";
              //     Rcpp::Rcout << "Next class_id is " << NextElem->class_id << "\n";
              //   }
              //   if (NextElem->type == "static") {
              //     int sum_lengths = 0;
              //     auto tempIt = it_var;
              //     while (tempIt->order <= NextElem->order) {
              //       --tempIt; // Move backwards
              //       if (tempIt->type == "static" && tempIt->expected_length) {
              //         sum_lengths += *tempIt->expected_length;
              //       }
              //     }
              //     std::string secondaryStart;
              //     std::string secondaryStop;
              //     if (NextElem->expected_length) {
              //       secondaryStart = NextElem->class_id + "|start-" + std::to_string(sum_lengths + 1); 
              //       secondaryStop = NextElem->class_id + "|start-1";
              //     }
              //     if (verbose) {
              //       Rcpp::Rcout << "Secondary Start: " << secondaryStart << ", Secondary Stop: " << secondaryStop << "\n";
              //     }
              //   }
              //   ++it_var; // Always proceed to the next element to avoid infinite loop
              // }
                 else {
                secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
                secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
              }
            } else if (prevType == "static") {
              if(*prevElem->expected_length < 15 && prevElem->global_class != "poly_tail"){
                auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
                int offset_result = (*it_var->expected_length + 1);
                secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(offset_result);
                if(pre_prevElem->type == "variable"){
                  secondaryStart = secondaryStart + "|chained_start";
                  }
                } else {
              secondaryStart = prevElem->class_id + "|none|left_terminal_linked";
                  }
              secondaryStop = prevElem->class_id + "|none|left_terminal_linked";
            }
          } else {
            if (prevType == "variable" && assignedPositions.find(prevElem->class_id) != assignedPositions.end()) {
              auto& prevAssigned = assignedPositions[prevElem->class_id];
                secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
                secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
            } else if (prevType == "static") {
              if(*prevElem->expected_length < 15 && prevElem->global_class != "poly_tail"){
                auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
                int offset_result = (*it_var->expected_length + 1);
                secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(offset_result);
                if(pre_prevElem->type == "variable"){
                  secondaryStart = secondaryStart + "|chained_start";
                }
              } else {
                secondaryStart = prevElem->class_id + "|none|left_terminal_linked";
              }
              auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
              if(nextElem->global_class == "poly_tail"){
                auto skipElem = readLayout.get<order_tag>().find(it_var->order + 2); 
                // move one over from the poly tail if there is one
                secondaryStop = skipElem->class_id+"|start-1|poly_skipped";
              } else {
                secondaryStop = prevElem->class_id + "|none|left_terminal_linked";
              }
            }
          }
          // Update assignedPositions
          if(assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
            auto& existingEntry = assignedPositions[it_var->class_id];
            if(std::get<0>(existingEntry) != primaryStart){
              std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
            }
            if(std::get<2>(existingEntry) != primaryStop){
              std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
            }
          } else {
            assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
          }
          if(verbose) {
            Rcpp::Rcout << "Parse Left: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
            Rcpp::Rcout << "Parse Left: Primary Stop "  << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
          }
        }
      }
      if (verbose) {
        Rcpp::Rcout << "Parse left varElem: " << it_var->class_id << "\n";
      }
    }
  }
}

//Parse static elements to the right of a read
void parse_static_right(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions, bool static_mode) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::min_element(segment_orders.begin(), segment_orders.end());
  for (auto it_var = std::make_reverse_iterator(var_elements.second); it_var != std::make_reverse_iterator(var_elements.first); ++it_var) {
    if(static_mode){
      if (it_var->type != "variable" || !it_var->expected_length || it_var->direction != direction) {
        continue;
      }
      if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
        continue;
      }
    } else {
      if(it_var -> global_class != "read"){
        continue;
      }
    }
    if (it_var->order >= read_order) {
      
      if (verbose) {
        Rcpp::Rcout << "Processing varElem: " << it_var->class_id << "\n";
      }
      
      if (it_var->order < segment_orders.back()) {
        auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
        std::string nextType = nextElem->type;
        
        std::string primaryStart;
        std::string primaryStop;
        std::string secondaryStart;
        std::string secondaryStop;
        
        //added in to try to ameliorate the issue here w/short static sequences
        if (nextElem != order_index.end()) {
          if(static_mode){
            primaryStart = nextElem->class_id + "|start-" + std::to_string(*it_var->expected_length);
          } else {
            auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
            std::string prevType = prevElem->type;
            primaryStart = prevElem->class_id + "|stop+1";
          }
            primaryStop = nextElem->class_id + "|start-1";
          // Check for chained variables
          if(static_mode){
            if (nextType == "variable" && assignedPositions.find(nextElem->class_id) != assignedPositions.end()) {
              auto& nextAssigned = assignedPositions[nextElem->class_id];
              primaryStart = primaryStart + "|chained_start";
              primaryStop = primaryStop + "|chained_stop";

                auto prevElem = readLayout.get<order_tag>().find(it_var->order-1);
                if(prevElem->type == "static"){
                  secondaryStart = prevElem->class_id + "|stop+1";
                  secondaryStop = prevElem->class_id + "|stop+"+ std::to_string(*prevElem->expected_length);
                }
               else {
                secondaryStart = std::get<0>(nextAssigned) + "|chained_start";
                secondaryStop = std::get<2>(nextAssigned) + "|chained_stop";
              }
            } else if (nextType == "static") {
              secondaryStart = nextElem->class_id + "|none|right_terminal_linked";
              if(*nextElem->expected_length < 15 && nextElem->global_class != "poly_tail"){
                auto post_nextElem = readLayout.get<order_tag>().find(it_var->order + 2);
                int offset_result = (*it_var->expected_length + 1);
                secondaryStop = post_nextElem->class_id + "|stop+" + std::to_string(offset_result);
              } else {
                secondaryStop = nextElem->class_id + "|none|right_terminal_linked";
              }
            }
          } else {
            if (nextType == "variable" && assignedPositions.find(nextElem->class_id) != assignedPositions.end()) {
              auto& nextAssigned = assignedPositions[nextElem->class_id];
              if(it_var->order+1 < read_order){
                auto NextElem = readLayout.get<order_tag>().find(it_var->order+1);
                if(NextElem->type == "static"){
                  secondaryStart = NextElem->class_id + "stop+1";
                  secondaryStop = NextElem->class_id + "stop+"+ std::to_string(*NextElem->expected_length);
                }
              } else {
                secondaryStart = std::get<0>(nextAssigned) + "|chained_start";
                secondaryStop = std::get<2>(nextAssigned) + "|chained_stop";
              }
            } else if (nextType == "static") {
              auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
              if(prevElem->global_class == "poly_tail"){
                auto skipElem = readLayout.get<order_tag>().find(it_var->order-2);
                secondaryStart = skipElem->class_id+"|stop+1|poly_skipped";
              } else {
                if(*prevElem->expected_length < 15 && prevElem->global_class != "poly_tail"){
                  auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
                  int offset_result = (*it_var->expected_length + 1);
                  secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(offset_result);
                  } else {
                    secondaryStart = prevElem->class_id + "|none|left_terminal_linked";
                  }
              }
              auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
              if(nextElem->global_class == "poly_tail"){
                auto skipElem = readLayout.get<order_tag>().find(it_var->order + 2); // move one over from the poly tail if there is one
                secondaryStop = skipElem->class_id+"|start-1|poly_skipped";
              } else {
                if(*nextElem->expected_length < 15 && nextElem->global_class != "poly_tail"){
                  auto post_nextElem = readLayout.get<order_tag>().find(it_var->order + 2);
                  int offset_result = (*it_var->expected_length + 1);
                  secondaryStop = post_nextElem->class_id + "|stop+" + std::to_string(offset_result);
                  } else {
                    secondaryStop = nextElem->class_id + "|none|right_terminal_linked";
                }
              }
            }
          }
          if(assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
            auto& existingEntry = assignedPositions[it_var->class_id];
            if(std::get<0>(existingEntry) != primaryStart){
              std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
            }
            if(std::get<2>(existingEntry) != primaryStop){
              std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
            }
          } else {
            assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
          }
          if(verbose) {
            Rcpp::Rcout << "Parse Right: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
            Rcpp::Rcout << "Parse Right: Primary Stop "  << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
          }
        }
      }
      if (verbose) {
        Rcpp::Rcout << "Parse right varElem: " << it_var->class_id << "\n";
      }
    }
  }
}

// VarScan function
void VarScan(ReadLayout& readLayout, bool verbose) {
  std::unordered_map<std::string, std::tuple<std::string, std::string, std::string, std::string>> assignedPositions;
  auto& global_class_index = readLayout.get<global_class_tag>();
  auto& order_index = readLayout.get<order_tag>();
  auto read_range = global_class_index.equal_range("read");
  
  // Separate forward and reverse orders
  std::vector<int> forward_orders, reverse_orders;
  for (auto it_order = order_index.begin(); it_order != order_index.end(); ++it_order) {
    if (it_order->direction == "forward") {
      forward_orders.push_back(it_order->order);
    } else {
      reverse_orders.push_back(it_order->order);
    }
  }
  
  // Process each read separately
  for (auto it_read = read_range.first; it_read != read_range.second; ++it_read) {
    std::string read_direction = it_read->direction;
    int read_order = it_read->order;
    
    // Subdivide orders based on read position
    std::vector<int> left_vector, right_vector;
    if (read_direction == "forward") {
      auto read_order_here = std::find(forward_orders.begin(), forward_orders.end(), read_order);
      left_vector = std::vector<int>(forward_orders.begin(), read_order_here);
      right_vector = std::vector<int>(read_order_here, forward_orders.end());
    } else {
      auto read_order_here = std::find(reverse_orders.begin(), reverse_orders.end(), read_order);
      left_vector = std::vector<int>(reverse_orders.begin(), read_order_here);
      right_vector = std::vector<int>(read_order_here, reverse_orders.end());
    }
    
    //Parse static elements
    parse_static_left(readLayout, left_vector, read_direction, verbose, assignedPositions, true);
    parse_static_right(readLayout, right_vector, read_direction, verbose, assignedPositions, true);
    //Parse read elements
    parse_static_left(readLayout, left_vector, read_direction, verbose, assignedPositions, false);
    parse_static_right(readLayout, right_vector, read_direction, verbose, assignedPositions, false);
  }
  for (const auto& entry : assignedPositions) {
    const auto& class_id = entry.first;
    const auto& positions = entry.second;
    // Retrieve the ReadElement from readLayout using class_id
    auto& id_index = readLayout.get<id_tag>();
    auto it = id_index.find(class_id);
    if (it != id_index.end()) {
      // Update the position_data field
      PositionInfo posInfo;
      posInfo.position_data = positions;
      id_index.modify(it, [&posInfo](ReadElement& elem) {
        elem.position_data = posInfo;
      });
    }
  }
  if (verbose) {
    Rcpp::Rcout << "Final Assigned Positions:\n";
    for (const auto& entry : assignedPositions) {
      const auto& class_id = entry.first;
      const auto& positions = entry.second;
      Rcpp::Rcout << "Class ID: " << class_id << "\n"
                  << "  Primary Start: " << std::get<0>(positions) << "\n"
                  << "  Secondary Start: " << std::get<1>(positions) << "\n"
                  << "  Primary Stop: " << std::get<2>(positions) << "\n"
                  << "  Secondary Stop: " << std::get<3>(positions) << "\n";
    }
  }
}

//StatCounts from a sigstring, tabulate forward and reverse elements
stat_elem StatCounts(const SigString& sigstring) {
  stat_elem counts;
  std::set<std::string> seen_forw_pass_ids;
  std::set<std::string> seen_rev_pass_ids;
  for (const auto& element : sigstring) {
    if (element.global_class == "poly_tail" || element.type != "static") continue;
    if (element.direction == "forward" && element.element_pass && *element.element_pass) {
      ++counts.forw_pass;
      // Increment unique counter if it's a new ID
      if (seen_forw_pass_ids.insert(element.class_id).second) {
        ++counts.unique_forw_pass;
      }
    } else if (element.direction == "reverse" && element.element_pass && *element.element_pass) {
      ++counts.rev_pass;
      // Increment unique counter if it's a new ID
      if (seen_rev_pass_ids.insert(element.class_id).second) {
        ++counts.unique_rev_pass;
      }
    }
  }
 for(const auto& element : sigstring){
   if ((element.global_class == "forw_primer" || element.global_class == "rev_primer") && element.element_pass && *element.element_pass) {
     if (element.direction == "forward") {
       ++counts.term_forw_elements;
     } else if (element.direction == "reverse") {
       ++counts.term_rev_elements;
     }
   }
 }
  return counts;
}

//concatenation scanning
void concat_scan(ReadData& readData, bool verbose) {
  stat_elem counts = StatCounts(readData.sigstring);
    if(verbose){
      Rcpp::Rcout << "Forward pass: " << counts.forw_pass << "\nReverse pass: " << counts.rev_pass << "\n";
      Rcpp::Rcout << "Forward UNIQUE pass: " << counts.unique_forw_pass << "\nReverse UNIQUE pass: " << counts.unique_rev_pass << "\n";
      Rcpp::Rcout << "Forward terminal pass: " << counts.term_forw_elements << "\nReverse terminal pass: " << counts.term_rev_elements << "\n";
   }
  if(counts.forw_pass == 0 && counts.rev_pass > 0){
    readData.string_info["type"].str_type = "R";
  } else if (counts.rev_pass == 0 && counts.forw_pass > 0){
    readData.string_info["type"].str_type = "F";
  }
  if(counts.forw_pass != counts.unique_forw_pass && counts.rev_pass == 0){
    readData.string_info["type"].str_type = "F_F";
  } else if (counts.rev_pass != counts.unique_rev_pass && counts.forw_pass == 0){
    readData.string_info["type"].str_type = "R_R";
  }
  if(counts.forw_pass > 0 && counts.rev_pass > 0){
    readData.string_info["type"].str_type = "FR_RF";
  }
  if(verbose){
    Rcpp::Rcout << "This sigstring is a type " << readData.string_info["type"].str_type << "!\n";
    Rcpp::Rcout << "This read is " << readData.string_info["type"].str_length << " bases long!\n";
  }
}

void concat_solve_hybrid(ReadData& readData, bool verbose) {
  auto& order_index = readData.sigstring.get<sig_order_tag>();
  for (auto it = order_index.begin(); it != order_index.end(); ++it) {
    if (it->global_class == "read") {
      int newStartPos = it->position.first;
      int newStopPos = it->position.second;
      if (newStartPos == -1) {
        auto preceding = it;
        if (preceding != order_index.begin()) {
          do {
            --preceding;
            if ((preceding->element_pass && *preceding->element_pass) || preceding->global_class == "poly_tail") {
              newStartPos = preceding->position.second + 1;
              break;
            }
          } while (preceding != order_index.begin());
        }
        if (newStartPos == -1) newStartPos = 1;
      }
      if (newStopPos == -1) {
        auto succeeding = it;
        ++succeeding;
        while (succeeding != order_index.end()) {
          if ((succeeding->element_pass && *succeeding->element_pass)||succeeding->global_class == "poly_tail") {
            newStopPos = succeeding->position.first - 1;
            break;
          }
          ++succeeding; 
        }
        if (newStopPos == -1) newStopPos = readData.string_info["type"].str_length;
      }
      // Apply the position updates
      order_index.modify(it, [newStartPos, newStopPos](SigElement& elem) {
        elem.position.first = newStartPos;
        elem.position.second = newStopPos;
      });
      if (verbose) {
        Rcpp::Rcout << "Read " << it->class_id << " adjusted to start: " << newStartPos << ", stop: " << newStopPos << "\n";
      }
    }
  }
}

std::vector<std::string> concat_solve_parallel(ReadData& readData, bool verbose) {
  std::string readType = readData.string_info["type"].str_type;
  std::vector<SigElement> elements;
  for (const auto& elem : readData.sigstring) {
    if(readType == "F_F"){
      if (elem.direction == "forward" && elem.type == "static") {
        elements.push_back(elem);
      }  
    } else if(readType == "R_R"){
      if(elem.direction == "reverse" && elem.type == "static") {
        elements.push_back(elem);
      }
    }
  }
  // First sort by order
  std::sort(elements.begin(), elements.end(),
    [](const SigElement& a, const SigElement& b) {
      return a.order < b.order;
    });
  // Then stable sort by start position
  std::stable_sort(elements.begin(), elements.end(),
    [](const SigElement& a, const SigElement& b) {
      return a.position.first < b.position.first;
    });
  // Detect where order decreases in next adjacent element
  std::vector<int> decreasePoints; // Indices where a decrease in order occurs
  for (size_t i = 0; i < elements.size() - 1; ++i) {
    if (elements[i].order > elements[i + 1].order) {
      // Detected a decrease in order, mark this point
      decreasePoints.push_back(i);
    }
  }
  // Split into two new sigstrings based on the first decrease point
  std::string sigString1, sigString2;
  int splitPoint = decreasePoints.empty() ? elements.size() : decreasePoints[0] + 1;
  // Build the first sigstring
  for (int i = 0; i < splitPoint; ++i) {
    if (i > 0) sigString1 += "|";
    sigString1 += elements[i].class_id + ":" + std::to_string(elements[i].edit_distance.value_or(0)) + ":" + std::to_string(elements[i].position.first) + ":" + std::to_string(elements[i].position.second);
  }
  // Append the new length position for the first string
  if (!elements.empty() && splitPoint > 0) {
    sigString1 += "<" + std::to_string(elements[splitPoint - 1].position.second) + ":undecided>";
  }
  // Build the second sigstring
  for (size_t i = splitPoint; i < elements.size(); ++i) {
    if (i > splitPoint) sigString2 += "|";
    sigString2 += elements[i].class_id + ":" + std::to_string(elements[i].edit_distance.value_or(0)) + ":" + std::to_string(elements[i].position.first) + ":" + std::to_string(elements[i].position.second);
  }
  // Append the new length position for the second string
  if (!elements.empty() && splitPoint < elements.size()) {
    sigString2 += "<" + std::to_string(elements.back().position.second) + ":undecided>";
  }
  // Verbose output to demonstrate the two new sigstrings
  if (verbose) {
    Rcpp::Rcout << "First new sigstring: " << sigString1 << "\n";
    Rcpp::Rcout << "Second new sigstring: " << sigString2 << "\n";
  }
  // Return the vector of sigstrings
  return {sigString1, sigString2};
}

using PositionFuncMap = std::unordered_map<std::string, PositionCalcFunc>;
PositionFuncMap createPositionFunctionMap(const ReadLayout& readLayout, bool verbose) {
    PositionFuncMap funcMap;
    for (const auto& element : readLayout) {
        if (element.type == "variable" && element.position_data) {
            const auto& positionData = element.position_data->position_data;
            auto primaryStart = std::get<0>(positionData);
            auto secondaryStart = std::get<1>(positionData);
            auto primaryStop = std::get<2>(positionData);
            auto secondaryStop = std::get<3>(positionData);
            auto parsePositionData = [verbose](const std::string& varElemClassId, const std::string& data, bool isPrimary) {
                std::stringstream ss(data);
                std::string refClassId, operation, additionalInfo;
                int offset = 0;
                std::getline(ss, refClassId, '|');
                std::getline(ss, operation, '|');
                if (!isPrimary) {
                    std::getline(ss, additionalInfo, '|'); // Only for secondary positions
                }
                if (verbose) {
                    Rcpp::Rcout << "Processing position data for variable element: " << varElemClassId << "\n";
                    Rcpp::Rcout << "  Reference class ID: " << refClassId << "\n";
                    Rcpp::Rcout << "  Operation: " << operation << "\n";
                    Rcpp::Rcout << "  Additional Information, if available: " << additionalInfo << "\n";
                }
                if (operation.find("+") != std::string::npos || operation.find("-") != std::string::npos) {
                    std::string offsetStr = operation.substr(operation.find_first_of("+-"));
                    offset = std::stoi(offsetStr);
                }
      return [varElemClassId, refClassId, offset, additionalInfo, &isPrimary, 
        operation, verbose](const ReadData& readData) -> int {
                    if (verbose) {
                        Rcpp::Rcout << "Made it into the lambda function!\n";
                    }

                    const auto& id_index = readData.sigstring.get<sig_id_tag>();
                    const auto& order_index = readData.sigstring.get<sig_order_tag>();
                    auto refClass = id_index.find(refClassId);
                    auto it = id_index.equal_range(refClassId).first;

                    if (it != id_index.end() && 
                        ((refClass->type == "static" && refClass->element_pass && *refClass->element_pass) || 
                        (refClass->type == "variable") || (refClass->global_class == "poly_tail"))) {
                        if (verbose) {
                            Rcpp::Rcout << "Made it into the loop and past the if statements!\n";
                        }
                        int startPos = it->position.first;
                        int stopPos = it->position.second;
                        int basePos = (operation.find("start") != std::string::npos) ? it->position.first : it->position.second;
                        int calculatedPos = basePos + offset;
                        if (verbose) {
                            Rcpp::Rcout << "  Variable Element: " << varElemClassId << "\n";
                            Rcpp::Rcout << "  Reference Element: " << refClassId << ", Base Position: " << basePos << "\n";
                            Rcpp::Rcout << "  Positional Information from start to stop: " << startPos << " and " << stopPos << "\n";
                            Rcpp::Rcout << "  Calculated Position: " << calculatedPos << "\n";
                        }
                        return calculatedPos;
                    } else if (!isPrimary && !additionalInfo.empty()) {
                      std::string readType = readData.string_info.at("type").str_type;
                      auto varClass = id_index.find(varElemClassId);
                      if(verbose){
                        string result = isPrimary ? "Primary" : "Secondary";
                        Rcpp::Rcout << "Flagged for some failures, checking which one it is!\n";
                        Rcpp::Rcout << "The read type is " << result << "\n";
                        Rcpp::Rcout << additionalInfo << " is the additional info!\n";
                      }
                      if (additionalInfo == "left_terminal_linked") {
                        if(verbose){
                          Rcpp::Rcout << "Detected a left_terminal linkage!\n";
                          Rcpp::Rcout << "Global class is " << varClass->global_class << "!\n";
                        }
                        if(varClass->global_class == "read"||varClass->global_class == "barcode"){
                          auto prevClass = readData.sigstring.get<sig_order_tag>().find(varClass->order - 1);
                          
                          int rel_pos_failure = static_cast<int>((1-(readData.string_info.at("type").str_length - prevClass->position.second) / static_cast<double>(readData.string_info.at("type").str_length)) * 100);
                          
                          int rel_edit_failure = 0;
                          double edit_distance = static_cast<double>(*prevClass->edit_distance);
                          int sequence_length = prevClass->position.second - prevClass->position.first + 1;
                          rel_edit_failure = static_cast<int>((1.0 - (edit_distance / sequence_length)) * 100);                                         if(rel_edit_failure >= 70 && verbose){
                            int new_terminal_pos = prevClass->position.second +1;
                          if(readData.string_info.at("type").str_type == "F" || readData.string_info.at("type").str_type == "R"){
                              Rcpp::Rcout << rel_pos_failure << " is the percent delta between EOS and the SNE!\n";
                              Rcpp::Rcout << rel_edit_failure << " is the percent error of the next element!\n";
                              Rcpp::Rcout << prevClass->edit_distance << " is the edit distance of the previous ordered element!\n";
                              Rcpp::Rcout << "Going to set the read start to the beginning of the next element (< 5% of read)!\n";
                              Rcpp::Rcout << "Setting new start to " << new_terminal_pos <<"!\n";
                            }
                            return new_terminal_pos;
                          } else
                          if(readData.string_info.at("type").str_type != "F" || readData.string_info.at("type").str_type != "R"){
                            return -1;
                          }
                            return 1;
                          } else {
                            return 0;
                          }
                      } else if (additionalInfo == "right_terminal_linked") {
                        if(verbose){
                          "Detected a right_terminal linkage!";
                          Rcpp::Rcout << "Global class is " << varClass->global_class << "!\n";
                        }
                        if(varClass->global_class == "read"||varClass->global_class == "barcode"){
                          auto nextClass = readData.sigstring.get<sig_order_tag>().find(varClass->order + 1);
                          int rel_pos_failure = static_cast<int>(((readData.string_info.at("type").str_length - nextClass->position.first) / static_cast<double>(readData.string_info.at("type").str_length)) * 100);
                          int rel_edit_failure = 0;
                          double edit_distance = static_cast<double>(*nextClass->edit_distance);
                          int sequence_length = nextClass->position.second - nextClass->position.first + 1;
                          rel_edit_failure = static_cast<int>((1.0 - (edit_distance / sequence_length)) * 100);
                        if(rel_edit_failure >= 70){
                          int new_terminal_pos = nextClass->position.first - 1;
                            if(verbose){
                              Rcpp::Rcout << rel_pos_failure << " is the percent delta between EOS and the SNE!\n";
                              Rcpp::Rcout << rel_edit_failure << " is the percent error of the next element!\n";
                              Rcpp::Rcout << nextClass->edit_distance << " is the edit distance of the next ordered element!\n";
                              Rcpp::Rcout << "Going to set the read end to the beginning of the next element (< 5% of read)!\n";
                              Rcpp::Rcout << "Setting new end to " << new_terminal_pos <<"!\n";
                            }
                            return new_terminal_pos;
                          } else
                          if(readData.string_info.at("type").str_type == "F" || readData.string_info.at("type").str_type == "R"){
                            if(readData.string_info.find("type") != readData.string_info.end()){
                              return readData.string_info.at("type").str_length; // End of the sequence
                            }
                        } else {
                          return -2;
                          }
                        }
                      }
                    }
                    return 0; // Default return value if conditions are not met
                };
            };
            funcMap[element.class_id] = {
                parsePositionData(element.class_id, primaryStart, true),
                parsePositionData(element.class_id, secondaryStart, false),
                parsePositionData(element.class_id, primaryStop, true),
                parsePositionData(element.class_id, secondaryStop, false)
            };
        }
    }
    if (verbose) {
        Rcpp::Rcout << "Position function map created with " << funcMap.size() << " entries.\n";
    }
    return funcMap;
}

void concat_solve_hybrid_v2(ReadData& readData, const PositionFuncMap& positionFuncMap, bool verbose) {
  if(verbose){
    Rcpp::Rcout << "Solving a hybrid concatenate!\n";
  }
  for (auto& element : readData.sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class != "read") {
      if(verbose) {
        Rcpp::Rcout << "Working on resolving: " << element.class_id << "\n";
        Rcpp::Rcout << "Now making the first pass...\n";
      }
      auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
      const auto& positionFuncs = positionFuncMap.at(element.class_id);
      int newStartPos = positionFuncs.primaryStartFunc(readData);
      if(newStartPos <= 0){
        newStartPos = positionFuncs.secondaryStartFunc(readData);
      }
      int newStopPos = positionFuncs.primaryStopFunc(readData);
      if(newStopPos <= 0){
        newStopPos = positionFuncs.secondaryStopFunc(readData);
      }
      
      if(newStartPos == -1 || newStopPos == -1){
        auto prevClass = readData.sigstring.get<sig_order_tag>().find(element.order - 1);
        if(verbose) {
          Rcpp::Rcout << "The previous class in order is " << prevClass->class_id << "\n";
          Rcpp::Rcout << "Now moving to the second pass...\n";
        }
        auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
        bool tmp_edit_pass = true;
        
        sig_id_index.modify(sig_id_index.iterator_to(*prevClass), [tmp_edit_pass](SigElement& elem) {
          elem.element_pass = tmp_edit_pass;
        });
        
        std::string element_pass = prevClass->element_pass ? "PASS\n" : "FAIL\n";
        
        if(verbose) {
          Rcpp::Rcout << "The previous class in order is " << prevClass->class_id << "\n";
          Rcpp::Rcout << "Testing to see if element pass works: " << element_pass;
          Rcpp::Rcout << "The new primary Start Pos of this is " << newStartPos << ".\n";
        }
        
        newStartPos = positionFuncs.primaryStartFunc(readData);
        if(newStartPos <= 0){
          newStartPos = positionFuncs.secondaryStartFunc(readData);
        }
        
        int newStopPos = positionFuncs.primaryStopFunc(readData);
        if(newStopPos <= 0){
          newStopPos = positionFuncs.secondaryStopFunc(readData);
        }
        
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });

        if(verbose){
          Rcpp::Rcout << "New pretend start pos: " << newStartPos << " to " << newStopPos << "\n";
        }
        
        bool new_edit_pass = false;
        sig_id_index.modify(sig_id_index.iterator_to(*prevClass), [new_edit_pass](SigElement& elem) {
          elem.element_pass = new_edit_pass;
        });
        if (verbose) {
          Rcpp::Rcout << "Updated positions for variable class_id: " << element.class_id
                      << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
        }} else
      if(newStartPos == -2 || newStopPos == -2){
        auto nextClass = readData.sigstring.get<sig_order_tag>().find(element.order + 1);
        if(verbose) {
          Rcpp::Rcout << "The previous class in order is " << nextClass->class_id << "\n";
          Rcpp::Rcout << "Now making the second pass...\n";
        }
        auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
        bool tmp_edit_pass = true;
        sig_id_index.modify(sig_id_index.iterator_to(*nextClass), [tmp_edit_pass](SigElement& elem) {
          elem.element_pass = tmp_edit_pass;
        });
        newStartPos = positionFuncs.primaryStartFunc(readData);
        std::string element_pass = nextClass->element_pass ? "pass\n" : "fail\n";
        
        if(verbose) {
          Rcpp::Rcout << "The previous class in order is " << nextClass->class_id << "\n";
          Rcpp::Rcout << "Testing to see if element pass works: " << element_pass;
          Rcpp::Rcout << "The new primary Start Pos of this is " << newStartPos << ".\n";
        }
        if(newStartPos <= 0){
          newStartPos = positionFuncs.secondaryStartFunc(readData);
        }
        int newStopPos = positionFuncs.primaryStopFunc(readData);
        if(newStopPos <= 0){
          newStopPos = positionFuncs.secondaryStopFunc(readData);
        }
        
        if(verbose){
          Rcpp::Rcout << "New start & stop pos pretending like your edit distance is 0: " << newStartPos << " to " << newStopPos << "\n";
        }
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
        
        bool new_edit_pass = false;
        sig_id_index.modify(sig_id_index.iterator_to(*nextClass), [new_edit_pass](SigElement& elem) {
          elem.element_pass = new_edit_pass;
        });
      }
      
      else {
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
      }
    } else {
      continue;
    }
  }
  for (auto& element : readData.sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class == "read") {
      if(verbose) {
        Rcpp::Rcout << "Dealing with everything that's IS a read!\n";
      }
      const auto& positionFuncs = positionFuncMap.at(element.class_id);
      int newStartPos = positionFuncs.primaryStartFunc(readData);
      if(newStartPos <= 0){
        newStartPos = positionFuncs.secondaryStartFunc(readData);
      }
      int newStopPos = positionFuncs.primaryStopFunc(readData);
      if(newStopPos <= 0){
        newStopPos = positionFuncs.secondaryStopFunc(readData);
      }
      //once for left_terminal linkages...and one for right.......yeah dude just copy paste code if it works the first time it's gotta work the n'th time amirite you moron
      if(newStartPos == -1 || newStopPos == -1){
        auto prevClass = readData.sigstring.get<sig_order_tag>().find(element.order - 1);
        if(verbose) {
          Rcpp::Rcout << "The previous class in order is " << prevClass->class_id << "\n";
          Rcpp::Rcout << "Now moving to the second pass...\n";
        }
        auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
        bool tmp_edit_pass = true;
        
        sig_id_index.modify(sig_id_index.iterator_to(*prevClass), [tmp_edit_pass](SigElement& elem) {
          elem.element_pass = tmp_edit_pass;
        });
        
        std::string element_pass = prevClass->element_pass ? "PASS\n" : "FAIL\n";
        
        if(verbose) {
          Rcpp::Rcout << "The previous class in order is " << prevClass->class_id << "\n";
          Rcpp::Rcout << "Testing to see if element pass works: " << element_pass;
          Rcpp::Rcout << "The new primary Start Pos of this is " << newStartPos << ".\n";
        }
        
        newStartPos = positionFuncs.primaryStartFunc(readData);
        if(newStartPos <= 0){
          newStartPos = positionFuncs.secondaryStartFunc(readData);
        }
        
        int newStopPos = positionFuncs.primaryStopFunc(readData);
        if(newStopPos <= 0){
          newStopPos = positionFuncs.secondaryStopFunc(readData);
        }
        
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
        
        if(verbose){
          Rcpp::Rcout << "New pretend start pos: " << newStartPos << " to " << newStopPos << "\n";
        }
        
        bool new_edit_pass = false;
        sig_id_index.modify(sig_id_index.iterator_to(*prevClass), [new_edit_pass](SigElement& elem) {
          elem.element_pass = new_edit_pass;
        });
        if (verbose) {
          Rcpp::Rcout << "Updated positions for variable class_id: " << element.class_id
                      << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
        }} else
      if(newStartPos == -2 || newStopPos == -2){
            auto nextClass = readData.sigstring.get<sig_order_tag>().find(element.order + 1);
            if(verbose) {
              Rcpp::Rcout << "The previous class in order is " << nextClass->class_id << "\n";
              Rcpp::Rcout << "Now making the second pass...\n";
            }
            auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
            bool tmp_edit_pass = true;
            sig_id_index.modify(sig_id_index.iterator_to(*nextClass), [tmp_edit_pass](SigElement& elem) {
              elem.element_pass = tmp_edit_pass;
            });
            newStartPos = positionFuncs.primaryStartFunc(readData);
            std::string element_pass = nextClass->element_pass ? "pass\n" : "fail\n";
            
            if(verbose) {
              Rcpp::Rcout << "The previous class in order is " << nextClass->class_id << "\n";
              Rcpp::Rcout << "Testing to see if element pass works: " << element_pass;
              Rcpp::Rcout << "The new primary Start Pos of this is " << newStartPos << ".\n";
            }
            if(newStartPos <= 0){
              newStartPos = positionFuncs.secondaryStartFunc(readData);
            }
            int newStopPos = positionFuncs.primaryStopFunc(readData);
            if(newStopPos <= 0){
              newStopPos = positionFuncs.secondaryStopFunc(readData);
            }
            
            if(verbose){
              Rcpp::Rcout << "New start & stop pos pretending like your edit distance is 0: " << newStartPos << " to " << newStopPos << "\n";
            }
            sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
              elem.position.first = newStartPos;
              elem.position.second = newStopPos;
            });
            
            bool new_edit_pass = false;
            sig_id_index.modify(sig_id_index.iterator_to(*nextClass), [new_edit_pass](SigElement& elem) {
              elem.element_pass = new_edit_pass;
            });
          }
      else {
        auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
            sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
              elem.position.first = newStartPos;
              elem.position.second = newStopPos;
            });
        }
    }
  }
}

//parsing a sigstring and adding stuff to it
void fillSigString(ReadData& readData, const ReadLayout& readLayout, 
  const std::string& signature, const PositionFuncMap& positionFuncMap, bool verbose) {
  //adding information to "type"
  auto read_info_pos = signature.find_last_of('<');
  std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2); // Remove '<' and '>'
  std::stringstream lengthTypeStream(lengthTypeStr);
  std::string lengthStr, id, type;
  std::getline(lengthTypeStream, lengthStr, ':');
  std::getline(lengthTypeStream, id, ':');
  std::getline(lengthTypeStream, type);
  
  int length = std::stoi(lengthStr);
  if(verbose) {
    Rcpp::Rcout << "Length is " << length << "\n";
    Rcpp::Rcout << "ID is " << id << "\n";
    Rcpp::Rcout << "Type is " << type << "\n";
  }  
  // Populate StringInfo with length and type
  readData.addStringInfo("type", length, id, type);

  // Process the signature and populate sigstring
  std::stringstream ss(signature);
  std::string token;
  auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
  while (std::getline(ss, token, '|')) {
    std::stringstream tokenStream(token);
    std::string id;
    int editDistance, startPos, endPos;
    
    std::getline(tokenStream, id, ':');
    tokenStream >> editDistance;
    tokenStream.ignore(1); // Ignore the colon
    tokenStream >> startPos;
    tokenStream.ignore(1); // Ignore the colon
    tokenStream >> endPos;
    
    auto& class_id_index = readLayout.get<id_tag>();
    auto it = class_id_index.find(id);
    if (it != class_id_index.end()) {
      // Calculate misalignment threshold if available
      boost::optional<bool> element_pass;
      if (it->misalignment_threshold && it->type == "static") {
        int threshold_floor = std::get<0>(*it->misalignment_threshold);
        element_pass = (editDistance <= threshold_floor);
      }
      // Create SigElement with complete information
      SigElement element(
          id,
          it->global_class,
          editDistance,
          std::make_pair(startPos, endPos),
          it->type,
          it->order,
          it->direction,
          element_pass
      );
      readData.sigstring.insert(std::move(element));
    }
  }
  for (const auto& readElement : readLayout) {
    if (sig_id_index.find(readElement.class_id) == sig_id_index.end()) {
      SigElement element(
          readElement.class_id,
          readElement.global_class,
          0, // default edit distance
          std::make_pair(-1, -1), // default positions
          readElement.type,
          readElement.order,
          readElement.direction
      );
      readData.sigstring.insert(std::move(element));
    }
  }
}

//displaying a sigstring
void displayOrderedSigString(const ReadData& readData) {
  auto printOrderedDirection = [&](const std::string& direction) {
    // Filter by direction and sort by order
    std::vector<SigElement> sorted_elements;
    auto& dir_index = readData.sigstring.get<sig_dir_tag>();
    auto dir_range = dir_index.equal_range(direction);
    for (auto it = dir_range.first; it != dir_range.second; ++it) {
      sorted_elements.push_back(*it);
    }
    std::sort(sorted_elements.begin(), sorted_elements.end(),
      [](const SigElement& a, const SigElement& b) { return a.order < b.order; });
    
    // Print elements
    Rcpp::Rcout << direction << " sequence layout:\n";
    for (const auto& element : sorted_elements) {
      std::string element_type = (element.type == "static") ? "static" : "variable";
      std::string pass_fail;
      if(element.type == "static" && element.global_class != "poly_tail"){
        pass_fail = (element.element_pass == true) ? ":pass" : ":fail";
        if(pass_fail == ":fail"){
          continue;
        }
      } else {
        pass_fail = "";
      }
      if(element.position.first == 0 || element.position.second == 0){
        continue;
      }
      std::string position_info = "|start:" + std::to_string(element.position.first) +
        "|stop:" + std::to_string(element.position.second);
      Rcpp::Rcout << "[" << element_type << ":" << element.class_id << pass_fail
                  << position_info << "]";
    }
    Rcpp::Rcout << "\n";
  };
  if(readData.string_info.at("type").str_type == "F"){
    printOrderedDirection("forward");
  } else if(readData.string_info.at("type").str_type == "R"){
    printOrderedDirection("reverse");
  } else {
    printOrderedDirection("forward");
    printOrderedDirection("reverse");
  }
}

std::string generateSigString(const ReadData& readData) {
  std::stringstream sigStringStream;
  // Function to append elements to the string stream based on direction
  auto appendElements = [&](const std::string& direction) {
    std::vector<SigElement> sorted_elements;
    auto& dir_index = readData.sigstring.get<sig_dir_tag>();
    auto dir_range = dir_index.equal_range(direction);
    for (auto it = dir_range.first; it != dir_range.second; ++it) {
      if (it->position.first <= 0 || it->position.second <= 0) continue; // Skip uninitialized elements
      sorted_elements.push_back(*it);
    }
    std::sort(sorted_elements.begin(), sorted_elements.end(),
      [](const SigElement& a, const SigElement& b) { return a.order < b.order; });
    for (size_t i = 0; i < sorted_elements.size(); ++i) {
      const auto& element = sorted_elements[i];
      if (i > 0) sigStringStream << "|";
      sigStringStream << element.class_id << ":" << element.edit_distance.value_or(0)
                      << ":" << element.position.first << ":" << element.position.second;
    }
  };
  // Process both forward and reverse elements
  appendElements("forward");
  if (readData.string_info.at("type").str_type == "R" || readData.string_info.at("type").str_type == "FR_RF") {
    if (!sigStringStream.str().empty()) sigStringStream << "+"; // Separator between directions if needed
    appendElements("reverse");
  }
  // Append the final part with length and type
  const auto& typeInfo = readData.string_info.at("type");
  sigStringStream << "<" << typeInfo.str_length << ":" << typeInfo.str_id << ":"<< typeInfo.str_type << ">";
  return sigStringStream.str();
}

void pos_scan(ReadData& readData, const PositionFuncMap& positionFuncMap, const ReadLayout& readLayout, bool verbose) {
  std::string readType = readData.string_info["type"].str_type;
  auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
  if(readType == "F" || readType == "R"){
    for (auto it = readData.sigstring.begin(); it != readData.sigstring.end(); ) {
      if ((readType == "F" && it->direction != "forward") || 
        (readType == "R" && it->direction != "reverse")) {
        it = readData.sigstring.erase(it); // Erase and move to the next element
      } else {
        ++it; // Move to the next element
      }
    }
    for (auto& element : readData.sigstring) {
      if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class != "read") {
        if(verbose) {
          Rcpp::Rcout << "Dealing with everything that's not a read!\n";
          Rcpp::Rcout << "Working on " << element.class_id << "\n";
        }
        const auto& positionFuncs = positionFuncMap.at(element.class_id);
        int newStartPos = positionFuncs.primaryStartFunc(readData);
        if(verbose){
          Rcpp::Rcout << newStartPos << "\n";
        }
        if(newStartPos <= 0){
          if(verbose){
            Rcpp::Rcout << "Triggered secondary function:\n";
          }
          newStartPos = positionFuncs.secondaryStartFunc(readData);
        }
        int newStopPos = positionFuncs.primaryStopFunc(readData);
        if(verbose){
          Rcpp::Rcout << newStopPos << "\n";
        }
        if(newStopPos <= 0){ 
          if(verbose){
            Rcpp::Rcout << "Triggered secondary function:\n";
          }
          newStopPos = positionFuncs.secondaryStopFunc(readData);
        } 
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
        if (verbose) {
          Rcpp::Rcout << "Updated positions for variable class_id: " << element.class_id
                      << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
        }
      } else {
        continue;
      }
    }
    for (auto& element : readData.sigstring) {
      if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class == "read") {
        if(verbose) {
          Rcpp::Rcout << "Dealing with everything that's IS a read!\n";
        } //repetitive code--please fix this if you know what's good for you, you goddamn idiot
        const auto& positionFuncs = positionFuncMap.at(element.class_id);
        int newStartPos = positionFuncs.primaryStartFunc(readData);
        if(newStartPos <= 0){
          newStartPos = positionFuncs.secondaryStartFunc(readData);
        }
        int newStopPos = positionFuncs.primaryStopFunc(readData);
        if(newStopPos <= 0){ 
          newStopPos = positionFuncs.secondaryStopFunc(readData);
        }
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
        if (verbose) {
          Rcpp::Rcout << "Updated positions for the read: " << element.class_id
                      << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
        }
      }
    }
    if(verbose){
      displayOrderedSigString(readData);
    }
  } else if(readType == "FR_RF"){
    // for (auto& element : readData.sigstring) {
    //   if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class != "read") {
    //     if(verbose) {
    //       Rcpp::Rcout << "Dealing with everything that's not a read!\n";
    //       Rcpp::Rcout << "Working on " << element.class_id << "\n";
    //     }
    //     const auto& positionFuncs = positionFuncMap.at(element.class_id);
    //     int newStartPos = positionFuncs.primaryStartFunc(readData);
    //     if(verbose){
    //       Rcpp::Rcout << newStartPos << "\n";
    //     }
    //     if(newStartPos <= 0){
    //       if(verbose){
    //         Rcpp::Rcout << "Triggered secondary func\n";
    //       }
    //       newStartPos = positionFuncs.secondaryStartFunc(readData);
    //     }
    //     int newStopPos = positionFuncs.primaryStopFunc(readData);
    //     if(verbose){
    //       Rcpp::Rcout << newStopPos << "\n";
    //     }
    //     if(newStopPos <= 0){ 
    //       if(verbose){
    //         Rcpp::Rcout << "Triggered secondary func\n";
    //       }
    //       newStopPos = positionFuncs.secondaryStopFunc(readData);
    //     } 
    //     sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
    //       elem.position.first = newStartPos;
    //       elem.position.second = newStopPos;
    //     });
    //     if (verbose) {
    //       Rcpp::Rcout << "Updated positions for variable class_id: " << element.class_id
    //                   << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
    //     }
    //   } else {
    //     continue;
    //   }
    // }
    // for (auto& element : readData.sigstring) {
    //   if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class == "read") {
    //     if(verbose) {
    //       Rcpp::Rcout << "Dealing with everything that's IS a read!\n";
    //     }
    //     const auto& positionFuncs = positionFuncMap.at(element.class_id);
    //     int newStartPos = positionFuncs.primaryStartFunc(readData);
    //     if(newStartPos <= 0){
    //       newStartPos = positionFuncs.secondaryStartFunc(readData);
    //     }
    //     int newStopPos = positionFuncs.primaryStopFunc(readData);
    //     if(newStopPos <= 0){ 
    //       newStopPos = positionFuncs.secondaryStopFunc(readData);
    //     }
    //     sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
    //       elem.position.first = newStartPos;
    //       elem.position.second = newStopPos;
    //     });
    //     if (verbose) {
    //       Rcpp::Rcout << "Updated positions for the read: " << element.class_id
    //                   << ", New Start: " << newStartPos << ", New Stop: " << newStopPos << "\n";
    //     }
    //   }
    // }
    // concat_solve_hybrid(readData, verbose);
    concat_solve_hybrid_v2(readData, positionFuncMap, verbose);
    // for (auto it = readData.sigstring.begin(); it != readData.sigstring.end(); ) {
    //   // Check if element_pass has a value and if the value is true
    //   if((it->element_pass && *(it->element_pass)) || it->type == "variable" || it->global_class == "poly_tail") {
    //     ++it; // Keep this element, move to the next
    //   }
    //   else {
    //     it = readData.sigstring.erase(it); // Erase this element
    //   }
    // }
    if(verbose){
      displayOrderedSigString(readData);
    }
  } else if(readType == "F_F"||readType == "R_R"){
    std::vector<std::string> new_sigs = concat_solve_parallel(readData, verbose);
    //for(const auto& sig: new_sigs){
      //Rcpp::Rcout << sig << "\n";
      //ReadData readData_new;
      //fillSigString(readData_new, readLayout, sig, positionFuncMap, verbose);
      //concat_scan(readData_new, verbose);
      //pos_scan(readData_new, positionFuncMap, readLayout, verbose);
      //break;
    //}
  }
}

Rcpp::CharacterVector SigRun(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigstrings, const PositionFuncMap& positionFuncMap, bool verbose) {
  std::vector<std::string> processed_sigstrings;
  for (const auto& sigstring : sigstrings) {
    ReadData readData;
    fillSigString(readData, readLayout, sigstring, positionFuncMap, verbose);
    concat_scan(readData, verbose);
    pos_scan(readData, positionFuncMap, readLayout, verbose);
    std::string processed_sigstring = generateSigString(readData); // Ensure this function is implemented
    processed_sigstrings.push_back(processed_sigstring);
  }
  return wrap(processed_sigstrings); // Convert std::vector<std::string> to Rcpp::CharacterVector
}

// [[Rcpp::export]]
Rcpp::CharacterVector sig_run(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, 
  const Rcpp::StringVector& rcpp_sigstrings, bool verbose = false) {
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold);
  VarScan(container, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(container, verbose);
  std::vector<std::string> sigstrings = Rcpp::as<std::vector<std::string>>(rcpp_sigstrings);
  return SigRun(container, sigstrings, positionFuncMap, verbose); // Capture and return the result of SigRun
}

