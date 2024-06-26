#include "rad.h"
using namespace std;
using namespace Rcpp;
//preparing the read layout container

ReadLayout container;
PositionFuncMap positionFuncMap;

ReadLayout prep_read_layout_cpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose) {
  ReadLayout container;
  std::unordered_map<std::string, std::tuple<int, int, int>> thresholdsMap;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  // Populate the thresholds map
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
  
  // Iterate over the read_layout and populate the ReadLayout container
  for (int i = 0; i < class_id.size(); ++i) {
    std::string cid = Rcpp::as<std::string>(class_id[i]);
    std::string class_value = Rcpp::as<std::string>(class_column[i]);
    
    boost::optional<std::tuple<int, int, int>> misalignment_threshold_opt;
    if (thresholdsMap.count(cid)) {
      misalignment_threshold_opt = boost::optional<std::tuple<int, int, int>>(thresholdsMap[cid]);
    } else {
      misalignment_threshold_opt = boost::none;
    }
    
    boost::optional<int> expected_length_opt;
    if (class_value != "read") {
      expected_length_opt = expected_length[i];
    }
    
    // Create a new ReadElement and insert it into the container
    ReadElement elem(
        cid,
        Rcpp::as<std::string>(seq[i]),
        expected_length_opt,
        Rcpp::as<std::string>(type[i]),
        order[i],
        Rcpp::as<std::string>(direction[i]),
        class_value,
        misalignment_threshold_opt
    );
    
    container.insert(elem);
    
    if(verbose){
      Rcpp::Rcout << "Inserted ReadElement:\n"
                  << "  class_id: " << elem.class_id << "\n"
                  << "  seq: " << elem.seq << "\n"
                  << "  expected_length: " << (elem.expected_length ? std::to_string(*elem.expected_length) : "none") << "\n"
                  << "  type: " << elem.type << "\n"
                  << "  order: " << elem.order << "\n"
                  << "  direction: " << elem.direction << "\n"
                  << "  global_class: " << elem.global_class << "\n"
                  << "  misalignment_threshold: "
                  << (elem.misalignment_threshold ? 
      std::to_string(std::get<0>(*elem.misalignment_threshold)) + ", " +
      std::to_string(std::get<1>(*elem.misalignment_threshold)) + ", " +
      std::to_string(std::get<2>(*elem.misalignment_threshold)) 
        : "none") << "\n";
    }
  }
  
  return container;
}

// [[Rcpp::export]]
void prep_read_layout_rcpp(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, bool verbose) {
  prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
}

void parse_static_left_new(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::max_element(segment_orders.begin(), segment_orders.end());
  
  for (auto it_var = var_elements.first; it_var != var_elements.second; ++it_var) {
    if (it_var->type != "variable" || it_var->direction != direction) {
      if (verbose) {
        Rcpp::Rcout << "Skipping " << it_var->class_id << "\n";
      }
      continue;
    }
    
    if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    
    if (it_var->order <= read_order && it_var->order > 1) {
      auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
      auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
      auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
      auto next_nextElem = readLayout.get<order_tag>().find(it_var->order + 2);
      
      bool prevElemExists = (prevElem != order_index.end());
      bool nextElemExists = (nextElem != order_index.end());
      bool prePrevElemExists = (pre_prevElem != order_index.end());
      bool nextNextElemExists = (next_nextElem != order_index.end());
      
      int nextElem_val = (nextElemExists && (nextElem->global_class == "stop")) ? 0 : 1;
      int prevElem_val = (prevElemExists && (prevElem->global_class == "start")) ? 0 : 1;
      int prevElem_stop_val = (prevElemExists && (prevElem->global_class == "start")) ? 1 : 0;
      
      
      if (verbose) {
        Rcpp::Rcout << "Processing varElem, global class: " << it_var->global_class << "\n";
        Rcpp::Rcout << "Processing varElem, class ID: " << it_var->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element: " << prevElem->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element, global_class: " << prevElem->global_class << "\n";
        if (prePrevElemExists) Rcpp::Rcout << "Pre-previous element: " << pre_prevElem->class_id << "\n";
        if (nextElemExists) Rcpp::Rcout << "Next element: " << nextElem->class_id << "\n";
        if (nextNextElemExists) Rcpp::Rcout << "Next-next element: " << next_nextElem->class_id << "\n";
      }
      
      int offset_result = 0;
      bool offset_available = true;
      bool is_read = false;
      if (it_var->expected_length && it_var->global_class != "read") {
        offset_result = (*it_var->expected_length);
      } else {
        offset_result = 0;
        offset_available = false;
        is_read = true;
      }
      if (verbose) {
        Rcpp::Rcout << "Offset result is " << std::to_string(offset_result) << "\n";
        Rcpp::Rcout << "Offset available: " << (offset_available ? "true" : "false") << "\n";
      }
      std::string primaryStart, primaryStop, secondaryStart, secondaryStop;
      if (prevElemExists) {
        primaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
        if (is_read) {
          if (verbose) {
            Rcpp::Rcout << "Found the read.\n";
            if (prevElemExists) Rcpp::Rcout << "prevElem exists: true\nprevElem class_id: " << prevElem->class_id << "\n";
            if (nextElemExists) Rcpp::Rcout << "nextElem exists: true\nnextElem class_id: " << nextElem->class_id << "\n";
          }
          bool prevPolyTail = prevElemExists && prevElem->global_class == "poly_tail";
          bool nextPolyTail = nextElemExists && nextElem->global_class == "poly_tail";
          bool prevShortStatic = prevElemExists && prevElem->type == "static" && *prevElem->expected_length <= 13 && prevElem->global_class != "poly_tail";
          bool nextShortStatic = nextElemExists && nextElem->type == "static" && *nextElem->expected_length <= 13 && nextElem->global_class != "poly_tail";
          if (prevPolyTail || prevShortStatic) {
            if (prePrevElemExists && pre_prevElem->type == "static") {
              if (prevPolyTail) {
                secondaryStart = pre_prevElem->class_id + "|stop+1|poly_skipped";
              } else {
                int prevElem_offset = *prevElem->expected_length;
                secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset) + "|skipped";
              }
            }
          }
          if (nextPolyTail || nextShortStatic) {
            if (nextNextElemExists && next_nextElem->type == "static") {
              secondaryStop = next_nextElem->class_id + "|start-1";
              if (nextPolyTail) {
                secondaryStop += "|poly_skipped";
              }
            }
          }
        } else {
          primaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result - prevElem_stop_val);
        }
        // Check for chained variables
        if (prevElemExists && prevElem->type == "variable" &! is_read) {
          auto& prevAssigned = assignedPositions[prevElem->class_id];
          primaryStart += "|chained_start";
          primaryStop += "|chained_stop";
          if (nextElemExists && nextElem->type == "static") {
            secondaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result);
            secondaryStop = nextElem->class_id + "|start-1";
          } else {
            secondaryStart = std::get<0>(prevAssigned) + "|chained_start";
            secondaryStop = std::get<2>(prevAssigned) + "|chained_stop";
          }
        } else if (prevElemExists && prevElem->type == "static" &! is_read) {
          if (nextElemExists && nextElem->type == "static" && (nextElem->expected_length >= 13) && nextElem->global_class != "poly_tail") {
            secondaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result);
            secondaryStop = nextElem->class_id + "|start-1";
          } else {
            if (nextNextElemExists && next_nextElem->type == "static" && nextElem->global_class != "read" && next_nextElem->global_class != "poly_tail") {
              int nextElem_offset = *nextElem->expected_length;
              secondaryStart = next_nextElem->class_id + "|start-" + std::to_string(nextElem_offset+offset_result) + "|skipped";
              secondaryStop = next_nextElem->class_id + "|start-" + std::to_string(nextElem_offset) + "|skipped";
            } else {
              if (prevElem->global_class == "poly_tail" && prePrevElemExists && pre_prevElem->global_class != "read") {
                secondaryStart = pre_prevElem->class_id + "|stop+1|poly_skipped";
              } else {
                secondaryStart = prevElem->class_id + "|stop+1|left_terminal_linked";
                secondaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result) + "|left_terminal_linked";
              }
            }
          }
        }
        // Update assignedPositions
        if (assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
          auto& existingEntry = assignedPositions[it_var->class_id];
          if (std::get<0>(existingEntry) != primaryStart) {
            std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
          }
          if (std::get<2>(existingEntry) != primaryStop) {
            std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
          }
        } else {
          assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
        }
        if (verbose) {
          Rcpp::Rcout << "Parse Left: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
          Rcpp::Rcout << "Parse Left: Primary Stop " << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
        }
      }
    }
    if (verbose) {
      Rcpp::Rcout << "Completed parsing left varElem: " << it_var->class_id << "\n";
    }
  }
}

void parse_static_right_new(const ReadLayout& readLayout, const std::vector<int>& segment_orders,
  const std::string& direction, bool verbose, std::unordered_map<std::string,
    std::tuple<std::string, std::string, std::string, std::string>>& assignedPositions) {
  auto& order_index = readLayout.get<order_tag>();
  auto& type_index = readLayout.get<type_tag>();
  auto var_elements = type_index.equal_range("variable");
  int read_order = *std::min_element(segment_orders.begin(), segment_orders.end());
  for (auto it_var = std::make_reverse_iterator(var_elements.second); it_var != std::make_reverse_iterator(var_elements.first); ++it_var) {
    if (it_var->type != "variable" || it_var->direction != direction) {
      if (verbose) {
        Rcpp::Rcout << "Skipping " << it_var->class_id << "\n";
      }
      continue;
    }
    if (std::find(segment_orders.begin(), segment_orders.end(), it_var->order) == segment_orders.end()) {
      continue;
    }
    
    if (it_var->order >= read_order) {
      auto nextElem = readLayout.get<order_tag>().find(it_var->order + 1);
      auto prevElem = readLayout.get<order_tag>().find(it_var->order - 1);
      auto next_nextElem = readLayout.get<order_tag>().find(it_var->order + 2);
      auto pre_prevElem = readLayout.get<order_tag>().find(it_var->order - 2);
      
      bool nextElemExists = (nextElem != order_index.end());
      bool prevElemExists = (prevElem != order_index.end());
      bool nextNextElemExists = (next_nextElem != order_index.end());
      bool prePrevElemExists = (pre_prevElem != order_index.end());
      
      int nextElem_val = (nextElemExists && (nextElem->global_class == "stop")) ? 0 : 1;
      int prevElem_val = (prevElemExists && (prevElem->global_class == "start")) ? 0 : 1;
      int nextElem_stop_val = (nextElemExists && (nextElem->global_class == "stop")) ? 1 : 0;
      
      if (verbose) {
        Rcpp::Rcout << "Processing varElem, global class: " << it_var->global_class << "\n";
        Rcpp::Rcout << "Processing varElem, class ID: " << it_var->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element: " << prevElem->class_id << "\n";
        if (prevElemExists) Rcpp::Rcout << "Previous element, global class: " << prevElem->global_class << "\n";
        if (prePrevElemExists) Rcpp::Rcout << "Pre-previous element: " << pre_prevElem->class_id << "\n";
        if (nextElemExists) Rcpp::Rcout << "Next element: " << nextElem->class_id << "\n";
        if (nextNextElemExists) Rcpp::Rcout << "Next-next element: " << next_nextElem->class_id << "\n";
      }
      int offset_result = 0;
      bool offset_available = true;
      bool is_read = false;
      if (it_var->expected_length && it_var->global_class != "read") {
        offset_result = (*it_var->expected_length);
      } else {
        offset_available = false;
        is_read = true;
      }
      if (verbose) {
        Rcpp::Rcout << "Offset result is " << std::to_string(offset_result) << "\n";
        Rcpp::Rcout << "Offset available: " << (offset_available ? "true" : "false") << "\n";
      }
      std::string primaryStart, primaryStop, secondaryStart, secondaryStop;
      if (prevElemExists && nextElemExists) {
        primaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
        primaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val+nextElem_stop_val);
        secondaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val) + "|left_terminal_linked";
        secondaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val) + "|right_terminal_linked";
        // Handle read element separately
        if (is_read) {
          if (verbose) {
            Rcpp::Rcout << "Found the read.\n";
            if (prevElemExists) Rcpp::Rcout << "prevElem exists: true\nprevElem class_id: " << prevElem->class_id << "\n";
            if (nextElemExists) Rcpp::Rcout << "nextElem exists: true\nnextElem class_id: " << nextElem->class_id << "\n";
          }
          // Check for poly_tail or short static element in previous or next elements
          bool prevPolyTail = prevElemExists && prevElem->global_class == "poly_tail";
          bool nextPolyTail = nextElemExists && nextElem->global_class == "poly_tail";
          bool prevShortStatic = prevElemExists && prevElem->type == "static" && *prevElem->expected_length <= 13 && prevElem->global_class != "poly_tail";
          bool nextShortStatic = nextElemExists && nextElem->type == "static" && *nextElem->expected_length <= 13 && nextElem->global_class != "poly_tail";
          if(verbose){
            if(prevPolyTail) Rcpp::Rcout << "Previous element is a poly tail!\n";
            if(nextPolyTail) Rcpp::Rcout << "Next element is a poly tail!\n";
          }
          if (prevPolyTail || prevShortStatic) {
            if (prePrevElemExists) {
              if (prevPolyTail) {
                secondaryStart = pre_prevElem->class_id + "|stop+1|poly_skipped";
              } else {
                int prevElem_offset = *prevElem->expected_length;
                secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset) + "|skipped";
              }
            }
          }
          if (nextPolyTail || nextShortStatic) {
            if (nextNextElemExists) {
              if (nextPolyTail) {
                secondaryStop = next_nextElem->class_id + "|stop+1|poly_skipped";
              } else {
                int nextElem_offset = *nextElem->expected_length;
                secondaryStop = next_nextElem->class_id + "|stop+" + std::to_string(nextElem_offset) + "|skipped";
              }
            }
          }
          else {
            secondaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val) + "|right_terminal_linked";
          }
        } else {
          primaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result);
        }
        // Check for chained variables
        if (nextElemExists && nextElem->type == "variable" && !is_read) {
          auto& nextAssigned = assignedPositions[nextElem->class_id];
          primaryStart += "|chained_start";
          primaryStop += "|chained_stop";
          if (prevElemExists && prevElem->type == "static") {
            secondaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
            secondaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result);
          } else {
            secondaryStart = std::get<0>(nextAssigned) + "|chained_start";
            secondaryStop = std::get<2>(nextAssigned) + "|chained_stop";
          }
        } else if (nextElemExists && nextElem->type == "static" && !is_read) {
          if (prevElemExists && prevElem->type == "static") {
            secondaryStart = prevElem->class_id + "|stop+" + std::to_string(prevElem_val);
            secondaryStop = prevElem->class_id + "|stop+" + std::to_string(offset_result);
          } else {
            if (prePrevElemExists && pre_prevElem->type == "static" && prevElem->global_class != "read") {
              int prevElem_offset = *prevElem->expected_length;
              secondaryStart = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset);
              secondaryStop = pre_prevElem->class_id + "|stop+" + std::to_string(prevElem_offset + offset_result);
            } else if (prePrevElemExists && pre_prevElem->type == "variable" && prevElem->global_class != "read") {
              if (prevElem->global_class == "poly_tail"||*prevElem->expected_length < 13) {
                secondaryStart = pre_prevElem->class_id + "|stop+1|skipped";
                secondaryStop = pre_prevElem->class_id + "|stop+" + std::to_string(offset_result) + "|skipped";
              } else {
                secondaryStart = nextElem->class_id + "|start-" + std::to_string(offset_result) + "|right_terminal_linked";
                secondaryStop = nextElem->class_id + "|start-" + std::to_string(nextElem_val) + "|right_terminal_linked";
              }
            }
          }
        }
        // Update assignedPositions
        if (assignedPositions.find(it_var->class_id) != assignedPositions.end()) {
          auto& existingEntry = assignedPositions[it_var->class_id];
          if (std::get<0>(existingEntry) != primaryStart) {
            std::get<1>(existingEntry) = primaryStart; // Update secondary start with new primary start
          }
          if (std::get<2>(existingEntry) != primaryStop) {
            std::get<3>(existingEntry) = primaryStop;  // Update secondary stop with new primary stop
          }
        } else {
          assignedPositions[it_var->class_id] = std::make_tuple(primaryStart, secondaryStart, primaryStop, secondaryStop);
        }
        if (verbose) {
          Rcpp::Rcout << "Parse Right: Primary Start " << primaryStart << " & Secondary Start: " << secondaryStart << "\n";
          Rcpp::Rcout << "Parse Right: Primary Stop " << primaryStop << " & Secondary Stop: " << secondaryStop << "\n";
        }
      }
    }
    if (verbose) {
      Rcpp::Rcout << "Completed parsing right varElem: " << it_var->class_id << "\n";
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
    
    if (verbose) {
      Rcpp::Rcout << "Read direction: " << read_direction << ", Read order: " << read_order << "\n";
      Rcpp::Rcout << "Left vector: ";
      for (const auto& order : left_vector) {
        Rcpp::Rcout << order << " ";
      }
      Rcpp::Rcout << "\nRight vector: ";
      for (const auto& order : right_vector) {
        Rcpp::Rcout << order << " ";
      }
      Rcpp::Rcout << "\n";
    }
    
    //Parse static elements
    parse_static_left_new(readLayout, left_vector, read_direction, verbose, assignedPositions); //, true);
    parse_static_right_new(readLayout, right_vector, read_direction, verbose, assignedPositions);//, true);
    //Parse read elements
    //parse_static_left_new(readLayout, left_vector, read_direction, verbose, assignedPositions); //, false);
   //parse_static_right_new(readLayout, right_vector, read_direction, verbose, assignedPositions);//, false);
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

// [[Rcpp::export]]
void single_test_run_new(Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold) {
  container = prep_read_layout_cpp(read_layout, misalignment_threshold, true);
  VarScan(container, true);
  using PositionFuncMap = std::unordered_map<std::string, PositionCalcFunc>;
  positionFuncMap = createPositionFunctionMap(container, true);
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

// Main function to filter sigstrings--TODO ADD FLAGS FOR FILTERING
std::vector<std::string> split(const std::string &s, char delimiter) {
std::vector<std::string> tokens;
std::string token;
std::istringstream tokenStream(s);
while (std::getline(tokenStream, token, delimiter)) {
  tokens.push_back(token);
}
return tokens;
}
void filter_sigstring(std::vector<std::string>& sigstrings) {
  std::regex readPattern(R"((read|rc_read):(\d+):(\d+):(\d+))");
  std::regex barcodePattern(R"((barcode(?:_\d+)?|rc_barcode(?:_\d+)?):(\d+):(\d+):(\d+))");
  std::smatch matches;
  for (auto& sigstring : sigstrings) {
    bool hasRead = false, hasBarcode = false, isValid = true;
    auto elements = split(sigstring, '|');
    for (const auto& element : elements) {
      if (std::regex_search(element, matches, readPattern)) {
        
        int startPos = std::stoi(matches.str(3));
        int stopPos = std::stoi(matches.str(4));
        //Rcpp::Rcout << "Caught at this filter because I think that " << stopPos << " is this and " << startPos << " is this!\n";
        if ((stopPos - startPos) <= 0 || (stopPos - startPos) < 100) {
          isValid = false;
          break;
        }
        hasRead = true;
      } else if (std::regex_search(element, matches, barcodePattern)) {
        int startPos = std::stoi(matches.str(3));
        int stopPos = std::stoi(matches.str(4));
        // Check barcode length
        if ((stopPos - startPos + 1) < 4) {
          //Rcpp::Rcout << "Caught at this filter because I think that " << stopPos << " is this and " << startPos << " is this!\n";
          
          isValid = false;
         break;
        }
        hasBarcode = true;
      }
    }
    if (!isValid || !(hasRead && hasBarcode)) {
      size_t endMetaPos = sigstring.rfind('>');
      if (endMetaPos != std::string::npos) {
        // Find the last colon before the '>' to identify where the 'F' or 'R' starts
        size_t lastColonPos = sigstring.rfind(':', endMetaPos);
        if (lastColonPos != std::string::npos) {
          // Replace the 'F' or 'R' with "undecided"
          std::string beforeType = sigstring.substr(0, lastColonPos + 1); // Keep everything before the type
          std::string afterType = sigstring.substr(endMetaPos); // Keep everything after the type
          sigstring = beforeType + "undecided" + afterType; // Combine with "undecided"
        }
      }
    }
  }
}

//using PositionFuncMap = std::unordered_map<std::string, PositionCalcFunc>;
//using positionFuncMap = std::unordered_map<std::string, PositionCalcFunc>;
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
      return [varElemClassId, refClassId, offset, additionalInfo, isPrimary, 
        operation, verbose](const ReadData& readData) -> int {
                    const auto& id_index = readData.sigstring.get<sig_id_tag>();
                    const auto& order_index = readData.sigstring.get<sig_order_tag>();
                    auto refClass = id_index.find(refClassId);
                    auto it = id_index.equal_range(refClassId).first;
                    if (verbose) {
                      Rcpp::Rcout << "Made it into the lambda function!\n";
                      Rcpp::Rcout << refClass->type << " is the type of the element.\n";
                      std::string result = isPrimary ? "Primary" : "Secondary";
                      Rcpp::Rcout << "The position is flagged as a " << result << "\n";
                    }
                    if (it != id_index.end() && 
                        ((refClass->type == "static" && refClass->element_pass && *refClass->element_pass) || 
                        (refClass->type == "variable") || (refClass->global_class == "poly_tail")) ||
                        (refClass->global_class == "start") || (refClass->global_class == "stop")) {
                      
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
                        if(startPos < 0 && stopPos < 0){
                          if(verbose){
                            Rcpp::Rcout << "Returning because it's a primary run dependent on an unannotated element!\n";
                          }
                          return -1;
                        } else {
                          return calculatedPos;
                        }
                    } else if (!isPrimary && !additionalInfo.empty()) {
                      std::string readType = readData.string_info.at("type").str_type;
                      auto varClass = id_index.find(varElemClassId);
                      if(verbose){
                        std::string result = isPrimary ? "Primary" : "Secondary";
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
                          rel_edit_failure = static_cast<int>((1.0 - (edit_distance / sequence_length)) * 100);                                         
                          if(rel_edit_failure >= 70 && verbose){
                            int new_terminal_pos;
                            if(operation == "start"){
                              new_terminal_pos = prevClass->position.second +1;
                            } else {
                              new_terminal_pos = prevClass->position.second + offset;
                            }
                          if(readData.string_info.at("type").str_type == "F" || readData.string_info.at("type").str_type == "R"){
                              Rcpp::Rcout << operation << " is the operation!\n";
                              Rcpp::Rcout << offset << " is the offset!\n";
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
                    return -1; // Default return value if conditions are not met
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

std::vector<std::string> processSigString(const std::string& sigstring) {
  std::vector<std::string> parts;
  size_t splitPos = sigstring.find("+");
  
  // Extracting the "<length:id:type>" part correctly
  size_t infoStartPos = sigstring.rfind("<") + 1;
  size_t infoEndPos = sigstring.find(">", infoStartPos);
  std::string infoPart = sigstring.substr(infoStartPos, infoEndPos - infoStartPos);
  // Parsing the length, ID, and type from the infoPart
  std::istringstream infoStream(infoPart);
  std::string lengthStr, id, type;
  std::getline(infoStream, lengthStr, ':');
  std::getline(infoStream, id, ':');
  std::getline(infoStream, type);
  if (splitPos != std::string::npos) {
    std::string firstPart = sigstring.substr(0, splitPos);
    std::string secondPart = sigstring.substr(splitPos + 1);
    // Finding the last ":" before "<" in each part to capture the correct length
    size_t lastColonInFirst = firstPart.rfind(":", firstPart.rfind("<"));
    std::string lengthFirst = firstPart.substr(lastColonInFirst + 1, firstPart.rfind("<") - lastColonInFirst - 1);
    size_t lastColonInSecond = secondPart.rfind(":", secondPart.rfind("<"));
    std::string lengthSecond = secondPart.substr(lastColonInSecond + 1, secondPart.rfind("<") - lastColonInSecond - 1);
    // Reconstructing the first and second parts with updated length and type
    size_t endFirstPart = firstPart.rfind("<");
    firstPart = firstPart.substr(0, endFirstPart) + "<" + lengthFirst + ":" + id + "+FR_RF:F>";
    size_t endSecondPart = secondPart.rfind("<");
    secondPart = secondPart.substr(0, endSecondPart) + "<" + lengthSecond + ":" + id + "+FR_RF:R>";
    parts.push_back(firstPart);
    parts.push_back(secondPart);
  } else {
    parts.push_back(sigstring); // Handling cases without "+"
  }
  return parts;
}

void concat_solve_hybrid(ReadData& readData, const PositionFuncMap& positionFuncMap, bool verbose) {
  std::string unrelated = "forw_primer:0:0:0|rev_primer:0:0:0+rc_rev_primer:0:0:0|rc_forw_primer:0:0:0<10:@theID:FR_RF>";
  std::vector<std::string> results = processSigString(unrelated);
  
  if(verbose){
    Rcpp::Rcout << "Solving a hybrid concatenate! First checking to see if the deconcatenating function works:\n";
    for (const auto& part : results) {
      Rcpp::Rcout << part << "\n";
   }
  }
  auto& sig_id_index = readData.sigstring.get<sig_id_tag>();
  for (auto& element : readData.sigstring) {
    if (element.type == "variable" && positionFuncMap.count(element.class_id) > 0 && element.global_class != "read") {
      if(verbose) {
        Rcpp::Rcout << "Working on resolving: " << element.class_id << "\n";
        Rcpp::Rcout << "Now making the first pass...\n";
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
      if(newStartPos == -1 || newStopPos == -1){
        auto prevClass = readData.sigstring.get<sig_order_tag>().find(element.order - 1);
        if(prevClass->global_class == "poly_tail"){
          if(verbose){
            Rcpp::Rcout << "The previous class in the order was a poly-tail. Skipping back one even farther...\n";
          }
           prevClass = readData.sigstring.get<sig_order_tag>().find(element.order -2);
        }
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
        }
      } else
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
    }
    else {
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
      if (newStartPos == -1 || newStopPos == -1){
        auto prevClass = readData.sigstring.get<sig_order_tag>().find(element.order - 1);
        if(prevClass->global_class == "poly_tail"||(prevClass->global_class == "other" && (prevClass->position.second- prevClass->position.first < 15))){
          if(verbose){
            Rcpp::Rcout << "The previous class in the order was a poly-tail, or too short. Skipping back one even farther...\n";
          }
          prevClass = readData.sigstring.get<sig_order_tag>().find(element.order -2);
        }
        
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
        if(newStopPos <= 0){
          if(verbose){
            Rcpp::Rcout << "Something's weird here. New Stop Pos is " << newStopPos << "\n";
          }
          if(newStopPos == -1){
            if(verbose){
              Rcpp::Rcout << "New Stop Pos is -1, meaning the right bound of the object still hasn't been set.\n";
            }
            auto nextClass = readData.sigstring.get<sig_order_tag>().find(element.order + 2);
            if(nextClass->global_class == "poly_tail"||nextClass->global_class == "rev_primer"||(nextClass->global_class == "other" && 
              (nextClass->position.second - nextClass->position.first < 15))){
              if(verbose){
                Rcpp::Rcout << "The next class in the order was extraneous or short. Skipping one even farther...\n";
              }
              nextClass = readData.sigstring.get<sig_order_tag>().find(element.order  + 3);
          }
            if(verbose){
              Rcpp::Rcout << "Last-ditch effort to save this one...\n";
              Rcpp::Rcout << "Your final hope rests on " << nextClass->class_id << "\n";
              Rcpp::Rcout << "Your final hope rests on the position" << nextClass->position.first << "\n";
            }
          newStopPos = nextClass->position.first-1;
        }
        sig_id_index.modify(sig_id_index.iterator_to(element), [newStartPos, newStopPos](SigElement& elem) {
          elem.position.first = newStartPos;
          elem.position.second = newStopPos;
        });
        }
        
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
        }
      } else
      if (newStartPos == -2 || newStopPos == -2){
            auto nextClass = readData.sigstring.get<sig_order_tag>().find(element.order + 1);
        if(nextClass->global_class == "poly_tail"||(nextClass->global_class == "other" && (nextClass->position.second - nextClass->position.first < 15))){
          if(verbose){
            Rcpp::Rcout << "The next class in the order was extraneous or short. Skipping one even farther...\n";
          }
          nextClass = readData.sigstring.get<sig_order_tag>().find(element.order  +2);
        }
            if(verbose) {
              Rcpp::Rcout << "The next class in order is " << nextClass->class_id << "\n";
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
  std::vector<std::string> sig_single_pass = generateSigString(readData);
  std::vector<std::string> deconcatenated = processSigString(sig_single_pass[0]);
  //std::string initial_pass = sig_single_pass.at(0);
  if(verbose){
    Rcpp::Rcout << "This is what a concatenate OG solution looks like:\n";
    Rcpp::Rcout << sig_single_pass[0] << "\n";
    Rcpp::Rcout << "And this is the deconcatenated version:\n";
    for (const auto& part : deconcatenated) {
      Rcpp::Rcout << part << "\n";
    }
  }
  std::string forward_deconcatenated = deconcatenated[0];
  std::string reverse_deconcatenated = deconcatenated[1];
  // Assuming 'readData' is an instance of ReadData already defined elsewhere
  if(!forward_deconcatenated.empty() || !reverse_deconcatenated.empty()) {
    readData.sigstring.clear();
    readData.secondary_sigstring = boost::none; // Reset optional secondary sigstring
  }
  // Process first sigstring
  if (!forward_deconcatenated.empty()) {
    fillSigString(readData, container, forward_deconcatenated, positionFuncMap, verbose);
  }
  if (!forward_deconcatenated.empty()) {
    if (!reverse_deconcatenated.empty()) {
      fillSigString(readData, container, reverse_deconcatenated, positionFuncMap, verbose);
    }
  }
  readData.string_info["type"].str_type = "F-R";
  readData.string_info["type"].split_types = std::make_pair("F","R");
}

void concat_solve_parallel(ReadData& readData, bool verbose) {
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
  std::string direction;
  if(readData.string_info.at("type").str_type == "F_F"){
    direction = "F";
  } else {
    direction = "R";
  }
  std::string id = readData.string_info.at("type").str_id;
  int splitPoint = decreasePoints.empty() ? elements.size() : decreasePoints[0] + 1;
  // Build the first sigstring
  for (int i = 0; i < splitPoint; ++i) {
    if (i > 0) sigString1 += "|";
    sigString1 += elements[i].class_id + ":" + std::to_string(elements[i].edit_distance.value_or(0)) + ":" + std::to_string(elements[i].position.first) + ":" + std::to_string(elements[i].position.second);
  }
  // Append the new length position for the first string
  if (!elements.empty() && splitPoint > 0) {
   sigString1 += "<" + std::to_string(elements[splitPoint - 1].position.second) + ":" + id +  ":" + direction + ">";
  }
   //Build the second sigstring
  for (size_t i = splitPoint; i < elements.size(); ++i) {
    if (i > splitPoint) sigString2 += "|";
    sigString2 += elements[i].class_id + ":" + std::to_string(elements[i].edit_distance.value_or(0)) + ":" + std::to_string(elements[i].position.first) + ":" + std::to_string(elements[i].position.second);
  }
  // Append the new length position for the second string
  if (!elements.empty() && splitPoint < elements.size()) {
   sigString2 += "<"+ std::to_string(elements.back().position.second) + ":" + id + ":" + direction + ">"; 
  }
  // Verbose output to demonstrate the two new sigstrings
  if (verbose) {
    Rcpp::Rcout << "First new sigstring after correction:\n" << sigString1 << "\n";
    Rcpp::Rcout << "Second new sigstring after correction:\n" << sigString2 << "\n";
  }
  // Assuming 'readData' is an instance of ReadData already defined elsewhere
  if(!sigString1.empty() || !sigString2.empty()) {
    // If either sigString is not empty, clear existing data to prepare for new data
    readData.sigstring.clear();
    readData.secondary_sigstring = boost::none; // Reset optional secondary sigstring
  }
  // Process first sigstring
  if (!sigString1.empty()) {
    fillSigString(readData, container, sigString1, positionFuncMap, verbose);
    pos_scan(readData, positionFuncMap, container, verbose);
    //sigString1 = generateSigString(readData); // Reassign after processing
  }
  // Process second sigstring, if present
  if (!sigString2.empty()) {
    // Check if the primary sigstring was populated, then use secondary sigstring
    if (!readData.sigstring.empty()) {
      // Initialize secondary sigstring and process it
      fillSigString(readData, container, sigString2, positionFuncMap, verbose); // Note: Assuming fillSigString can accept a flag to fill secondary sigstring
      pos_scan(readData, positionFuncMap, container, verbose); // Assuming pos_scan can also handle secondary sigstring
      //sigString2 = generateSigString(readData); // Generate from secondary, assuming overload of generateSigString for secondary sigstring
    } else {
      // If primary sigstring was not used, fill it now
      fillSigString(readData, container, sigString2, positionFuncMap, verbose);
      pos_scan(readData, positionFuncMap, container, verbose);
      //sigString2 = generateSigString(readData); // Generate from primary as fallback
    }
  }
  if (verbose) {
    Rcpp::Rcout << "First new sigstring: " << sigString1 << "\n";
    Rcpp::Rcout << "Second new sigstring: " << sigString2 << "\n";
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
  
  SigString* targetSigString = &readData.sigstring; // Default to primary sigstring
  if (!readData.sigstring.empty() && !readData.secondary_sigstring) {
    // Condition 2: sigstring is populated and secondary_sigstring is not, populate secondary_sigstring.
    readData.secondary_sigstring = SigString(); // Initialize secondary_sigstring
    targetSigString = &(*readData.secondary_sigstring); // Target secondary sigstring
  } else if (!readData.sigstring.empty() && readData.secondary_sigstring) {
    // Condition 3: Both sigstring and secondary_sigstring are populated, empty sigstring and refill it with the new sigstring data.
    readData.sigstring.clear(); // Empty primary sigstring
    targetSigString = &readData.sigstring; // Target primary sigstring
  }
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
      targetSigString->insert(std::move(element));
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
      targetSigString->insert(std::move(element));
    }
  }
}

using df_map = std::unordered_map<std::string, NanoRead>;
void df_to_unordered_map(const Rcpp::DataFrame& df, df_map& map, bool verbose) {
  Rcpp::StringVector ids = df["id"];
  Rcpp::StringVector quals = df["qual"];
  Rcpp::StringVector seqs = df["seq"];
  
  for(int i = 0; i < ids.size(); ++i) {
    std::string id = Rcpp::as<std::string>(ids[i]);
    std::string qual = Rcpp::as<std::string>(quals[i]);
    std::string seq = Rcpp::as<std::string>(seqs[i]);
    map[id] = {qual, seq};
  }
  if (verbose) {
    Rcpp::Rcout << "Map populated. Entries: " << map.size() << "\n";
  }
}

using df_map_v2 = std::unordered_map<std::string, NanoRead_v2>;
void df_to_unordered_map_refs(const Rcpp::DataFrame& df, df_map_v2& map, bool verbose) {
  Rcpp::StringVector ids = df["id"];
  Rcpp::StringVector quals = df["qual"];
  Rcpp::StringVector seqs = df["seq"];
  for(int i = 0; i < ids.size(); ++i) {
    std::string id = Rcpp::as<std::string>(ids[i]);
    map[id] = {static_cast<std::size_t>(i), static_cast<std::size_t>(i)};
  }
  if (verbose) {
    Rcpp::Rcout << "Map populated. Entries: " << map.size() << "\n";
  }
}

void clearMap(df_map& map, bool verbose) {
  if (verbose) {
    Rcpp::Rcout << "Clearing map. Entries before: " << map.size();
  }
  map.clear();
  if (verbose) {
    Rcpp::Rcout << ". Entries after: " << map.size() << "\n";
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame sig_extraction(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold,
  Rcpp::DataFrame& df, const Rcpp::CharacterVector& processed_sigstrings, bool verbose) {
  // Convert R DataFrame to C++ ReadLayout structure
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  // Process variable scanning and create position function map
  VarScan(readLayout, verbose);
  positionFuncMap = createPositionFunctionMap(readLayout, verbose);
  // Convert processed_sigstrings to std::vector<std::string>
  std::vector<std::string> sigstrings = Rcpp::as<std::vector<std::string>>(processed_sigstrings);
  // Extract all unique class IDs, including both forward and reverse complements
  std::set<std::string> uniqueClassIds;
  auto& class_id_index = readLayout.get<id_tag>();
  for (const auto& elem : class_id_index) {
    if (elem.type == "variable") {
      if (elem.class_id.substr(0, 3) == "rc_") {
        uniqueClassIds.insert(elem.class_id.substr(3));
      } else {
        uniqueClassIds.insert(elem.class_id);
      }
    }
  }
  // Create a map for fast lookup of sequences by ID
  std::unordered_map<std::string, std::string> seq_map;
  Rcpp::StringVector ids = df["id"];
  Rcpp::StringVector seqs = df["seq"];
  for (int i = 0; i < ids.size(); ++i) {
    seq_map[Rcpp::as<std::string>(ids[i])] = Rcpp::as<std::string>(seqs[i]);
  }
  // Prepare the result list
  Rcpp::List result(sigstrings.size());
  std::vector<std::string> column_names(uniqueClassIds.begin(), uniqueClassIds.end());
  for (int i = 0; i < sigstrings.size(); ++i) {
    Rcpp::List read_result;
    ReadData readData;
    std::string signature = sigstrings[i];
    // Parse the signature
    auto read_info_pos = signature.find_last_of('<');
    std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2);
    std::stringstream lengthTypeStream(lengthTypeStr);
    std::string lengthStr, read_id, type;
    std::getline(lengthTypeStream, lengthStr, ':');
    std::getline(lengthTypeStream, read_id, ':');
    std::getline(lengthTypeStream, type);
    // Debugging statements
    if(verbose){
      Rcpp::Rcout << "Signature: " << signature << std::endl;
      Rcpp::Rcout << "Length: " << lengthStr << ", Read ID: " << read_id << ", Type: " << type << std::endl;
    }
    auto plus_pos = read_id.find('+');
    if (plus_pos != std::string::npos) {
      read_id = read_id.substr(0, plus_pos);
    }
    if(type == "undecided"){
      continue;
    }
    int length = std::stoi(lengthStr);
    readData.addStringInfo("type", length, read_id, type);
    SigString* targetSigString = &readData.sigstring;
    std::stringstream ss(signature.substr(0, read_info_pos));
    std::string token;
    
    while (std::getline(ss, token, '|')) {
      std::string id;
      int editDistance, startPos, endPos;
      std::stringstream tokenStream(token);
      
      std::getline(tokenStream, id, ':');
      tokenStream >> editDistance;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> startPos;
      tokenStream.ignore(1); // Ignore the colon
      tokenStream >> endPos;
      
      auto it = class_id_index.find(id);
      if (it != class_id_index.end() && it->type == "variable") {
        auto seq_it = seq_map.find(read_id);
        if (seq_it == seq_map.end()) {
          continue;
        }
        std::string seq = seq_it->second;
        std::string extracted_seq;
        if (startPos > 0 && endPos >= startPos && (seq.length() >= endPos)) {
          extracted_seq = seq.substr(startPos - 1, endPos - startPos + 1); // Corrected substring extraction
          if (it->direction != "forward") {
            extracted_seq = revcomp_cpp(extracted_seq);
          }
        } else {
          continue;
        }
        std::string classId = id;
        if (classId.substr(0, 3) == "rc_") {
          classId = classId.substr(3);
        }
        
        read_result[classId] = extracted_seq;
      }
    }
    result[i] = read_result;
  }
  // Convert the list to a data.table
  Rcpp::List data_table_list;
  data_table_list["id"] = processed_sigstrings;
  for (const auto& classId : uniqueClassIds) {
    Rcpp::CharacterVector column(sigstrings.size(), NA_STRING);
    for (int i = 0; i < sigstrings.size(); ++i) {
      Rcpp::List read_result = result[i];
      if (read_result.containsElementNamed(classId.c_str())) {
        column[i] = Rcpp::as<std::string>(read_result[classId]);
      }
    }
    data_table_list[classId] = column;
  }
  // Create the data.table
  data_table_list.attr("class") = "data.table";
  data_table_list.attr("row.names") = Rcpp::seq(1, sigstrings.size());
  return data_table_list;
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

std::vector<std::string> generateSigString(ReadData& readData) {  // Note: passing readData by non-const reference
  std::vector<std::string> sigStrings; // Vector to hold one or two sigstrings
  
  auto& typeInfo = readData.string_info["type"];  // Getting a modifiable reference to typeInfo
  
  auto generateFromSigString = [&](const SigString& sigstring, const std::string& typeOverride) -> std::string {
    std::stringstream sigStringStream;
    auto appendElements = [&](const std::string& direction) {
      std::vector<SigElement> sorted_elements;
      auto& dir_index = sigstring.get<sig_dir_tag>();
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
    appendElements("forward");
    if (typeOverride == "R" || typeOverride == "FR_RF") {
      if (!sigStringStream.str().empty()) sigStringStream << "+";
      appendElements("reverse");
    }
    sigStringStream << "<" << typeInfo.str_length << ":" << typeInfo.str_id << ":" << typeOverride << ">";
    return sigStringStream.str();
  };
  // Logic to handle F-R type
  if(typeInfo.str_type == "F-R") {
    sigStrings.push_back(generateFromSigString(readData.sigstring, "F"));
    if(readData.secondary_sigstring) {
      sigStrings.push_back(generateFromSigString(*readData.secondary_sigstring, "R"));
    }
  } else {
    // Default behavior for other types
    sigStrings.push_back(generateFromSigString(readData.sigstring, typeInfo.str_type));
    if(readData.secondary_sigstring) {
      sigStrings.push_back(generateFromSigString(*readData.secondary_sigstring, typeInfo.str_type));
    }
  }
  return sigStrings;
}

void pos_scan(ReadData& readData, const PositionFuncMap& positionFuncMap, const ReadLayout& readLayout, 
  bool verbose) {
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
    concat_solve_hybrid(readData, positionFuncMap, verbose);
  } else if(readType == "F_F"||readType == "R_R"){
    concat_solve_parallel(readData, verbose);
  }
}

Rcpp::CharacterVector process_sigstrings(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigs, const PositionFuncMap& positionFuncMap, bool verbose, int nthreads = 1) {
  omp_set_num_threads(nthreads);
  std::vector<std::string> processed_sigstrings;
  #pragma omp parallel
  {
    std::vector<std::string> local_processed_sigstrings;
    #pragma omp for nowait
    for (int i = 0; i < sigs.size(); ++i) {
      const auto& sigstring = sigs[i];
      ReadData readData;
      fillSigString(readData, readLayout, sigstring, positionFuncMap, verbose);
      concat_scan(readData, verbose);
      pos_scan(readData, positionFuncMap, readLayout, verbose);
      std::vector<std::string> generated_sigstrings = generateSigString(readData);
      filter_sigstring(generated_sigstrings);
      local_processed_sigstrings.insert(local_processed_sigstrings.end(), 
        generated_sigstrings.begin(), 
        generated_sigstrings.end());
    }
    #pragma omp critical
      processed_sigstrings.insert(processed_sigstrings.end(), local_processed_sigstrings.begin(), local_processed_sigstrings.end());
  }
  return wrap(processed_sigstrings);
}

std::vector<std::string> process_sigstrings_cpp(const ReadLayout& readLayout, 
  const std::vector<std::string>& sigs, const PositionFuncMap& positionFuncMap, bool verbose, int nthreads = 1) {
  omp_set_num_threads(nthreads);
  std::vector<std::string> processed_sigstrings;
#pragma omp parallel
{
  std::vector<std::string> local_processed_sigstrings;
#pragma omp for nowait
  for (int i = 0; i < sigs.size(); ++i) {
    const auto& sigstring = sigs[i];
    ReadData readData;
    fillSigString(readData, readLayout, sigstring, positionFuncMap, verbose);
    concat_scan(readData, verbose);
    pos_scan(readData, positionFuncMap, readLayout, verbose);
    std::vector<std::string> generated_sigstrings = generateSigString(readData);
    filter_sigstring(generated_sigstrings);
    local_processed_sigstrings.insert(local_processed_sigstrings.end(), 
      generated_sigstrings.begin(), 
      generated_sigstrings.end());
  }
#pragma omp critical
  processed_sigstrings.insert(processed_sigstrings.end(), local_processed_sigstrings.begin(), local_processed_sigstrings.end());
}
return processed_sigstrings;
}

// [[Rcpp::export]]
Rcpp::CharacterVector sigrun(const Rcpp::DataFrame& read_layout, const Rcpp::DataFrame& misalignment_threshold, 
  const Rcpp::StringVector& sigstrings, int nthreads = 1, bool verbose = false) {
  container = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  VarScan(container, verbose);
  positionFuncMap = createPositionFunctionMap(container, verbose);
  std::vector<std::string> sigs = Rcpp::as<std::vector<std::string>>(sigstrings);
  return process_sigstrings(container, sigs, positionFuncMap, verbose, nthreads);
}