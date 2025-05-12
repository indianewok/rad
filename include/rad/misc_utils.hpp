#pragma once
#include "rad_headers.h"

namespace basic_stats {
    double mean(const std::vector<double>& x) {
        if (x.empty()) return 0.0;
        return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    }

    double sd(const std::vector<double>& x) {
        if (x.size() <= 1) return 0.0;
        double m = mean(x);
        double sum = 0.0;
        for (const auto& val : x) {
            double diff = val - m;
            sum += diff * diff;
        }
        return std::sqrt(sum / (x.size() - 1));
    }
};

namespace path_utils {
    namespace bfs = boost::filesystem;
    namespace bdl = boost::dll;
  
    inline bfs::path find_resource_root() {
      // 1) Where the exe lives
      bfs::path exe = bdl::program_location();
      bfs::path dir = exe.parent_path();
  
      // 2) Climb upwards looking for a "resources" folder
      while (true) {
        bfs::path cand = dir / "resources";
        if (bfs::exists(cand) && bfs::is_directory(cand)) {
          return bfs::canonical(cand);
        }
        if (!dir.has_parent_path()) break;
        dir = dir.parent_path();
      }
  
      // 3) Fallback: ./resources in your CWD
      bfs::path alt = bfs::current_path() / "resources";
      if (bfs::exists(alt) && bfs::is_directory(alt)) {
        return bfs::canonical(alt);
      }
  
      throw std::runtime_error("Cannot locate 'resources' folder up the tree or in CWD");
    }
};

namespace seq_utils {

    std::string revcomp(const std::string& seq) {
        if (seq.empty()) return seq;
        std::string rc = seq;
        std::reverse(rc.begin(), rc.end());
        std::transform(rc.begin(), rc.end(), rc.begin(), [](char c) {
            switch (c) {
                case 'A': return 'T';
                case 'T': return 'A';
                case 'G': return 'C';
                case 'C': return 'G';
                case 'N': return 'N';
                case 'n': return 'n';
                default: return c;
            }
        });
        return rc;
    }
    
    std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\n\r\f\v");
        if (start == std::string::npos)
            return "";
        size_t end = s.find_last_not_of(" \t\n\r\f\v");
        return s.substr(start, end - start + 1);
    }

    bool is_rc(const std::string& str) {
        return str.size() >= 3 && str.substr(0, 3) == "rc_";
    }

    std::string remove_rc(const std::string& str) {
        if (str.size() >= 3 && str.substr(0, 3) == "rc_") {
            return str.substr(3); // Return string without the "rc_" prefix
        }
        return str; // Return original string if no prefix found
    }

    std::string substr_w_padding(const std::string& str, int start, int stop, int pad) {
        int start_pos = std::max(1, start - pad) -1 ;
        int stop_pos = std::min(static_cast<int>(str.length()), stop + pad);
        int length = stop_pos - start_pos;
        if(length <= 0) {
            return ""; // Return empty string if length is non-positive
        }
        return str.substr(start_pos, length);
    }

    static std::vector<std::pair<size_t,std::string>>extract(const std::string& seq, const std::string& pattern, size_t length, int offset){
        std::vector<std::pair<size_t,std::string>> results;
        if (pattern.empty() || length == 0 || seq.empty())
            return results;
        size_t patLen = pattern.size();
        // search for every match of pattern
        for (size_t pos = seq.find(pattern);
            pos != std::string::npos;
            pos = seq.find(pattern, pos + 1))
        {
            // compute start index (can be before or after the match)
            ptrdiff_t start = offset >= 0 ? static_cast<ptrdiff_t>(pos + patLen) + offset  : static_cast<ptrdiff_t>(pos) + offset;

            // ensure we stay in bounds
            if (start >= 0
            && static_cast<size_t>(start) + length <= seq.size())
            {
                results.emplace_back(
                    pos,
                    seq.substr(static_cast<size_t>(start), length)
                );
            }
        }
        return results;
    }

    std::vector<std::string> circ_kmerize(const std::string& seq, size_t k = 16) {
        if (k == 0) 
            throw std::invalid_argument("k must be > 0");
        size_t n = seq.size();
        if (n == 0) 
            return {};  // nothing to do on an empty sequence
        // Case 1: normal k‑mers
        if (n >= k) {
            std::vector<std::string> kmers;
            kmers.reserve(n - k + 1);
            for (size_t i = 0; i + k <= n; ++i) {
                kmers.emplace_back(seq.substr(i, k));
            }
            return kmers;
        }
    
        // Case 2: k > n -> circular wrap
        std::vector<std::string> kmers;
        kmers.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            std::string kmer;
            kmer.reserve(k);
            for (size_t j = 0; j < k; ++j) {
                // wrap index via modulo
                kmer.push_back(seq[(i + j) % n]);
            }
            kmers.emplace_back(std::move(kmer));
        }
        return kmers;
    }

};

namespace streaming_utils { 
    std::vector<std::string> import_text(const std::string &path, size_t max_lines = SIZE_MAX) {
        // 1) open the right kind of stream
        std::unique_ptr<std::istream> in;
        if (path.size() > 3 && path.substr(path.size() - 3) == ".gz") {
            in = std::make_unique<igzstream>(path.c_str());
        } else {
            in = std::make_unique<std::ifstream>(path);
        }
        if (!in || !*in) {
            throw std::runtime_error("Failed to open file: " + path);
        }
        // 2) read lines
        std::vector<std::string> out;
        out.reserve(std::min(max_lines, size_t(10)));
        std::string line;
        size_t count = 0;
        while (count < max_lines && std::getline(*in, line)) {
            out.push_back(std::move(line));
            ++count;
        }
        return out;
    }
};

namespace whitelist_utils {
    // Map of kit names to their respective whitelist paths
    static std::unordered_map<std::string,std::string> kit_wl_paths = {
        //splitseq barcodes
        {"splitseq_bc1", "resources/wl/splitseq_bc1_bitlist.csv.gz"},
        {"splitseq_bc2", "resources/wl/splitseq_bc2_bitlist.csv.gz"},
        //cellranger kits and whitelists--this is just the total list that I could pull together

        //10x 3' kits
        { "10x_3v1", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_3v2", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_3v3", "resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz" },
        { "10x_3v3.1", "resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz" },
        { "10x_3HTv3.1","resources/wl/3M-february-2018-3v3.txt_bitlist.csv.gz" },
        { "10x_3v4", "resources/wl/3M-3pgex-may-2023.txt_bitlist.csv.gz" },

        //10x 5' kits
        { "10x_5v1", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_5v2", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_5HTv2", "resources/wl/737K-august-2016_bitlist.csv.gz" },
        { "10x_5v3", "resources/wl/3M-5pgex-jan-2023.txt_bitlist.csv.gz" },

        //10x Visium kits (and Xenium? need to check)
        { "10x_Vis_V1", "resources/wl/visium-v1_v2_bitlist.csv.gz" },
        { "10x_Vis_V2", "resources/wl/visium-v1_v2_bitlist.csv.gz" },
        { "10x_Vis_V3", "resources/wl/visium-v3_v4_bitlist.csv.gz" },
        { "10x_Vis_V4", "resources/wl/visium-v3_v4_bitlist.csv.gz" },
        { "10x_Vis_V5", "resources/wl/visium-v5_bitlist.csv.gz" }
    };

    // split “kit:path” or “kit:kit” or “path” into 1–2 specs
    static std::vector<std::string> parse_whitelist_specs(std::string const &field) {
        std::vector<std::string> specs;
        std::stringstream ss(field);
        std::string part;
        while (std::getline(ss, part, ':')) {
            if (!part.empty()) specs.push_back(part);
        }
        // keep 1 or 2
        if (specs.size() > 2) specs.resize(2);
        return specs;
    }

    inline std::string kit_to_path(const std::string &field) {
        namespace bfs = boost::filesystem;
    
        // 1) Split into up to two specs
        auto specs = whitelist_utils::parse_whitelist_specs(field);
        std::vector<std::string> resolved;
        resolved.reserve(specs.size());
    
        for (auto &spec : specs) {
            // 2a) if it’s a known kit key, map into resources/wl/…
            if (auto it = kit_wl_paths.find(spec); it != kit_wl_paths.end()) {
                bfs::path rel(it->second); // e.g. "resources/wl/foo.csv.gz"
    
                // strip leading "resources/" if present
                auto iter = rel.begin();
                if (iter!=rel.end() && *iter=="resources") ++iter;
    
                // recombine the remainder
                bfs::path sub;
                for (; iter!=rel.end(); ++iter) sub /= *iter;
    
                // now prepend the true resources folder
                bfs::path root = path_utils::find_resource_root();
                bfs::path full = root / sub;
    
                if (!bfs::exists(full))
                    throw std::runtime_error("Whitelist resource not found: " + full.string());
                resolved.push_back(bfs::canonical(full).string());
            }
            else {
                // 2b) otherwise treat spec as a path
                bfs::path p(spec);
                if (p.is_relative()) p = bfs::absolute(p);
                if (!bfs::exists(p))
                    throw std::runtime_error("Whitelist file not found: " + p.string());
                resolved.push_back(bfs::canonical(p).string());
            }
        }
    
        // 3) re-join with ":" if there were two pieces
        if (resolved.size() == 2) {
            return resolved[0] + ":" + resolved[1];
        } else {
            return resolved[0];
        }
    }
    /*
    inline std::string kit_to_path(const std::string &field) {
        namespace bfs = boost::filesystem;
        // 1) Peel off any "prefix:" to isolate the kit key
        std::string kit = field;
        if (auto pos = field.find(':'); pos!=std::string::npos) {
          kit = field.substr(pos+1);
        }
    
        // 2) If it's one of our known kits, build resources/wl/... from the relative map entry
        auto it = kit_wl_paths.find(kit);
        if (it != kit_wl_paths.end()) {
          bfs::path orig_rel(it->second); // e.g. "resources/wl/foo.csv.gz" or "wl/foo.csv.gz"
    
          // strip leading "resources/" if present
          auto iter = orig_rel.begin();
          if (iter!=orig_rel.end() && *iter=="resources") {
            ++iter;
          }

          // recombine what remains
          bfs::path rel;
          for (; iter!=orig_rel.end(); ++iter) {
            rel /= *iter;
          }
    
          // now join against the real resources folder
          bfs::path root = path_utils::find_resource_root();
          bfs::path full = root / rel;
    
          if (!bfs::exists(full))
            throw std::runtime_error("Whitelist resource not found: " + full.string());
          return bfs::canonical(full).string();
        }
    
        // 3) Otherwise treat the field as a filesystem path and canonicalize it
        bfs::path p(field);
        if (p.is_relative()) {
          p = bfs::absolute(p);
        }
        if (!bfs::exists(p))
          throw std::runtime_error("Whitelist file not found: " + p.string());
        return bfs::canonical(p).string();
      }
    std::string kit_to_path(const std::string &field) {
        auto kit = field;
        if (auto p = field.find(':'); p != std::string::npos){
            kit = field.substr(p+1);
        }
        if (auto it = kit_wl_paths.find(kit); it != kit_wl_paths.end()){
            return it->second;
        }
        return field; 
    }
    */

    // Check if it's a kit
    bool is_kit(const std::string &name){
        if(kit_wl_paths.find(name) != kit_wl_paths.end()){
            return true;
        }
        return false;
    }

    //check if the whitelist we're importing is a bitlist or a whitelist
    bool check_if_bitlist(const std::string &field, bool verbose, size_t N = 10) {
        const std::string path = kit_to_path(field);

        //skip the header and keep it moving
        auto lines = streaming_utils::import_text(path, N + 1);
        if (lines.empty()) {
            if(verbose) std::cerr << "[check_if_bitlist] ERROR: no lines read from “" << path << "”\n";
            return false;
        }
        if(verbose) std::cout << "[check_if_bitlist] debug: first " << lines.size() << " raw lines from “" << path << "”:\n";

        for (size_t i = 0; i < lines.size(); ++i) {
            if(verbose) std::cout << "  [" << i << "]: " << lines[i] << "\n";
        }

        lines.erase(lines.begin());

        if(verbose) std::cout << "[check_if_bitlist] debug: next up to " << N << " content lines (after header):\n";

        for (size_t i = 0; i < lines.size() && i < N; ++i) {
            if(verbose) std::cout << "  [" << i << "]: " << lines[i] << "\n";
        }
    
        // 2) extract up to N tokens (before first comma/tab)
        std::vector<std::string> tokens;
        tokens.reserve(N);
        for (auto &ln : lines) {
          if (tokens.size() >= N) break;
          if (ln.empty()) continue;
          auto p = ln.find_first_of(",\t");
          tokens.push_back(p == std::string::npos ? ln : ln.substr(0, p));
        }
    
        // nothing real to check?
        if (tokens.empty()){
            return false;
        }
        // 3) ensure every character in each token is a digit
        for (auto &t : tokens) {
          for (unsigned char c : t) {
            if (!std::isdigit(c))
              return false;
          }
        }
        return true;  // all N tokens passed the digit test
      }

    inline const std::string& get_whitelist_path(const std::string &kit) {
        auto it = kit_wl_paths.find(kit);
        if (it == kit_wl_paths.end())
            throw std::invalid_argument("Unknown kit name: " + kit);
        return it->second;
    }

    inline void set_whitelist_path(const std::string &kit, const std::string &path) {
        kit_wl_paths[kit] = path;
    }

    inline bool remove_whitelist_path(const std::string &kit) {
        return kit_wl_paths.erase(kit) > 0;
    }

};

namespace config_utils {

    //config for read layout
    static std::unordered_map<std::string, std::string> layout_files = {
        {"five_prime", "resources/read_layout/five_prime_read_layout.csv" },
        {"three_prime", "resources/read_layout/three_prime_read_layout.csv" },
        {"splitseq", "resources/read_layout/splitseq_read_layout.csv" }
    };

    // check if the layout is a custom one
    inline bool check_if_custom_rl(const std::string &type) {
        return layout_files.find(type) == layout_files.end();
    }
    
    // get read layout path from pre-existing elements
    /*
    inline const std::string& get_read_layout(const std::string &type) {
        auto it = layout_files.find(type);
        if (it == layout_files.end())
            throw std::invalid_argument("Unknown layout type: " + type);
        return it->second;
    }
    */
   inline std::string get_read_layout(const std::string &type) {
    namespace bfs = boost::filesystem;
    // 1) If `type` looks like an absolute or contains a path‐sep, treat it as a real file:
    bfs::path p(type);
    if (p.is_absolute() ||
        type.find(bfs::path::preferred_separator) != std::string::npos)
    {
      if (!bfs::exists(p))
        throw std::runtime_error("Layout file not found: " + type);
      return bfs::canonical(p).string();
    }

    // 2) Otherwise look up our map
    auto it = layout_files.find(type);
    if (it == layout_files.end())
      throw std::invalid_argument("Unknown layout key: " + type);

    // 3) Take the stored relative‐to‐resources path
    //    e.g. it->second == "splitseq_read_layout.csv"   OR
    //         it->second == "read_layout/splitseq_read_layout.csv" OR
    //         it->second == "resources/read_layout/splitseq_read_layout.csv"
    bfs::path orig_rel(it->second);

    // strip any leading "resources/" component:
    auto iter = orig_rel.begin();
    if (iter!=orig_rel.end() && *iter == "resources") {
      ++iter;
    }
    // now build a new path from whatever remains
    bfs::path rel;
    for (; iter!=orig_rel.end(); ++iter) {
      rel /= *iter;
    }

    // 4) Join against the real resources folder
    bfs::path root = path_utils::find_resource_root();
    bfs::path full = root / rel;

    if (!bfs::exists(full))
      throw std::runtime_error("Missing layout resource: " + full.string());
    return bfs::canonical(full).string();
  }

    inline void save_read_layout(const std::string &type, const std::string &path) {
        layout_files[type] = path;
    }
    
    inline bool remove_read_layout(const std::string &type) {
        return layout_files.erase(type) > 0;
    }

};