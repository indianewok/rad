#include "include/rad/rad_headers.h"

// ─────────────────────────────────────────────────────────────────────────────
// CLI + helpers (unchanged)
// ─────────────────────────────────────────────────────────────────────────────
struct Options {
    std::string input  = "-";      // "-" = stdin
    std::string output = "-";      // "-" = stdout (uncompressed)
    int threads        = 0;        // 0 = OMP default
    int chunk_size     = 50000;    // reads per chunk
    std::string delim  = ":";      // separator between QNAME/CB/UB
    bool keep_umi      = true;     // include UB if present
};

static void usage(const char* prog) {
    std::cerr
      << "Usage: " << prog << " [options]\n"
      << "  -i, --input   FILE     Input FASTQ(.gz) or '-' for stdin (file or directory ok)\n"
      << "  -o, --output  FILE     Output FASTQ(.gz) or '-' for stdout\n"
      << "  -t, --threads INT      OpenMP threads (default: OMP default)\n"
      << "  -k, --chunk   INT      Reads per chunk (default: 50000)\n"
      << "  -d, --delimiter STR    Delimiter between fields (default: ':')\n"
      << "      --keep-umi BOOL    Include UB in ID (default: true)\n";
}

static inline bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() && std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}

struct FqRec {
    std::string name;  // kseq name (no leading '@'), stops at first space/tab
    std::string seq;
    std::string qual;  // if absent, we synthesize 'I'*
    std::string aux;   // everything after the first space/tab from header line (tabs preserved)
};

// Find CB/UB values from a tab-separated tail using zero-copy string_view.
static inline std::string_view find_tag(std::string_view tail, std::string_view key) {
    size_t pos = tail.find(key);
    if (pos == std::string_view::npos) return {};
    pos += key.size();
    size_t end = tail.find('\t', pos);
    return (end == std::string_view::npos) ? tail.substr(pos) : tail.substr(pos, end - pos);
}

// Build new ID: QNAME[:CB][:UB] with a single allocation.
static inline std::string make_new_id(const std::string& qname,
                                      std::string_view cb,
                                      std::string_view ub,
                                      const std::string& delim,
                                      bool keep_umi)
{
    size_t need = qname.size();
    if (!cb.empty()) need += delim.size() + cb.size();
    if (keep_umi && !ub.empty()) need += delim.size() + ub.size();

    std::string out;
    out.reserve(need);
    out.append(qname);
    if (!cb.empty()) { out.append(delim); out.append(cb.data(), cb.size()); }
    if (keep_umi && !ub.empty()) { out.append(delim); out.append(ub.data(), ub.size()); }
    return out;
}

// Output sink: plain or .gz, or stdout (zlib for .gz)
class FastqWriter {
public:
    explicit FastqWriter(const std::string& path) {
        if (path == "-" || path.empty()) { use_stdout = true; return; }
        use_gz = ends_with(path, ".gz");
        if (use_gz) {
            gz = gzopen(path.c_str(), "wb");
            if (!gz) throw std::runtime_error("Failed to open gzip output: " + path);
        } else {
            out.open(path, std::ios::out | std::ios::binary);
            if (!out) throw std::runtime_error("Failed to open output: " + path);
        }
    }
    ~FastqWriter() { if (gz) gzclose(gz); if (out.is_open()) out.close(); }

    inline void write_record(const std::string& name,
                             const std::string& seq,
                             const std::string& qual) {
        if (use_stdout) {
            std::cout << '@' << name << '\n' << seq << '\n' << "+\n" << qual << '\n';
            return;
        }
        if (use_gz) {
            gzputc(gz, '@'); gzwrite(gz, name.data(), (unsigned)name.size()); gzputc(gz, '\n');
            gzwrite(gz, seq.data(),  (unsigned)seq.size());  gzputc(gz, '\n');
            gzputs(gz, "+\n");
            gzwrite(gz, qual.data(), (unsigned)qual.size()); gzputc(gz, '\n');
        } else {
            out << '@' << name << '\n' << seq << '\n' << "+\n" << qual << '\n';
        }
    }
private:
    bool use_stdout = false, use_gz = false;
    gzFile gz = nullptr;
    std::ofstream out;
};

// ─────────────────────────────────────────────────────────────────────────────
// MAIN: uses your file_streaming + read_streaming to support pigz/fd + gz/zlib
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    Options opt;
    static struct option long_opts[] = {
        {"input",     required_argument, 0, 'i'},
        {"output",    required_argument, 0, 'o'},
        {"threads",   required_argument, 0, 't'},
        {"chunk",     required_argument, 0, 'k'},
        {"delimiter", required_argument, 0, 'd'},
        {"keep-umi",  required_argument, 0,  1 },
        {"help",      no_argument,       0, 'h'},
        {0,0,0,0}
    };

    int c, li = 0;
    while ((c = getopt_long(argc, argv, "i:o:t:k:d:h", long_opts, &li)) != -1) {
        switch (c) {
            case 'i': opt.input  = optarg; break;
            case 'o': opt.output = optarg; break;
            case 't': opt.threads = std::max(0, atoi(optarg)); break;
            case 'k': opt.chunk_size = std::max(1, atoi(optarg)); break;
            case 'd': opt.delim = optarg; break;
            case 'h': usage(argv[0]); return 0;
            case 1: {
                std::string v = optarg ? std::string(optarg) : "true";
                std::transform(v.begin(), v.end(), v.begin(), ::tolower);
                opt.keep_umi = !(v == "0" || v == "false" || v == "no");
                break;
            }
            default: usage(argv[0]); return 1;
        }
    }

    if (opt.threads > 0) omp_set_num_threads(opt.threads);

    // If input is stdin ("-"), we cannot run pigz; fall back to zlib path via your reader.
    // Otherwise, your file_streaming picks pigz for .gz if available, else gz/zlib, else plain.
    std::unique_ptr<file_streaming> files_ptr;
    std::unique_ptr<read_streaming> reader_ptr;

    try {
        files_ptr = std::make_unique<file_streaming>(opt.input, std::max(1, opt.threads));
        reader_ptr = std::make_unique<read_streaming>(*files_ptr);
    } catch (const std::exception& e) {
        std::cerr << "Error opening input: " << e.what() << "\n";
        return 1;
    }

    FastqWriter writer(opt.output);

    std::vector<FqRec> chunk;
    chunk.reserve(opt.chunk_size);

    std::atomic<size_t> total{0};

    // Stream → chunk → transform headers (parallel) → ordered write
    while (true) {
        chunk.clear();

        // Fill a chunk from the reader
        for (int i = 0; i < opt.chunk_size; ++i) {
            auto rec = reader_ptr->next_sequence();
            if (!rec) break;

            FqRec r;
            r.name = std::move(rec->id);
            r.seq  = std::move(rec->seq);

            if (rec->is_fastq && !rec->qual.empty()) {
                r.qual = std::move(rec->qual);
            } else {
                continue;
            }

            r.aux = std::move(rec->comment); // comment contains the tabbed tags tail
            chunk.push_back(std::move(r));
        }
        if (chunk.empty()) break;

        // Transform headers in parallel (read-only over aux/name)
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < chunk.size(); ++i) {
            auto& rec = chunk[i];

            // Tail is tab-separated tags
            std::string_view tail(rec.aux.data(), rec.aux.size());

            // Zero-copy extract CB/UB
            std::string_view cb = find_tag(tail, "CB:Z:");
            std::string_view ub = find_tag(tail, "UB:Z:");

            // Build final ID; if tags missing/empty, they’re simply omitted
            std::string newid = make_new_id(rec.name, cb, ub, opt.delim, opt.keep_umi);
            rec.name.swap(newid);
        }

        // Ordered write
        for (const auto& r : chunk) writer.write_record(r.name, r.seq, r.qual);

        total += chunk.size();
        if (total % (size_t(10) * opt.chunk_size) == 0) {
            std::cerr << "[fq_header_rewriter] processed " << total.load() << " reads\n";
        }
    }

    std::cerr << "[fq_header_rewriter] done. total reads: " << total.load() << "\n";
    return 0;
}
