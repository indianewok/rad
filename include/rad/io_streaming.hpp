#pragma once
#include "rad_headers.h"

//simplest solution i could figure out 
namespace kseq_fd {
    #ifndef KSEQ_INIT_INT_READ
    #define KSEQ_INIT_INT_READ
    KSEQ_INIT(int, read)
    #endif
}

namespace kseq_gz {
    #ifndef KSEQ_INIT_GZFILE_GZREAD  
    #define KSEQ_INIT_GZFILE_GZREAD
    KSEQ_INIT(gzFile, gzread)
    #endif
}

/**
 * @brief Parallel decompression reading using pigz
 * 
 * Manages pigz subprocess for multi-threaded decompression of gzipped files.
 */
class pigz_reading {
public:
    /**
     * @brief Open a pigz decompression pipe for reading
     * @param input_path Path to compressed file (.gz, .fq.gz, .fastq.gz, etc.)
     * @param threads Number of decompression threads (default: 4)
     * @return FILE* pointer to pipe, or nullptr on failure
     */
    inline FILE* open_pigz_read_pipe(const std::string& input_path, int threads = 4) {
        // Allow disabling via environment variable
        if (const char* no = std::getenv("RAD_NO_PIGZ"); no && *no) {
            return nullptr;
        }

        #ifndef PIGZ_PATH
        #  define PIGZ_PATH "pigz"
        #endif
        
        std::string pigz = PIGZ_PATH;
        if (const char* env = std::getenv("RAD_PIGZ"); env && *env) {
            pigz = env;
        }

        // Use pigz -dc for decompression to stdout
        // -d = decompress, -c = write to stdout, -p = threads
        std::string cmd = "'" + pigz + "' -dc -p " + std::to_string(threads) + 
                         " '" + input_path + "'";
        
        FILE* fp = popen(cmd.c_str(), "r");
        if (!fp) {
            return nullptr;
        }
        
        // Large read buffer (16 MB) to reduce read() syscalls
        static thread_local std::vector<char> read_buf(16u << 20);
        setvbuf(fp, read_buf.data(), _IOFBF, read_buf.size());
        
        return fp;
    }

    /**
     * @brief Close pigz decompression pipe
     * @param fp FILE pointer from open_pigz_read_pipe
     * @return Exit status of pigz process
     */
    inline int close_pigz_read_pipe(FILE* fp) {
        if (!fp) return -1;
        return pclose(fp);
    }

    /**
     * @brief Check if pigz is available on the system
     * @return true if pigz is found in PATH
     */
    static bool is_pigz_available() {
        #ifndef PIGZ_PATH
        #  define PIGZ_PATH "pigz"
        #endif
        
        std::string pigz = PIGZ_PATH;
        if (const char* env = std::getenv("RAD_PIGZ"); env && *env) {
            pigz = env;
        }
        
        std::string cmd = "which '" + pigz + "' > /dev/null 2>&1";
        return system(cmd.c_str()) == 0;
    }
};

/**
 * @brief File streaming with automatic parallel decompression
 * 
 * Automatically uses pigz for .gz files when available, falls back to gzip.
 * Handles both single files and directories of FASTQ/FASTA files.
 */
class file_streaming {
public:
    struct file_pointer {
        FILE* fp_pigz;      // Pigz pipe handle or regular file handle
        gzFile fp_gzip;     // Fallback gzip handle
        kseq_fd::kseq_t* seq_fd;     // kseq for FILE* (pigz path or uncompressed)
        kseq_gz::kseq_t* seq_gz;     // kseq for gzFile (fallback path)
        int fd;             // File descriptor for kseq
        std::string path;
        bool using_pigz;
        pigz_reading pigz_reader;
        
        inline kseq_fd::kseq_t* get_seq_fd() { 
            return seq_fd; 
        }

        inline kseq_gz::kseq_t* get_seq_gz() { 
            return seq_gz; 
        }

        // one-step read that picks the right backend
        inline int kseq_read_any() {
            if (using_pigz) {
                return seq_fd ? kseq_fd::kseq_read(seq_fd) : -1;  // -1 = EOF, -2 = bad FASTQ
            } else {
                return seq_gz ? kseq_gz::kseq_read(seq_gz) : -1;
            }
        }

        // unified field accessors (avoid duplicating branches)
        inline const char* name_s() const {
            return using_pigz ? seq_fd->name.s : seq_gz->name.s;
        }
        inline const char* comment_s() const {
            return using_pigz ? (seq_fd->comment.l ? seq_fd->comment.s : nullptr)
                            : (seq_gz->comment.l ? seq_gz->comment.s : nullptr);
        }
        inline const char* seq_s() const {
            return using_pigz ? seq_fd->seq.s : seq_gz->seq.s;
        }
        inline int qual_len() const {
            return using_pigz ? seq_fd->qual.l : seq_gz->qual.l;
        }
        inline const char* qual_s() const {
            return using_pigz ? (seq_fd->qual.l ? seq_fd->qual.s : nullptr)
                            : (seq_gz->qual.l ? seq_gz->qual.s : nullptr);
        }

        file_pointer(const std::string& file_path, int pigz_threads = 4) 
            : fp_pigz(nullptr), fp_gzip(nullptr), seq_fd(nullptr), seq_gz(nullptr),
              fd(-1), path(file_path), using_pigz(false) {
            
            // Check if file is compressed
            bool is_compressed = has_gz_extension(file_path);
            
            if (is_compressed) {
                // Try pigz first for parallel decompression
                if (pigz_reading::is_pigz_available()) {
                    fp_pigz = pigz_reader.open_pigz_read_pipe(file_path, pigz_threads);
                    if (fp_pigz) {
                        fd = fileno(fp_pigz);
                        seq_fd = kseq_fd::kseq_init(fd);
                        using_pigz = true;
                        return;
                    }
                }

                // Fallback to single-threaded gzip
                fp_gzip = gzopen(file_path.c_str(), "r");
                if (fp_gzip) {
                    gzbuffer(fp_gzip, 1 << 22);  // 4 MiB buffer
                    seq_gz = kseq_gz::kseq_init(fp_gzip);
                    using_pigz = false;
                }
            } else {
                // Uncompressed file - use regular file descriptor
                fp_pigz = fopen(file_path.c_str(), "r");
                if (fp_pigz) {
                    fd = fileno(fp_pigz);
                    seq_fd =  kseq_fd::kseq_init(fd);
                    setvbuf(fp_pigz, nullptr, _IOFBF, 1 << 22);  // 4 MiB buffer
                    using_pigz = true;  // Using FILE*, not gzFile
                }
            }
        }

        ~file_pointer() {
            if (seq_fd) kseq_destroy(seq_fd);
            if (seq_gz) kseq_destroy(seq_gz);
            if (using_pigz && fp_pigz) {
                pigz_reader.close_pigz_read_pipe(fp_pigz);
            } else if (fp_gzip) {
                gzclose(fp_gzip);
            }
        }

        bool is_valid() const {
            return (using_pigz && fp_pigz && seq_fd) || (fp_gzip && seq_gz);
        }

    private:
        static bool has_gz_extension(const std::string& path) {
            return path.size() >= 3 && 
                   path.substr(path.size() - 3) == ".gz";
        }
    };

    std::unique_ptr<file_pointer> current_file;
    std::filesystem::directory_iterator dir_iter;
    bool dir_status;
    std::string path;
    int pigz_threads;

    explicit file_streaming(const std::string& input_path, int threads = 4) 
        : path(input_path), pigz_threads(threads) {
        dir_status = is_directory(path);
        if (dir_status) {
            dir_iter = std::filesystem::directory_iterator(path);
            open_next_file();
        } else if (is_fastqa(path)) {
            current_file = std::make_unique<file_pointer>(path, pigz_threads);
            if (!current_file->is_valid()) {
                throw std::runtime_error("Failed to open sequence file: " + path);
            }
        } else {
            throw std::invalid_argument("Invalid input path: " + path);
        }
    }

    bool open_next_file() {
        if (!dir_status) return false;
        while (dir_iter != std::filesystem::directory_iterator()) {
            const auto& entry = *dir_iter++;
            if (!is_file(get_file_path(entry))) continue;
            
            std::string file_path = entry.path().string();
            if (is_fastqa(file_path)) {
                current_file = std::make_unique<file_pointer>(file_path, pigz_threads);
                if (current_file->is_valid()) {
                    return true;
                }
            }
        }
        return false;
    }

    static bool has_suffix(const std::string& str, const std::string& suffix) {
        return str.size() >= suffix.size() &&
               str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    static bool is_fastqa(const std::string& path) {
        static const std::vector<std::string> extensions = {
            ".fastq", ".fq", ".fasta", ".fa",
            ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz"
        };
        
        for (const auto& ext : extensions) {
            if (has_suffix(path, ext)) return true;
        }
        return false;
    }

    static bool is_directory(const std::string& path) {
        return std::filesystem::is_directory(path);
    }

    static bool is_file(const std::string& path) {
        return std::filesystem::is_regular_file(path);
    }

    static std::string get_file_path(const std::filesystem::directory_entry& entry) {
        return entry.path().string();
    }

    std::string current_file_path() const {
        return current_file ? current_file->path : "";
    }
};

/**
 * @brief Read streaming for FASTQ/FASTA sequences
 */
class read_streaming {
public:
    struct sequence {
        std::string id;
        std::string comment;
        std::string seq;
        std::string qual;
        bool is_fastq;
    };

    file_streaming& files;
    
    explicit read_streaming(file_streaming& f) : files(f) {}
    
    std::optional<sequence> next_sequence() {
        while (files.current_file || files.open_next_file()) {
            auto& fp = files.current_file;
            int l = fp->kseq_read_any();
            if (l >= 0) {
                sequence r;
                r.id = fp->name_s();
                r.seq = fp->seq_s();
                r.is_fastq = fp->qual_len() > 0;

                if (const char* c = fp->comment_s()){
                    r.comment = c;
                }
                if (r.is_fastq){
                    r.qual = fp->qual_s();
                }
                return r;
            }

            // Move to next file if in directory mode
            if (files.dir_status) {
                files.current_file.reset();
            } else {
                break;
            }
        }
        return std::nullopt;
    }
};

// ====== minimal memory probe (no extra headers/files) ======
namespace readmem_utils {
    // enable/disable logging quickly
    constexpr bool ENABLE = true;

    // --- size helpers using capacity() (no fudge) ---
    static inline size_t bytes_of(const std::string& s) noexcept { return s.capacity(); }

    // specialize for your record type: read_streaming::sequence
    static inline size_t bytes_of_rec(const read_streaming::sequence& r) noexcept {
        return bytes_of(r.id) + bytes_of(r.comment) + bytes_of(r.seq) + bytes_of(r.qual) + sizeof(r);
    }

    template <class Vec>
    static inline size_t bytes_of_vec(const Vec& v) noexcept {
        using T = typename Vec::value_type;
        size_t total = v.capacity() * sizeof(T);     // vector buffer itself
        for (auto const& e : v) total += bytes_of_rec(e); // payloads (strings' capacities)
        return total;
    }

    // --- global in-flight counters (bytes/chunks) ---
    struct InFlight {
        std::atomic<size_t> bytes{0};
        std::atomic<size_t> chunks{0};
        std::atomic<size_t> peak_bytes{0};
        std::atomic<size_t> peak_chunks{0};

        void add(size_t b) noexcept {
            const size_t nb = bytes.fetch_add(b, std::memory_order_relaxed) + b;
            const size_t nc = chunks.fetch_add(1, std::memory_order_relaxed) + 1;
            size_t pb = peak_bytes.load(std::memory_order_relaxed);
            while (nb > pb && !peak_bytes.compare_exchange_weak(pb, nb, std::memory_order_relaxed)) {}
            size_t pc = peak_chunks.load(std::memory_order_relaxed);
            while (nc > pc && !peak_chunks.compare_exchange_weak(pc, nc, std::memory_order_relaxed)) {}
        }
        void sub(size_t b) noexcept {
            bytes.fetch_sub(b, std::memory_order_relaxed);
            chunks.fetch_sub(1, std::memory_order_relaxed);
        }
    };
    inline InFlight& inflight() { static InFlight g; return g; }

    struct ScopedInFlight {
        size_t b{0}; bool active{false};
        explicit ScopedInFlight(size_t bytes) : b(bytes), active(true) { inflight().add(b); }
        ~ScopedInFlight(){ if (active) inflight().sub(b); }
        ScopedInFlight(const ScopedInFlight&) = delete;
        ScopedInFlight& operator=(const ScopedInFlight&) = delete;
    };

    static inline double rss_gib() {
        FILE* f = std::fopen("/proc/self/statm", "r");
        if (!f) return -1.0;
        long long size_pages=0, res_pages=0;
        if (std::fscanf(f, "%lld %lld", &size_pages, &res_pages) != 2) { std::fclose(f); return -1.0; }
        std::fclose(f);
        const double page_kb = double(::sysconf(_SC_PAGESIZE)) / 1024.0;
        const double rss_kb = res_pages * page_kb;
        return rss_kb / (1024.0 * 1024.0);
    }

    // --- tiny logger ---
    static inline void log_chunk(std::size_t chunk_bytes, std::size_t chunk_id=0) {
        if (!ENABLE) return;
        const double mb = double(chunk_bytes) / (1024.0*1024.0);
        const double inflight_gib = double(inflight().bytes.load(std::memory_order_relaxed)) / (1024.0*1024.0*1024.0);
        const std::size_t inflight_chunks = inflight().chunks.load(std::memory_order_relaxed);
        const double rss = rss_gib();
        std::fprintf(stderr,
            "[radmem] chunk=%zu size_mb=%.2f inflight_chunks=%zu inflight_gib=%.3f rss_gib=%.3f\n",
            chunk_id, mb, inflight_chunks, inflight_gib, rss);
    }
}


/**
 * @brief Chunk-based streaming with parallel decompression
 * 
 * Processes input files in chunks, using pigz for parallel decompression
 * and OpenMP tasks for parallel processing.
 */
template<typename T, typename Function> 
class chunk_streaming {
public:
    explicit chunk_streaming(size_t chunk_size, int pigz_threads = 4, int max_in_flight = -1) 
        : chunk_size_(chunk_size), pigz_threads_(pigz_threads), 
          max_in_flight_(max_in_flight), seqs_processed_(0), peak_in_flight_(0) {}

    void modify_chunk_size(size_t new_size) {
        chunk_size_ = new_size;
    }

    size_t get_chunk_size() const {
        return chunk_size_;
    }

    int get_mif() const {
        return max_in_flight_;
    }
  
    void process_chunks(const std::string& input_path, 
                    Function chunk_func, 
                    int num_threads, 
                    int64_t max_seqs = -1) {
            file_streaming files(input_path, pigz_threads_);
            read_streaming reader(files);

            std::mutex reader_mutex;
            std::atomic<bool>   more_data{true};
            std::atomic<size_t> seqs_processed{0};

            #pragma omp parallel num_threads(num_threads)
            {
                while (true) {
                    std::vector<T> chunk;
                    {
                        std::lock_guard<std::mutex> lock(reader_mutex);
                        // if another thread already signaled "no more data", bail
                        if (!more_data.load(std::memory_order_acquire)) {
                            break;
                        }
                        chunk.reserve(chunk_size_);
                        while (chunk.size() < chunk_size_) {
                            if (max_seqs >= 0 && seqs_processed.load(std::memory_order_relaxed) >= static_cast<size_t>(max_seqs)) {
                                more_data.store(false, std::memory_order_release);
                                break;
                            }
                            // try to read
                            auto rec = reader.next_sequence();
                            if (!rec) {
                                more_data.store(false, std::memory_order_release);
                                break;
                            }
                            chunk.push_back(std::move(*rec));
                            seqs_processed.fetch_add(1, std::memory_order_relaxed);
                        }
                    } 
                    if (chunk.empty()) {
                        if (!more_data.load(std::memory_order_acquire)){
                            break;
                        }
                        continue;
                    }

                    const size_t __chunk_bytes = readmem_utils::bytes_of_vec(chunk);
                    readmem_utils::ScopedInFlight __guard(__chunk_bytes);
                    readmem_utils::log_chunk(__chunk_bytes /*, optional_chunk_id */);
                    chunk_func(chunk, input_path);
                }
            }
        }


    private:
    size_t chunk_size_;
    int pigz_threads_;
    int max_in_flight_;  // -1 = unlimited, >0 = limited
    std::atomic<size_t> seqs_processed_;
    std::atomic<int> peak_in_flight_;
};

// ============================================================================
// OUTPUT STREAMING (PIGZ COMPRESSION)
// ============================================================================

class pigz_writing {
public:
    inline FILE* open_pigz_pipe(const std::string& out_path, int threads, int level=1) {
        if (const char* no = std::getenv("RAD_NO_PIGZ"); no && *no) return nullptr;

    #ifndef PIGZ_PATH
    #  define PIGZ_PATH "pigz"
    #endif
        std::string pigz = PIGZ_PATH;
        if (const char* env = std::getenv("RAD_PIGZ"); env && *env) pigz = env;

        std::string cmd = "'" + pigz + "' -c -p " + std::to_string(threads) +
                        " -" + std::to_string(level) + " > '" + out_path + "'";
        FILE* fp = popen(cmd.c_str(), "w");
        if (!fp) return nullptr;
        // 16MB stdio buffer to reduce write() syscalls
        static std::vector<char> buf(16u << 20);
        setvbuf(fp, buf.data(), _IOFBF, buf.size());
        return fp;
    }

    inline void pigz_write(FILE* fp, const char* data, size_t n) {
        if (!n) return;
        size_t w = std::fwrite(data, 1, n, fp);
        if (w != n) {
            throw std::runtime_error("pigz_write: partial write");
        }
    }

    inline int close_pigz_pipe(FILE* fp) {
        if (!fp) return -1;
        std::fflush(fp);
        return pclose(fp);
    }
};

namespace bio = boost::iostreams;
class sigstring_writing {
public:
    enum class format { 
        FASTQA, 
        SIGSTRING, 
        CSV 
    };

private:
    // Common
    format fmt;
    bool compress_;
    bool use_pigz_ = false;

    // --- pigz path ---
    pigz_writing pigz_;
    FILE* pigz_fp_ = nullptr;
    std::string pigz_buf_;
    size_t pigz_flush_threshold_ = (8u << 20); // 8MB

    // --- zlib/boost path ---
    std::ofstream file_;
    bio::filtering_ostream out_; // chain: [gzip?] -> file_

    inline void pigz_flush_if_big_() {
        if (pigz_buf_.size() >= pigz_flush_threshold_) {
            pigz_.pigz_write(pigz_fp_, pigz_buf_.data(), pigz_buf_.size());
            pigz_buf_.clear();
            // pigz_buf_.shrink_to_fit();
        }
    }

    template <typename T>
    inline void write_one_(const T& item) {
        switch (fmt) {
            case format::FASTQA:
                if (use_pigz_) { pigz_buf_.append(item.to_fastqa()); pigz_flush_if_big_(); }
                else           { out_ << item.to_fastqa(); }
                break;
            case format::SIGSTRING:
                if (use_pigz_) { pigz_buf_.append(item.to_sigstring()); pigz_buf_.push_back('\n'); pigz_flush_if_big_(); }
                else           { out_ << item.to_sigstring() << "\n"; }
                break;
            case format::CSV:
                if (use_pigz_) { pigz_buf_.append(item.to_csv()); pigz_flush_if_big_(); }
                else           { out_ << item.to_csv(); }
                break;
        }
    }

public:
    /// @param output_path   base filename ("myout"), will add ".gz" if compress==true and it isn't present
    /// @param output_format which of the three formats to emit
    /// @param compress      if true, writes gzipped output (pigz if available; else zlib/boost)
    /// @param append        if true, opens in append mode; otherwise truncates
    /// @param pigz_threads  pigz worker threads (ignored if not using pigz)
    /// @param level         gzip compression level [1..9]

    sigstring_writing(std::string output_path, format output_format, bool compress, bool append, int pigz_threads, int level=1)
        : fmt(output_format), compress_(compress)
    {
        std::ios_base::openmode mode =
            std::ios::out | std::ios::binary | (append ? std::ios::app : std::ios::trunc);

        if (compress_) {
            // ensure .gz suffix for either pigz or zlib
            if (output_path.size() < 3 || output_path.substr(output_path.size()-3) != ".gz")
                output_path += ".gz";

            // Try pigz first
            pigz_fp_ = pigz_.open_pigz_pipe(output_path, pigz_threads, level);
            if (pigz_fp_) {
                use_pigz_ = true;
                // nothing else to open; pigz writes the file
                return;
            }

            // Fallback: existing Boost gzip chain (single-threaded zlib)
            file_.open(output_path, mode);
            if (!file_.is_open()) throw std::runtime_error("Failed to open file: " + output_path);

            bio::gzip_params params;
            params.level = level;
            out_.push(bio::gzip_compressor(params));
            out_.push(file_);
        } else {
            // Plain text (no compression)
            file_.open(output_path, mode);
            if (!file_.is_open()) throw std::runtime_error("Failed to open file: " + output_path);
            out_.push(file_);
        }
    }

    template<typename T>
    void operator()(const std::vector<T>& chunk) {
        if (chunk.empty()) return;
        if (use_pigz_) {
            // Build big uncompressed buffer for pigz (few syscalls, fewer context switches)
            for (auto const& item : chunk) write_one_(item);
        } else {
            // Your legacy path (ostream through Boost chain)
            for (auto const& item : chunk) write_one_(item);
        }
    }

    // Optional: adjust pigz flush size (in bytes)
    void set_pigz_flush_threshold(size_t bytes) { pigz_flush_threshold_ = bytes; }

    ~sigstring_writing() {
        try {
            if (use_pigz_) {
                // flush tail and close pipe
                if (!pigz_buf_.empty()) {
                    pigz_.pigz_write(pigz_fp_, pigz_buf_.data(), pigz_buf_.size());
                    pigz_buf_.clear();
                }
                (void)pigz_.close_pigz_pipe(pigz_fp_);
                pigz_fp_ = nullptr;
                use_pigz_ = false;
            } else {
                // flush Boost chain + file
                out_.reset();  // pops filters and flushes
                file_.flush();
                file_.close();
            }
        } catch (...) {
        }
    }

    // in sigstring_writing (public)
    void operator()(const std::string& slab) {
        if (slab.empty()) return;
        if (use_pigz_) {
            pigz_buf_.append(slab);
            pigz_flush_if_big_();
        } else {
            out_.write(slab.data(), slab.size());
        }
    }

     // reuses the existing std::string overload
    void operator()(const std::vector<std::string>& slabs) {
        if (slabs.empty()) return;
        for (auto const& s : slabs) {
            (*this)(s); 
        }
    }

};

// ============================================================================
// PARALLEL WRITER
// ============================================================================

class parallel_writer {
public:
    parallel_writer() : running(true), jobs_processed(0), total_write_time_ms(0), total_queue_time_ms(0) {
        writer_thread = std::thread([this]() {
            this->process_jobs();
        });
    }
    
    ~parallel_writer() {
        stop();
    }

    template<typename T>
    void write_all(sigstring_writing& sig_writer, sigstring_writing& csv_writer, sigstring_writing& fastqa_writer, std::vector<T>& data) {
        if (data.empty()) return;
        
        auto job_function = [
            sig_ptr=&sig_writer,
            csv_ptr=&csv_writer,
            fastq_ptr=&fastqa_writer,
            data_vec=std::move(data)
        ]() mutable {
            (*sig_ptr)(data_vec);
            (*csv_ptr)(data_vec);
            (*fastq_ptr)(data_vec);
            std::vector<T>().swap(data_vec);
        };
        
        auto* job = new WriteJob();
        job->func = std::move(job_function);
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = data.size();
        job->queued_time = std::chrono::high_resolution_clock::now();
        
        while (!write_queue.push(job)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    template<typename T>
    void write(sigstring_writing& fastqa_writer, std::vector<T>& data) {
        if (data.empty()) return;
        
        auto job_function = [
            fastq_ptr=&fastqa_writer,
            data_vec=std::move(data)
        ]() mutable {
            (*fastq_ptr)(data_vec);
            std::vector<T>().swap(data_vec);
        };
        
        auto* job = new WriteJob();
        job->func = std::move(job_function);
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = data.size();
        job->queued_time = std::chrono::high_resolution_clock::now();
        
        while (!write_queue.push(job)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
    
    void write_slabs(sigstring_writing& fastqa_writer, std::vector<std::string>&& slabs) {
        size_t total_bytes = 0;
        for (auto const& s : slabs) total_bytes += s.size();
        if (total_bytes == 0) return;
        
        auto* job = new WriteJob();
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = total_bytes;
        job->queued_time = std::chrono::high_resolution_clock::now();
        job->func = [fastq_ptr=&fastqa_writer, slabs_vec=std::move(slabs)]() mutable {
            (*fastq_ptr)(slabs_vec);
            std::vector<std::string>().swap(slabs_vec);
        };
        
        while (!write_queue.push(job)) std::this_thread::yield();
    }
    
    void write_raw_string(sigstring_writing& writer, std::string&& data) {
        if (data.empty()) return;
        
        size_t data_size = data.size();
        
        auto job_function = [
            writer_ptr=&writer,
            data_str=std::move(data)
        ]() mutable {
            (*writer_ptr)(data_str);
            std::string().swap(data_str);
        };
        
        auto* job = new WriteJob();
        job->func = std::move(job_function);
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = data_size;
        job->queued_time = std::chrono::high_resolution_clock::now();
        
        while (!write_queue.push(job)) {
            std::this_thread::yield();
        }
    }

    void stop() {
        if (!running.load()) return;
        
        running.store(false, std::memory_order_release);
        
        if (writer_thread.joinable()) {
            writer_thread.join();
        }
        
        std::cout << "\n[Writer Thread] Statistics:" << std::endl;
        std::cout << "  Total jobs processed: " << jobs_processed.load() << std::endl;
        std::cout << "  Total write time: " << total_write_time_ms / 1000.0 << " seconds" << std::endl;
        if (jobs_processed.load() > 0) {
            std::cout << "  Average write time per chunk: " << total_write_time_ms / jobs_processed.load() << " ms" << std::endl;
            std::cout << "  Average queue time per chunk: " << total_queue_time_ms / jobs_processed.load() << " ms" << std::endl;
        }
    }
    
    double get_total_write_time_ms() const {
        return total_write_time_ms;
    }
    
    size_t get_jobs_processed() const {
        return jobs_processed.load();
    }

private:
    struct WriteJob {
        std::function<void()> func;
        size_t chunk_id;
        size_t size;
        std::chrono::time_point<std::chrono::high_resolution_clock> queued_time;
    };
    
    boost::lockfree::queue<WriteJob*, boost::lockfree::capacity<64>> write_queue;
    std::atomic<bool> running;
    std::atomic<size_t> next_chunk_id{1};
    std::atomic<size_t> jobs_processed{0};
    std::thread writer_thread;
    
    double total_write_time_ms;
    double total_queue_time_ms;
    
    void process_jobs() {
        while (running.load(std::memory_order_acquire)) {
            WriteJob* job = nullptr;
            
            if (write_queue.pop(job)) {
                auto start_time = std::chrono::high_resolution_clock::now();
                auto queue_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    start_time - job->queued_time).count();
                total_queue_time_ms += queue_time_ms;
                
                job->func();
                
                auto end_time = std::chrono::high_resolution_clock::now();
                auto write_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    end_time - start_time).count();
                total_write_time_ms += write_time_ms;
                
                jobs_processed.fetch_add(1, std::memory_order_relaxed);
                
                std::cout << "\n[io_writer] chunk=" << job->chunk_id
                          << " bytes=" << job->size
                          << " queued_ms=" << queue_time_ms
                          << " wrote_ms=" << write_time_ms
                          << std::endl;
                
                delete job;
            } else {
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            }
        }
        
        // Drain remaining jobs
        WriteJob* job = nullptr;
        while (write_queue.pop(job)) {
            job->func();
            delete job;
            jobs_processed.fetch_add(1, std::memory_order_relaxed);
        }
    }
};