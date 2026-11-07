#pragma once
#include "rad_headers.h"
#include <cerrno>
#include <exception>
#include <system_error>

// needed a way for IO to work on .gz and non .gz files, and needed a PIGZ/non-PIGZ solution for files.
// finangled a way to deal with this by just opening two kseq instantiations because you can't trade them off once opened.
// so just open two, kseq_fd and kseq_gz, and use the appropriate one based on whether pigz is being used or not/filetype.
// this is all wrapped up in file_streaming and file_pointer classes below. Learned a lot about reader threads and async IO.
// basically, the smart way to do things is to have pigz do the decompression in parallel with multiple threads (`pigz_reading`) 
// and then have a reader class (`file_streaming`) that manages the file pointers and kseq instances 
// for reading stuff through kseq.

/**
 * @brief kseq intialization for regular file
 */
namespace kseq_fd {
    #ifndef KSEQ_INIT_INT_READ
    #define KSEQ_INIT_INT_READ
    KSEQ_INIT(int, read)
    #endif
}

/**
 * @brief kseq initialization for gzFile
 */
namespace kseq_gz {
    #ifndef KSEQ_INIT_GZFILE_GZREAD  
    #define KSEQ_INIT_GZFILE_GZREAD
    KSEQ_INIT(gzFile, gzread)
    #endif
}

/**
 * @brief Parallel decompression reading using pigz
 * 
 * Manages pigz subprocess for multi-threaded decompression of gzipped files
 */
class pigz_reading {
public:
    static std::string shell_quote_posix(const std::string& s) {
        // POSIX /bin/sh-safe single-quote escaping:
        // abc'd -> 'abc'"'"'d'
        std::string out;
        out.reserve(s.size() + 2);
        out.push_back('\'');
        for (char c : s) {
            if (c == '\'') {
                out.append("'\"'\"'");
            } else {
                out.push_back(c);
            }
        }
        out.push_back('\'');
        return out;
    }

    static std::string sanitize_user_path(std::string p) {
        // Trim ASCII whitespace (incl. CR/LF) commonly introduced by CSV/TSV parsing.
        const auto start = p.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        const auto end = p.find_last_not_of(" \t\r\n");
        p = p.substr(start, end - start + 1);

        // Strip UTF-8 BOM if present.
        if (p.size() >= 3) {
            const unsigned char b0 = static_cast<unsigned char>(p[0]);
            const unsigned char b1 = static_cast<unsigned char>(p[1]);
            const unsigned char b2 = static_cast<unsigned char>(p[2]);
            if (b0 == 0xEF && b1 == 0xBB && b2 == 0xBF) {
                p.erase(0, 3);
            }
        }

        // Strip surrounding quotes if the whole field is quoted.
        if (p.size() >= 2) {
            const char f = p.front();
            const char l = p.back();
            if ((f == '"' && l == '"') || (f == '\'' && l == '\'')) {
                p = p.substr(1, p.size() - 2);
            }
        }

        // Expand leading "~/" to $HOME (shell would normally do this).
        if (!p.empty() && p[0] == '~' && (p.size() == 1 || p[1] == '/')) {
            if (const char* home = std::getenv("HOME"); home && *home) {
                p = std::string(home) + p.substr(1);
            }
        }
        return p;
    }

    /**
     * @brief Open a pigz decompression pipe for reading
     * @param input_path Path to compressed file (.gz, .fq.gz, .fastq.gz, etc.)
     * @param threads Number of decompression threads (default: 4)
     * @return FILE* pointer to pipe, or nullptr on failure
     */
    inline FILE* open_pigz_read_pipe(const std::string& input_path, int threads = 4) {
        // Resolve pigz (PATH-first; honors RAD_PIGZ / RAD_NO_PIGZ). "" => skip.
        const std::string& pigz = resolve_pigz();
        if (pigz.empty()) {
            return nullptr;
        }

        const std::string path = sanitize_user_path(input_path);
        if (path.empty()) {
            return nullptr;
        }

        // Avoid invoking pigz on a path that doesn't exist; this also makes it
        // much easier to diagnose hidden characters (e.g. CR, BOM) in the path.
        boost::system::error_code ec;
        if (!boost::filesystem::is_regular_file(path, ec)) {
            return nullptr;
        }

        // Use pigz -dc for decompression to stdout
        // -d = decompress, -c = write to stdout, -p = threads
        std::string cmd = shell_quote_posix(pigz) + " -dc -p " + std::to_string(threads) +
                          " " + shell_quote_posix(path);
        
        FILE* fp = popen(cmd.c_str(), "r");
        if (!fp) {
            return nullptr;
        }
        
        //removed because if you set a static thread_local buffer, it messes with multiple files opened in different threads.
        // will cause issues streaming from multiple files in parallel
        // had originally set it up trying to debug, TSPMTFO

        // Large read buffer (16 MB) to reduce read() syscalls
        
        //static thread_local std::vector<char> read_buf(16u << 20);
       // setvbuf(fp, read_buf.data(), _IOFBF, read_buf.size());
        
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
     * @brief Resolve the pigz executable to use, preferring one on PATH (so a
     * conda or system pigz "just works"), then the compile-time PIGZ_PATH build
     * hint. Honors RAD_PIGZ (explicit override) and RAD_NO_PIGZ (disable).
     * Returns "" when pigz is unavailable or disabled. Resolved once + cached.
     */
    static const std::string& resolve_pigz() {
        static const std::string cached = []() -> std::string {
            if (const char* no = std::getenv("RAD_NO_PIGZ"); no && *no) return "";
            auto available = [](const std::string& exe) -> bool {
                if (exe.empty()) return false;
                std::string c = "command -v " + shell_quote_posix(exe) + " > /dev/null 2>&1";
                return std::system(c.c_str()) == 0;
            };
            // explicit override wins (if it resolves)
            if (const char* env = std::getenv("RAD_PIGZ"); env && *env) {
                return available(env) ? std::string(env) : std::string();
            }
            // pigz on PATH -- the common case (conda run-dep, system install)
            if (available("pigz")) return "pigz";
            // fall back to the compile-time build hint if it still resolves
            #ifdef PIGZ_PATH
            { std::string p = PIGZ_PATH; if (p != "pigz" && available(p)) return p; }
            #endif
            return "";
        }();
        return cached;
    }

    /**
     * @brief Check if pigz is available (on PATH, via RAD_PIGZ, or PIGZ_PATH)
     * @return true if a usable pigz was found and is not disabled
     */
    static bool is_pigz_available() {
        return !resolve_pigz().empty();
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
/**
 * @brief File pointer structure for kseq reading
 * @param fp_pigz Pigz pipe handle or regular file handle
 * @param fp_gzip Fallback gzip handle
 * @param seq_fd kseq for FILE* (pigz path or uncompressed)
 * @param seq_gz kseq for gzFile (fallback path)
 * @param fd File descriptor for kseq
 * @param path File path
 * @param using_pigz Flag indicating if pigz is used
 * @param pigz_reader pigz reading manager
 */
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

        /**
         * @brief Constructor to open file pointer
         * @param file_path Path to input file
         * @param pigz_threads Number of pigz threads for decompression
         * @note Automatically detects .gz files and uses pigz/gzip as needed
         */
        file_pointer(const std::string& file_path, int pigz_threads = 4,
                     bool allow_pigz = true)
            : fp_pigz(nullptr), fp_gzip(nullptr), seq_fd(nullptr), seq_gz(nullptr),
              fd(-1), path(file_path), using_pigz(false) {
            path = pigz_reading::sanitize_user_path(path);
            if (path.empty()) {
                return;
            }
            
            // Check if file is compressed
            bool is_compressed = has_gz_extension(path);
            
            if (is_compressed) {
                // Try pigz first for parallel decompression
                if (allow_pigz && pigz_reading::is_pigz_available()) {
                    fp_pigz = pigz_reader.open_pigz_read_pipe(path, pigz_threads);
                    if (fp_pigz) {
                        fd = fileno(fp_pigz);
                        seq_fd = kseq_fd::kseq_init(fd);
                        using_pigz = true;
                        return;
                    }
                }

                // Fallback to single-threaded gzip
                fp_gzip = gzopen(path.c_str(), "r");
                if (fp_gzip) {
                    gzbuffer(fp_gzip, 1 << 22);  // 4 MiB buffer
                    seq_gz = kseq_gz::kseq_init(fp_gzip);
                    using_pigz = false;
                }
            } else {
                // Uncompressed file - use regular file descriptor
                fp_pigz = fopen(path.c_str(), "r");
                if (fp_pigz) {
                    fd = fileno(fp_pigz);
                    seq_fd =  kseq_fd::kseq_init(fd);
                    setvbuf(fp_pigz, nullptr, _IOFBF, 1 << 22);  // 4 MiB buffer
                    using_pigz = true;  // Using FILE*, not gzFile
                }
            }
        }

        /**
         * @brief Destructor to clean up resources
         */
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
    boost::filesystem::directory_iterator dir_iter;
    bool dir_status;
    std::string path;
    int pigz_threads;
    bool allow_pigz;

    /**
     * @brief Constructor to initialize file streaming
     * @param input_path Path to input file or directory
     * @param threads Number of pigz threads for decompression
     * @throws std::runtime_error if file cannot be opened or path is invalid
     */
    explicit file_streaming(const std::string& input_path, int threads = 4,
                            bool use_pigz = true)
        : pigz_threads(threads), allow_pigz(use_pigz) {
        path = pigz_reading::sanitize_user_path(input_path);
        dir_status = is_directory(path);
        if (dir_status) {
            dir_iter = boost::filesystem::directory_iterator(path);
            open_next_file();
        } else if (is_fastqa(path)) {
            current_file = std::make_unique<file_pointer>(
                path, pigz_threads, allow_pigz);
            if (!current_file->is_valid()) {
                throw std::runtime_error("Failed to open sequence file: " + path);
            }
        } else {
            throw std::invalid_argument("Invalid input path: " + path);
        }
    }

    bool open_next_file() {
        if (!dir_status) return false;
        while (dir_iter != boost::filesystem::directory_iterator()) {
            const auto& entry = *dir_iter++;
            if (!is_file(get_file_path(entry))) continue;
            
            std::string file_path = entry.path().string();
            if (is_fastqa(file_path)) {
                current_file = std::make_unique<file_pointer>(
                    file_path, pigz_threads, allow_pigz);
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
        return boost::filesystem::is_directory(path);
    }

    static bool is_file(const std::string& path) {
        return boost::filesystem::is_regular_file(path);
    }

    static std::string get_file_path(const boost::filesystem::directory_entry& entry) {
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
    /**
     * @brief structure of an incoming read
     * @param id original read id `string`
     * @param comment optional comment `string`
     * @param seq  `string` sequence
     * @param qual  `string` qual string for fastq
     * @param is_fastq `bool` flag indicating if read is from FASTQ
     */
    struct sequence {
        std::string id;
        std::string comment;
        std::string seq;
        std::string qual;
        bool is_fastq;
    };

    file_streaming& files;
    int read_error_code_ = 0;
    std::string read_error_path_;
    uint64_t empty_sequences_skipped_ = 0;
    
    explicit read_streaming(file_streaming& f) : files(f) {}

    bool read_failed() const {
        return read_error_code_ < -1;
    }

    int read_error_code() const {
        return read_error_code_;
    }

    const std::string& read_error_path() const {
        return read_error_path_;
    }

    uint64_t empty_sequences_skipped() const {
        return empty_sequences_skipped_;
    }

    void note_empty_sequence(file_streaming::file_pointer& fp) {
        ++empty_sequences_skipped_;
        if (empty_sequences_skipped_ <= 10) {
            const char* name = fp.name_s();
            std::cerr << "Warning: skipping zero-length sequence record"
                      << (name ? " '" + std::string(name) + "'" : "")
                      << " in " << files.current_file_path() << "\n";
        } else if (empty_sequences_skipped_ == 11) {
            std::cerr
                << "Warning: additional zero-length sequence warnings "
                   "suppressed\n";
        }
    }

    /**
     * @brief Get a view of the next sequence without copying read metadata.
     *
     * The view remains valid only until the next read from this object.  This
     * is intended for serial consumers that finish processing a sequence
     * before advancing kseq.
     */
    std::optional<std::string_view> next_sequence_view() {
        while (files.current_file || files.open_next_file()) {
            auto& fp = files.current_file;
            int l = fp->kseq_read_any();
            if (l == 0) {
                note_empty_sequence(*fp);
                continue;
            }
            if (l > 0) {
                return std::string_view(
                    fp->seq_s(), static_cast<size_t>(l));
            }
            if (l < -1) {
                read_error_code_ = l;
                read_error_path_ = files.current_file_path();
                return std::nullopt;
            }

            if (files.dir_status) {
                files.current_file.reset();
            } else {
                break;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief get the next sequence from the stream
     * @return sequence object, std::nullopt if end of stream
     */
    std::optional<sequence> next_sequence() {
        while (files.current_file || files.open_next_file()) {
            auto& fp = files.current_file;
            int l = fp->kseq_read_any();
            if (l == 0) {
                note_empty_sequence(*fp);
                continue;
            }
            if (l > 0) {
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
            if (l < -1) {
                read_error_code_ = l;
                read_error_path_ = files.current_file_path();
                return std::nullopt;
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

/**
 * @brief Memory utilities for read streaming, designed for monitoring memory usage of read streaming
 */
namespace readmem_utils {

    /**
     * @brief Get the memory bytes used by a string
     * @param s Input string
     * @return Number of bytes used by the string
     */
    static inline size_t bytes_of(const std::string& s) noexcept { 
        return s.capacity(); 
    }

    /**
     * @brief Get the memory bytes used by a read_streaming::sequence
     * @param r Input read_streaming::sequence
     * @return Number of bytes used by the sequence
     */
    static inline size_t bytes_of_rec(const read_streaming::sequence& r) noexcept {
        return bytes_of(r.id) + bytes_of(r.comment) + bytes_of(r.seq) + bytes_of(r.qual) + sizeof(r);
    }

    /**
     * @brief Get the memory bytes used by a vector of read_streaming::sequence
     * @tparam Vec Vector type
     * @param v Input vector
     * @return Number of bytes used by the vector
     */
    template <class Vec>
    static inline size_t bytes_of_vec(const Vec& v) noexcept {
        using T = typename Vec::value_type;
        size_t total = v.capacity() * sizeof(T);     // vector buffer itself
        for (auto const& e : v) total += bytes_of_rec(e); // payloads (strings' capacities)
        return total;
    }

    /**
     * @brief In-flight memory tracking for read streaming, tracks bytes and chunks in progress with `std::atomic<size_t>` counters
     * @param bytes Number of bytes in flight
     * @param chunks Number of chunks in flight
     * @param peak_bytes Peak bytes in flight
     * @param peak_chunks Peak chunks in flight
     */
    struct data_in_flight {
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

    inline data_in_flight& inflight() { 
        static data_in_flight g; 
        return g; 
    }

    //
    struct scoped_in_flight {
        size_t b{0}; bool active{false};
        explicit scoped_in_flight(size_t bytes) : b(bytes), active(true) { inflight().add(b); }
        ~scoped_in_flight(){ if (active) inflight().sub(b); }
        scoped_in_flight(const scoped_in_flight&) = delete;
        scoped_in_flight& operator=(const scoped_in_flight&) = delete;
    };

    static inline double rss_gib() {
    #if defined(__APPLE__)
        // macOS: use Mach task info
        mach_task_basic_info info;
        mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) != KERN_SUCCESS) {
            return -1.0;
        }
        const double rss_bytes = static_cast<double>(info.resident_size);
        return rss_bytes / (1024.0 * 1024.0 * 1024.0);
    #else
        // Linux: read statm
        FILE* f = std::fopen("/proc/self/statm", "r");
        if (!f) return -1.0;
        long long size_pages=0, res_pages=0;
        if (std::fscanf(f, "%lld %lld", &size_pages, &res_pages) != 2) { std::fclose(f); return -1.0; }
        std::fclose(f);
        const double page_kb = double(::sysconf(_SC_PAGESIZE)) / 1024.0;
        const double rss_kb = res_pages * page_kb;
        return rss_kb / (1024.0 * 1024.0);
    #endif
    }

    // --- tiny logger ---
    static inline void log_chunk(std::size_t chunk_bytes, std::size_t chunk_id=0) {
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
                    int64_t max_seqs = -1,
                    bool verbose = false,
                    size_t log_every_chunks = 0,
                    bool fail_on_read_error = false) {
            file_streaming files(input_path, pigz_threads_);
            read_streaming reader(files);

            std::mutex reader_mutex;
            std::atomic<bool>   more_data{true};
            std::atomic<size_t> seqs_processed{0};
            std::atomic<size_t> chunks_seen{0};
            std::atomic<bool> worker_failed{false};
            std::mutex error_mutex;
            std::exception_ptr worker_error;

            // Exceptions must not leave an OpenMP structured block. Preserve
            // the first failure and report it to the caller after all workers
            // have stopped; this also carries asynchronous output errors.
            #pragma omp parallel num_threads(num_threads)
            {
              try {
                while (!worker_failed.load(std::memory_order_acquire)) {
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

                    if (worker_failed.load(std::memory_order_acquire)) break;
                    const size_t __chunk_bytes = readmem_utils::bytes_of_vec(chunk);
                    const size_t __chunk_id = next_chunk_id_.fetch_add(1, std::memory_order_relaxed);
                    readmem_utils::scoped_in_flight __guard(__chunk_bytes);
                    readmem_utils::log_chunk(__chunk_bytes, __chunk_id);
                    chunk_func(chunk, input_path);

                    if (verbose) {
                        size_t cseen = chunks_seen.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (log_every_chunks == 0 || (cseen % log_every_chunks == 0)) {
                            std::fprintf(stderr, "[chunk_streaming] chunk=%zu size=%zu seqs=%zu\n",
                                         __chunk_id, __chunk_bytes, chunk.size());
                        }
                    }
                }
              } catch (...) {
                worker_failed.store(true, std::memory_order_release);
                more_data.store(false, std::memory_order_release);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
              }
            }
            if (worker_error) std::rethrow_exception(worker_error);
            if (fail_on_read_error && reader.read_failed()) {
                throw std::runtime_error("Failed while reading sequence input: " +
                    reader.read_error_path() + " (parser code " +
                    std::to_string(reader.read_error_code()) + ")");
            }
        }


    private:
    size_t chunk_size_;
    int pigz_threads_;
    int max_in_flight_;  // -1 = unlimited, >0 = limited
    std::atomic<size_t> seqs_processed_;
    std::atomic<int> peak_in_flight_;
    std::atomic<size_t> next_chunk_id_{1};
};

/**
 * @brief pigz writing for parallel compression
 */
class pigz_writing {
public:
    inline FILE* open_pigz_pipe(const std::string& out_path, int threads, int level=1) {
        // Resolve pigz (PATH-first; honors RAD_PIGZ / RAD_NO_PIGZ). "" => skip.
        const std::string& pigz = pigz_reading::resolve_pigz();
        if (pigz.empty()) return nullptr;

        std::string cmd = "'" + pigz + "' -c -p " + std::to_string(threads) +
                        " -" + std::to_string(level) + " > '" + out_path + "'";
        FILE* fp = popen(cmd.c_str(), "w");
        if (!fp){
            return nullptr;
        } 
            // 16MB stdio buffer to reduce write() syscalls

            // removed because in debug mode, the buffer gets shared between all of the thread outputs
            // so then everything gets scrambled
            // TSPMO
        //static std::vector<char> buf(16u << 20);
        //setvbuf(fp, buf.data(), _IOFBF, buf.size());
        return fp;
    }

    inline void pigz_write(FILE* fp, const char* data, size_t n) {
        if (!n) return;
        size_t written = 0;
        while (written < n) {
            errno = 0;
            const size_t count = std::fwrite(data + written, 1, n - written, fp);
            written += count;
            if (written == n) break;
            const int error = errno ? errno : EIO;
            if (error == EINTR) {
                std::clearerr(fp);
                continue;
            }
            throw std::system_error(error, std::generic_category(), "writing compressed output");
        }
    }

    inline int close_pigz_pipe(FILE* fp) {
        if (!fp) return -1;
        std::fflush(fp);
        return pclose(fp);
    }
};

namespace boost_stream = boost::iostreams;
class sigstring_writing {
public:
    enum class format { FASTQA, SIGSTRING, CSV };

private:
    format fmt;
    bool compress_;
    bool use_pigz_ = false;
    bool rc_umi_ = true;
    bool finished_ = false;
    FILE* stream_fp_ = nullptr; // Non-owning; the caller closes its saved stdout.

    pigz_writing pigz_;
    FILE* pigz_fp_ = nullptr;
    std::string pigz_buf_;
    size_t pigz_flush_threshold_ = (8u << 20);

    std::ofstream file_;
    boost_stream::filtering_ostream out_;

    static void flush_stream_(FILE* stream) {
        while (true) {
            errno = 0;
            if (std::fflush(stream) == 0) return;
            const int error = errno ? errno : EIO;
            if (error == EINTR) {
                std::clearerr(stream);
                continue;
            }
            throw std::system_error(error, std::generic_category(), "flushing sequence output");
        }
    }

    static void write_stream_(FILE* stream, const char* data, size_t size) {
        size_t written = 0;
        while (written < size) {
            errno = 0;
            const size_t count = std::fwrite(data + written, 1, size - written, stream);
            written += count;
            if (written == size) break;
            const int error = errno ? errno : EIO;
            if (error == EINTR) {
                std::clearerr(stream);
                continue;
            }
            throw std::system_error(error, std::generic_category(), "writing sequence output");
        }
    }

    void write_bytes_(const std::string& data) {
        if (finished_) throw std::logic_error("writing to a finished output stream");
        if (data.empty()) return;
        if (stream_fp_) {
            write_stream_(stream_fp_, data.data(), data.size());
        } else if (use_pigz_) {
            pigz_buf_.append(data);
            if (pigz_buf_.size() >= pigz_flush_threshold_) {
                pigz_.pigz_write(pigz_fp_, pigz_buf_.data(), pigz_buf_.size());
                pigz_buf_.clear();
            }
        } else {
            out_.write(data.data(), static_cast<std::streamsize>(data.size()));
            if (!out_) throw std::system_error(EIO, std::generic_category(), "writing output file");
        }
    }

    template <typename T> void write_one_(const T& item) {
        switch (fmt) {
            case format::FASTQA:
                write_bytes_(item.to_fastqa(rc_umi_));
                break;
            case format::SIGSTRING: {
                auto text = item.to_sigstring();
                text.push_back('\n');
                write_bytes_(text);
                break;
            }
            case format::CSV:
                write_bytes_(item.to_csv());
                break;
        }
    }

public:
    sigstring_writing(std::string output_path, format output_format, bool compress,
                      bool append, int pigz_threads, int level = 1)
        : fmt(output_format), compress_(compress) {
        const auto mode = std::ios::out | std::ios::binary |
                          (append ? std::ios::app : std::ios::trunc);
        if (compress_) {
            if (output_path.size() < 3 || output_path.substr(output_path.size() - 3) != ".gz")
                output_path += ".gz";
            pigz_fp_ = pigz_.open_pigz_pipe(output_path, pigz_threads, level);
            if (pigz_fp_) {
                use_pigz_ = true;
                return;
            }
        }
        file_.open(output_path, mode);
        if (!file_.is_open()) throw std::runtime_error("Failed to open file: " + output_path);
        if (compress_) {
            boost_stream::gzip_params params;
            params.level = level;
            out_.push(boost_stream::gzip_compressor(params));
        }
        out_.push(file_);
        out_.exceptions(std::ios::badbit | std::ios::failbit);
    }

    // The data stream is separate from std::cout, which may carry diagnostics.
    // Flush complete chunks so small stdout runs are visible promptly.
    sigstring_writing(FILE* stream, format output_format)
        : fmt(output_format), compress_(false), stream_fp_(stream) {
        if (!stream_fp_) throw std::invalid_argument("null sequence output stream");
    }

    void set_rc_umi(bool value) { rc_umi_ = value; }
    void set_pigz_flush_threshold(size_t bytes) { pigz_flush_threshold_ = bytes; }

    template <typename T> void operator()(const std::vector<T>& chunk) {
        if (chunk.empty()) return;
        for (const auto& item : chunk) write_one_(item);
        if (stream_fp_) flush_stream_(stream_fp_);
    }

    void operator()(const std::string& slab) {
        if (slab.empty()) return;
        write_bytes_(slab);
        if (stream_fp_) flush_stream_(stream_fp_);
    }

    void operator()(const std::vector<std::string>& slabs) {
        for (const auto& slab : slabs) (*this)(slab);
    }

    // Call explicitly after the asynchronous writer has joined. Destructors
    // only provide fallback cleanup; successful callers must observe failures.
    void finish() {
        if (finished_) return;
        finished_ = true;
        if (stream_fp_) {
            flush_stream_(stream_fp_);
            return;
        }

        std::exception_ptr failure;
        const auto remember_failure = [&] {
            if (!failure) failure = std::current_exception();
        };
        if (use_pigz_) {
            FILE* pipe = pigz_fp_;
            pigz_fp_ = nullptr;
            use_pigz_ = false;
            try {
                if (!pigz_buf_.empty())
                    pigz_.pigz_write(pipe, pigz_buf_.data(), pigz_buf_.size());
            } catch (...) { remember_failure(); }
            pigz_buf_.clear();
            try { flush_stream_(pipe); }
            catch (...) { remember_failure(); }
            errno = 0;
            const int status = pclose(pipe);
            const int close_error = errno;
            if (status == -1) {
                try {
                    throw std::system_error(close_error ? close_error : EIO,
                                            std::generic_category(), "closing output compressor");
                } catch (...) { remember_failure(); }
            } else if (status != 0) {
                try {
                    throw std::runtime_error("Output compressor failed (wait status " +
                                             std::to_string(status) + ")");
                } catch (...) { remember_failure(); }
            }
        } else {
            try { out_.flush(); }
            catch (...) { remember_failure(); }
            try { out_.reset(); }
            catch (...) { remember_failure(); }
            if (file_.is_open()) {
                try {
                    file_.flush();
                    if (!file_) throw std::system_error(EIO, std::generic_category(), "flushing output file");
                } catch (...) { remember_failure(); }
                try {
                    file_.close();
                    if (file_.fail()) throw std::system_error(EIO, std::generic_category(), "closing output file");
                } catch (...) { remember_failure(); }
            }
        }
        if (failure) std::rethrow_exception(failure);
    }

    ~sigstring_writing() noexcept {
        try { finish(); } catch (...) {}
    }
};

// Completed chunks enter a bounded queue. Output failures are caught in the
// writer thread and rethrown by the producer or explicit stop(), never allowed
// to terminate the worker silently while producers wait on a full queue.
class parallel_writer {
public:
    parallel_writer() {
        writer_thread = std::thread([this] { process_jobs(); });
    }

    ~parallel_writer() noexcept {
        try { stop(); } catch (...) {}
    }

    template <typename T> void write(sigstring_writing& output, std::vector<T>& data) {
        if (data.empty()) return;
        const size_t size = data.size();
        auto job = std::make_unique<WriteJob>();
        job->func = [sink = &output, values = std::move(data)]() mutable {
            (*sink)(values);
        };
        enqueue(std::move(job), size);
    }

    void write_slabs(sigstring_writing& output, std::vector<std::string>&& slabs) {
        size_t size = 0;
        for (const auto& slab : slabs) size += slab.size();
        if (!size) return;
        auto job = std::make_unique<WriteJob>();
        job->func = [sink = &output, values = std::move(slabs)]() mutable {
            (*sink)(values);
        };
        enqueue(std::move(job), size);
    }

    void write_raw_string(sigstring_writing& output, std::string&& data) {
        if (data.empty()) return;
        const size_t size = data.size();
        auto job = std::make_unique<WriteJob>();
        job->func = [sink = &output, value = std::move(data)]() mutable {
            (*sink)(value);
        };
        enqueue(std::move(job), size);
    }

    template <typename T>
    void write_debug(sigstring_writing* sig, sigstring_writing* csv,
                     sigstring_writing* fastqa, std::vector<T>& data) {
        if (data.empty()) return;
        const size_t size = data.size();
        auto job = std::make_unique<WriteJob>();
        job->func = [sig, csv, fastqa, values = std::move(data)]() mutable {
            (*sig)(values);
            (*csv)(values);
            (*fastqa)(values);
        };
        enqueue(std::move(job), size);
    }

    bool failed() const noexcept { return failed_.load(std::memory_order_acquire); }

    void rethrow_if_failed() const {
        if (!failed()) return;
        std::exception_ptr error;
        {
            std::lock_guard<std::mutex> lock(error_mutex);
            error = writer_error;
        }
        if (error) std::rethrow_exception(error);
    }

    void stop() {
        std::lock_guard<std::mutex> lock(stop_mutex);
        running.store(false, std::memory_order_release);
        if (writer_thread.joinable()) writer_thread.join();
        discard_pending();
        rethrow_if_failed();
        if (statistics_printed) return;
        statistics_printed = true;
        std::cout << "\n[Writer Thread] Statistics:" << std::endl;
        std::cout << "  Total jobs processed: " << jobs_processed.load() << std::endl;
        std::cout << "  Total write time: " << total_write_time_ms.load() / 1000.0 << " seconds" << std::endl;
        if (jobs_processed.load()) {
            std::cout << "  Average write time per chunk: "
                      << total_write_time_ms.load() / jobs_processed.load() << " ms" << std::endl;
            std::cout << "  Average queue time per chunk: "
                      << total_queue_time_ms.load() / jobs_processed.load() << " ms" << std::endl;
        }
    }

    double get_total_write_time_ms() const { return total_write_time_ms.load(); }
    size_t get_jobs_processed() const { return jobs_processed.load(); }

private:
    struct WriteJob {
        std::function<void()> func;
        size_t chunk_id = 0;
        size_t size = 0;
        std::chrono::time_point<std::chrono::high_resolution_clock> queued_time;
    };
    struct ProducerGuard {
        std::atomic<size_t>& active;
        explicit ProducerGuard(std::atomic<size_t>& value) : active(value) {
            active.fetch_add(1, std::memory_order_acq_rel);
        }
        ~ProducerGuard() { active.fetch_sub(1, std::memory_order_release); }
    };

    boost::lockfree::queue<WriteJob*, boost::lockfree::capacity<64>> write_queue;
    std::atomic<bool> running{true};
    std::atomic<bool> failed_{false};
    std::atomic<size_t> active_producers{0};
    std::atomic<size_t> next_chunk_id{1};
    std::atomic<size_t> jobs_processed{0};
    std::thread writer_thread;
    mutable std::mutex error_mutex;
    std::exception_ptr writer_error;
    std::mutex stop_mutex;
    bool statistics_printed = false;
    std::atomic<double> total_write_time_ms{0.0};
    std::atomic<double> total_queue_time_ms{0.0};

    void enqueue(std::unique_ptr<WriteJob> job, size_t size) {
        ProducerGuard producer(active_producers);
        rethrow_if_failed();
        if (!running.load(std::memory_order_acquire))
            throw std::logic_error("enqueueing output after writer shutdown");
        job->chunk_id = next_chunk_id.fetch_add(1, std::memory_order_relaxed);
        job->size = size;
        job->queued_time = std::chrono::high_resolution_clock::now();
        while (!write_queue.push(job.get())) {
            rethrow_if_failed();
            if (!running.load(std::memory_order_acquire))
                throw std::logic_error("enqueueing output after writer shutdown");
            std::this_thread::yield();
        }
        job.release(); // Queue owns it, including if failure races this enqueue.
        rethrow_if_failed();
    }

    void discard_pending() noexcept {
        WriteJob* job = nullptr;
        while (write_queue.pop(job)) delete job;
    }

    void record_failure() noexcept {
        {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!writer_error) writer_error = std::current_exception();
        }
        failed_.store(true, std::memory_order_release);
        running.store(false, std::memory_order_release);
    }

    void process_jobs() noexcept {
        try {
            while (true) {
                if (failed_.load(std::memory_order_acquire)) break;
                WriteJob* raw_job = nullptr;
                if (write_queue.pop(raw_job)) {
                    std::unique_ptr<WriteJob> job(raw_job);
                    const auto start = std::chrono::high_resolution_clock::now();
                    const auto queued_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        start - job->queued_time).count();
                    total_queue_time_ms.store(total_queue_time_ms.load(std::memory_order_relaxed) + queued_ms,
                                              std::memory_order_relaxed);
                    job->func();
                    const auto wrote_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start).count();
                    total_write_time_ms.store(total_write_time_ms.load(std::memory_order_relaxed) + wrote_ms,
                                              std::memory_order_relaxed);
                    jobs_processed.fetch_add(1, std::memory_order_relaxed);
                    std::cout << "\n[io_writer] chunk=" << job->chunk_id
                              << " size_mib=" << std::fixed << std::setprecision(2)
                              << double(job->size) / (1024.0 * 1024.0)
                              << " queued_ms=" << queued_ms << " wrote_ms=" << wrote_ms << std::endl;
                } else {
                    if (!running.load(std::memory_order_acquire) &&
                        active_producers.load(std::memory_order_acquire) == 0) break;
                    std::this_thread::sleep_for(std::chrono::microseconds(100));
                }
            }
        } catch (...) {
            record_failure();
        }
        // A producer may have observed healthy state just before the failure.
        // Wait for those enqueues and free their buffers without further writes.
        if (failed_.load(std::memory_order_acquire)) {
            while (active_producers.load(std::memory_order_acquire)) {
                discard_pending();
                std::this_thread::yield();
            }
            discard_pending();
        }
    }
};
