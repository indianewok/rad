#pragma once
#include "rad_headers.h"

// Initialize kseq for gzipped file reading
KSEQ_INIT(gzFile, gzread)

class file_streaming {
public:

    struct file_pointer {
        gzFile fp;
        kseq_t* seq;
        std::string path;
        
        file_pointer(const std::string& file_path) : path(file_path) {
            fp = gzopen(path.c_str(), "r");
            if (fp) seq = kseq_init(fp);
            else seq = nullptr;
        }
        
        ~file_pointer() {
            if (seq) kseq_destroy(seq);
            if (fp) gzclose(fp);
        }
    };

    std::unique_ptr<file_pointer> current_file;
    std::filesystem::directory_iterator dir_iter;
    bool dir_status;
    std::string path;

    file_streaming(const std::string& input_path) : path(input_path) {
        dir_status = is_directory(path);
        if (dir_status) {
            dir_iter = std::filesystem::directory_iterator(path);
            open_next_file();
        } else if (is_fastqa(path)) {
            current_file = std::make_unique<file_pointer>(path);
            if (!current_file->fp || !current_file->seq) {
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
                current_file = std::make_unique<file_pointer>(file_path);
                if (current_file->fp && current_file->seq) {
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

class read_streaming {
public:
    //structure of a sequence that gets streamed in 
    struct sequence {
        std::string id; // read id (whatever comes out of the sequencer)
        std::string comment; // any comments that come in
        std::string seq; // actual sequence
        std::string qual; // qual score 
        bool is_fastq; // is this a fastq read
    };

    file_streaming& files;
    read_streaming(file_streaming& f) : files(f) {}
    std::optional<sequence> next_sequence() {
        while (files.current_file || files.open_next_file()) {
            int l = kseq_read(files.current_file->seq);
            if (l >= 0) {
                sequence seq;
                seq.id = files.current_file->seq->name.s;
                seq.seq = files.current_file->seq->seq.s;
                seq.is_fastq = files.current_file->seq->qual.l > 0;
                if (files.current_file->seq->comment.l) {
                    seq.comment = files.current_file->seq->comment.s;
                }
                if (seq.is_fastq) {
                    seq.qual = files.current_file->seq->qual.s;
                }
                return seq;
            } else if (files.dir_status) {
                files.current_file.reset();
            } else {
                break;
            }
        }
        return std::nullopt;
    }
};

template<typename T, typename Function> 
class chunk_streaming {
public:
    // Construct with the max number of sequences per chunk
    explicit chunk_streaming(size_t chunk_size) : chunk_size(chunk_size), seqs_processed(0)
    {}

    // Change the chunk size on the fly
    void modify_chunk_size(size_t new_size) {
        chunk_size = new_size;
    }

    /// Query the current chunk size
    size_t get_chunk_size() const {
        return chunk_size;
    }

    /**
     *  Read sequences from `input_path` and for each chunk
     *  of up to chunk_size call the chunk_func on one of
     *  an OpenMP thread pool.
     *
     *  @param input_path   FASTQ or similar to stream from
     *  @param chunk_func   Fn(vector<T> const&, string const&)
     *  @param num_threads  # of OpenMP worker threads
     *  @param max_seqs     stop after this many total reads (or <0 for no limit)
     */
    void process_chunks(const std::string& input_path, Function chunk_func, int num_threads, int64_t max_seqs = -1) {
        // Open the file/stream handles
        file_streaming files(input_path);
        read_streaming reader(files);
        // pre-allocate one chunk buffer
        std::vector<T> chunk;
        chunk.reserve(chunk_size);
        // Launch OpenMP parallel region
        #pragma omp parallel num_threads(num_threads)
        #pragma omp single
        {
            while (true) {
                chunk.clear();
                // Fill up to chunk_size or until max_seqs reached
                while (chunk.size() < chunk_size
                       && (max_seqs < 0 || seqs_processed < (size_t)max_seqs))
                {
                    auto seq = reader.next_sequence();
                    if (!seq) break;
                    chunk.push_back(std::move(*seq));
                    ++seqs_processed;
                }
                if (chunk.empty()) {
                    // no more data -> exit loop
                    break;
                }
                // Make a copy for the task
                auto this_chunk = std::move(chunk);
                // Dispatch an OpenMP task to handle it
                #pragma omp task firstprivate(this_chunk)
                {
                    chunk_func(this_chunk, input_path);
                }
            }
            // Wait for all outstanding tasks to complete
            #pragma omp taskwait
        }
    }

private:
    size_t chunk_size; //< max reads per chunk
    std::atomic<size_t> seqs_processed; //< global counter
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
        format fmt;
        std::ofstream file_;
        bio::filtering_ostream out_; // the chain: [gzip_compressor?] -> file_
    
    public:
        /// @param output_path   base filename ("myout"), or with .gz
        /// @param output_format which of the three formats to emit
        /// @param compress      if true, writes gzipped output; otherwise plaintext
        /// @param append        if true, opens in append mode; otherwise truncates
        /// @param level         gzip compression level [1..9] (1=fastest, 9=best)
        sigstring_writing(std::string output_path,format output_format, bool compress, bool append,int level = 1): fmt(output_format) {
            std::ios_base::openmode mode = 
                std::ios::out | std::ios::binary | (append ? std::ios::app : std::ios::trunc);
    
            if (compress) {
                // ensure .gz suffix
                if (output_path.size() < 3 ||
                    output_path.substr(output_path.size()-3) != ".gz")
                    output_path += ".gz";
    
                file_.open(output_path, mode);
                if (!file_.is_open())
                    throw std::runtime_error("Failed to open file: " + output_path);
    
                bio::gzip_params params;
                params.level = level;        // set the compression level
                out_.push(bio::gzip_compressor(params));
                out_.push(file_);
            }
            else {
                file_.open(output_path, mode);
                if (!file_.is_open())
                    throw std::runtime_error("Failed to open file: " + output_path);
                // no compressor ⇒ write straight to file_
                out_.push(file_);
            }
        }
    
        template<typename T> void operator()(const std::vector<T>& chunk) {
                for (auto const& item : chunk) {
                    switch (fmt) {
                        case format::FASTQA:
                            out_ << item.to_fastqa();
                            break;
                        case format::SIGSTRING:
                            out_ << item.to_sigstring() << "\n";
                            break;
                        case format::CSV:
                            out_ << item.to_csv();
                            break;
                    }
            }
        }
        ~sigstring_writing() {
            // flushing:
            out_.reset();    // first flush and pop all filters
            file_.flush();   // then flush the underlying file
        }
    };

class parallel_writer {
    public:
        parallel_writer() : running(true), queue_size(0), total_write_time_ms(0), total_queue_time_ms(0) {
            writer_thread = std::thread([this]() {
                this->process_jobs();
            });
        }
        
        ~parallel_writer() {
            stop();
        }
        
        template<typename T>
        void write(sigstring_writing& sig_writer, sigstring_writing& csv_writer, sigstring_writing& fastqa_writer, std::vector<T>& data) {
            // Check if empty
            if (data.empty()) return;
            
            // Create a type-erased job using std::function
            auto job_function = [
                sig_ptr=&sig_writer, 
                csv_ptr=&csv_writer, 
                fastq_ptr=&fastqa_writer, 
                data_vec=std::move(data)  // Move the data
            ]() mutable {
                (*sig_ptr)(data_vec);
                (*csv_ptr)(data_vec);
                (*fastq_ptr)(data_vec);
                
                // Clear data to free memory
                std::vector<T>().swap(data_vec);
            };
            
            // Create a job
            auto job = std::make_shared<WriteJob>();
            job->func = std::move(job_function);
            job->chunk_id = next_chunk_id++;
            job->size = data.size();
            
            // Control queue size to limit memory usage
            constexpr size_t MAX_QUEUE_SIZE = 10;
            
            // Add job to the queue
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                
                // Wait if queue gets too large (backpressure)
                queue_cv.wait(lock, [this, MAX_QUEUE_SIZE]() {
                    return queue_size < MAX_QUEUE_SIZE || !running;
                });
                
                if (!running) return;
                
                // Add job to queue
                job->queued_time = std::chrono::high_resolution_clock::now();
                write_queue.push_back(job);
                queue_size++;
            }
            
            queue_cv.notify_one();
        }
        
        void stop() {
            if (!running) return;
            
            {
                std::lock_guard<std::mutex> lock(queue_mutex);
                running = false;
            }
            queue_cv.notify_all();
            
            if (writer_thread.joinable()) {
                writer_thread.join();
            }
            
            // Print write statistics
            std::cout << "\n[Writer Thread] Statistics:" << std::endl;
            std::cout << "  Total jobs processed: " << jobs_processed << std::endl;
            std::cout << "  Total write time: " << total_write_time_ms / 1000.0 << " seconds" << std::endl;
            if (jobs_processed > 0) {
                std::cout << "  Average write time per chunk: " << total_write_time_ms / jobs_processed << " ms" << std::endl;
                std::cout << "  Average queue time per chunk: " << total_queue_time_ms / jobs_processed << " ms" << std::endl;
            }
        }
        
        // Get total write time in milliseconds
        double get_total_write_time_ms() const {
            return total_write_time_ms;
        }
        
        // Get number of jobs processed
        size_t get_jobs_processed() const {
            return jobs_processed;
        }
        
    private:
        struct WriteJob {
            std::function<void()> func; //  function that does the writing
            size_t chunk_id;
            size_t size;
            std::chrono::time_point<std::chrono::high_resolution_clock> queued_time;
        };
        
        std::vector<std::shared_ptr<WriteJob>> write_queue;
        std::mutex queue_mutex;
        std::condition_variable queue_cv;
        std::thread writer_thread;
        std::atomic<bool> running;
        std::atomic<size_t> queue_size;
        std::atomic<size_t> next_chunk_id{1};
        std::atomic<size_t> jobs_processed{0};
        
        // Use regular doubles with mutex for thread safety
        double total_write_time_ms;
        double total_queue_time_ms;
        std::mutex timing_mutex;
        
        // For logging write times
        std::mutex log_mutex;
        
        void process_jobs() {
            while (running || queue_size > 0) {
                std::shared_ptr<WriteJob> job;
                
                // Get a job from the queue
                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    
                    if (write_queue.empty()) {
                        if (!running) break;
                        
                        // Wait for more jobs
                        queue_cv.wait(lock, [this]() {
                            return !write_queue.empty() || !running;
                        });
                        
                        if (write_queue.empty() && !running) break;
                    }
                    
                    if (!write_queue.empty()) {
                        job = write_queue.front();
                        write_queue.erase(write_queue.begin());
                        queue_size--;
                        
                        // Notify that queue has space
                        lock.unlock();
                        queue_cv.notify_one();
                    }
                }
                
                // Process the job
                if (job) {
                    // Calculate queue time
                    auto start_time = std::chrono::high_resolution_clock::now();
                    auto queue_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        start_time - job->queued_time).count();
                    // Update queue time with mutex protection
                    {
                        std::lock_guard<std::mutex> lock(timing_mutex);
                        total_queue_time_ms += queue_time_ms;
                    }
                    
                    // Execute the job function (which does the writing)
                    job->func();
                    
                    // Calculate write time
                    auto end_time = std::chrono::high_resolution_clock::now();
                    auto write_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time).count();
                    
                    // Update write time with mutex protection
                    {
                        std::lock_guard<std::mutex> lock(timing_mutex);
                        total_write_time_ms += write_time_ms;
                    }
                    
                    jobs_processed++;
                    
                    // Log the write time
                    {
                        std::lock_guard<std::mutex> lock(log_mutex);
                        std::cout << "[io_writer] chunk " << job->chunk_id << " written in " << write_time_ms << " ms" << std::endl;
                    }
                    // Release the job (shared_ptr will clean up the memory)
                    job.reset();
                }
            }
        }
    };

namespace io_bc_split_utils {
    // Extract barcode from FASTQ header line
    std::string extract_header_bcs(const std::string& header) {
        size_t cb_pos = header.find("CB:Z:");
        if (cb_pos == std::string::npos) return "";
        
        size_t start = cb_pos + 5;
        size_t end = header.find(' ', start);
        if (end == std::string::npos) end = header.length();
        
        return header.substr(start, end - start);
    }

    // Fix messy bcs for filename safety
    std::string fix_messy_header_bcs(const std::string& barcode) {
        std::string sanitized = barcode;
        std::replace(sanitized.begin(), sanitized.end(), '/', '_');
        std::replace(sanitized.begin(), sanitized.end(), '\\', '_');
        std::replace(sanitized.begin(), sanitized.end(), ':', '_');
        std::replace(sanitized.begin(), sanitized.end(), '*', '_');
        std::replace(sanitized.begin(), sanitized.end(), '?', '_');
        std::replace(sanitized.begin(), sanitized.end(), '"', '_');
        std::replace(sanitized.begin(), sanitized.end(), '<', '_');
        std::replace(sanitized.begin(), sanitized.end(), '>', '_');
        std::replace(sanitized.begin(), sanitized.end(), '|', '_');
        
        if (sanitized.empty() || std::all_of(sanitized.begin(), sanitized.end(), ::isspace)) {
            return "unknown";
        }
        return sanitized;
    }

    /**
     * @brief Load whitelist barcodes from CSV file (second column after header)
     */
    std::unordered_set<std::string> load_whitelist_barcodes(const std::string& whitelist_path, bool verbose = false) {
        std::unordered_set<std::string> whitelist;
        
        if (whitelist_path.empty() || !boost::filesystem::exists(whitelist_path)) {
            if (verbose) {
                std::cout << "[split_fastqas] No whitelist file found, processing all barcodes" << std::endl;
            }
            return whitelist; // Empty set means no filtering
        }
        
        if (verbose) {
            std::cout << "[split_fastqas] Loading whitelist from: " << whitelist_path << std::endl;
        }
        
        std::ifstream file(whitelist_path);
        std::string line;
        bool first_line = true;
        
        while (std::getline(file, line)) {
            if (first_line) {
                first_line = false;
                continue; // Skip header
            }
            
            if (line.empty()) continue;
            
            // Get second column (barcode) - skip first column
            size_t first_comma = line.find(',');
            if (first_comma == std::string::npos) continue; // No comma found, skip line
            
            size_t second_comma = line.find(',', first_comma + 1);
            std::string barcode;
            
            if (second_comma != std::string::npos) {
                // Extract between first and second comma
                barcode = line.substr(first_comma + 1, second_comma - first_comma - 1);
            } else {
                // Extract from first comma to end of line
                barcode = line.substr(first_comma + 1);
            }
            
            // Trim whitespace
            barcode.erase(0, barcode.find_first_not_of(" \t\r\n"));
            barcode.erase(barcode.find_last_not_of(" \t\r\n") + 1);
            
            if (!barcode.empty()) {
                whitelist.insert(barcode);
                // Also add reverse complement to whitelist for matching
                whitelist.insert(seq_utils::revcomp(barcode));
            }
        }
        
        if (verbose) {
            std::cout << "[split_fastqas] Loaded " << whitelist.size() / 2 << " barcodes (with reverse complements)" << std::endl;
        }
        
        return whitelist;
    }

    /**
     * @brief Enhanced file writer that manages multiple barcode files efficiently
     */
    class parallel_barcode_writer {
    private:
        std::string output_dir;
        std::unordered_map<std::string, std::unique_ptr<boost::iostreams::filtering_ostream>> writers;
        std::unordered_map<std::string, std::unique_ptr<std::ofstream>> base_streams;
        std::unordered_map<std::string, size_t> access_counter; // Track access for LRU
        mutable std::mutex writers_mutex;
        size_t max_open_files;
        bool verbose;
        int compression_level;
        size_t access_tick;
        
    public:
        parallel_barcode_writer(const std::string& output_directory, bool verbose_mode = false, 
                               size_t max_files = 500, int comp_level = 1) 
            : output_dir(output_directory), max_open_files(max_files), verbose(verbose_mode), 
              compression_level(comp_level), access_tick(0) {
            boost::filesystem::path output_path(output_dir);
            if (!boost::filesystem::exists(output_path)) {
                boost::filesystem::create_directories(output_path);
            }
            
            if (verbose) {
                //std::cout << "[split_fastqas] Using compression level: " << compression_level  << " (1=fastest, 6=default, 9=best)" << std::endl;
            }
        }
        
        void write_records(const std::unordered_map<std::string, std::vector<std::string>>& barcode_records) {
            std::lock_guard<std::mutex> lock(writers_mutex);
            
            // Close some files if we have too many open
            if (writers.size() >= max_open_files) {
                close_lru_files(max_open_files / 3); // Close 1/3 of files
            }
            
            for (const auto& [barcode, records] : barcode_records) {
                // Skip empty record sets
                if (records.empty()){
                    continue;
                }
                // Get or create writer for this barcode
                if (writers.find(barcode) == writers.end()) {
                    create_writer(barcode);
                }
                
                // Update access time for LRU
                access_counter[barcode] = ++access_tick;
                
                // Write all records for this barcode
                auto& writer = *writers[barcode];
                for (const auto& record : records) {
                    writer << record;
                }
                writer.flush();
            }
        }
        
        // Explicit cleanup method
        void close_all() {
            std::lock_guard<std::mutex> lock(writers_mutex);
            
            // Safely close all writers
            for (auto& [barcode, writer] : writers) {
                try {
                    if (writer) {
                        writer->reset();
                    }
                } catch (const std::exception& e) {
                    if (verbose) {
                        std::cerr << "[split_fastqas] Warning: Error closing writer for " 
                                  << barcode << ": " << e.what() << std::endl;
                    }
                }
            }
            
            writers.clear();
            base_streams.clear();
            access_counter.clear();
        }
        
        ~parallel_barcode_writer() {
            close_all();
        }
        
    private:
        void create_writer(const std::string& barcode) {
            try {
                std::string output_filename = (boost::filesystem::path(output_dir) / (barcode + ".fastq.gz")).string();
                
                auto base_stream = std::make_unique<std::ofstream>(output_filename, std::ios::binary | std::ios::app);
                if (!base_stream->is_open()) {
                    throw std::runtime_error("Failed to create output file: " + output_filename);
                }
                
                auto filtering_stream = std::make_unique<boost::iostreams::filtering_ostream>();
                
                // Configure compression parameters
                boost::iostreams::gzip_params gzip_params;
                gzip_params.level = compression_level;
                
                filtering_stream->push(boost::iostreams::gzip_compressor(gzip_params));
                filtering_stream->push(*base_stream);
                
                base_streams[barcode] = std::move(base_stream);
                writers[barcode] = std::move(filtering_stream);
                access_counter[barcode] = ++access_tick;
                
                if (verbose && writers.size() % 100 == 0) {
                    //std::cout << "[split_fastqas] Created writer for barcode: " << barcode << " (total: " << writers.size() << ")" << std::endl;
                }
            } catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "[split_fastqas] Error creating writer for " << barcode 
                              << ": " << e.what() << std::endl;
                }
                throw;
            }
        }
        
        void close_lru_files(size_t num_to_close) {
            if (writers.empty() || num_to_close == 0) return;
            
            // Create vector of (barcode, access_time) pairs
            std::vector<std::pair<std::string, size_t>> access_times;
            access_times.reserve(writers.size());
            
            for (const auto& [barcode, writer] : writers) {
                access_times.emplace_back(barcode, access_counter[barcode]);
            }
            
            // Sort by access time (oldest first)
            std::sort(access_times.begin(), access_times.end(),
                     [](const auto& a, const auto& b) { return a.second < b.second; });
            
            // Close the oldest files
            size_t closed = 0;
            for (const auto& [barcode, access_time] : access_times) {
                if (closed >= num_to_close) break;
                
                try {
                    auto writer_it = writers.find(barcode);
                    if (writer_it != writers.end() && writer_it->second) {
                        writer_it->second->reset();
                        writers.erase(writer_it);
                    }
                    
                    base_streams.erase(barcode);
                    access_counter.erase(barcode);
                    closed++;
                    
                    if (verbose && closed % 50 == 0) {
                        //std::cout << "[split_fastqas] Closed " << closed << " LRU files" << std::endl;
                    }
                } catch (const std::exception& e) {
                    if (verbose) {
                        std::cerr << "[split_fastqas] Warning: Error closing LRU file " 
                                  << barcode << ": " << e.what() << std::endl;
                    }
                }
            }
        }
    };

    /**
     * @brief Parallel barcode splitter using kseq with chunk_streaming
     */
    size_t split_fastqa_file(const std::string& input_fastq_path, 
                             const std::string& output_dir_path,
                             int num_threads,
                             const std::string& whitelist_path = "",
                             bool verbose = false,
                             size_t chunk_size = 50000,
                             int compression_level = 1) {
        
        if (verbose) {
            std::cout << "[split_fastqas] Starting parallel barcode splitting for: " << input_fastq_path << std::endl;
            std::cout << "[split_fastqas] Using " << num_threads << " threads with chunk size " << chunk_size << std::endl;
        }
        
        // Load whitelist
        auto whitelist = load_whitelist_barcodes(whitelist_path, verbose);
        bool use_whitelist = !whitelist.empty();
        
        // Create parallel barcode writer with specified compression level
        parallel_barcode_writer bc_writer(output_dir_path, verbose, 500, compression_level);
        
        // Shared data structures (thread-safe)
        std::unordered_set<std::string> seen_barcodes;
        std::mutex seen_barcodes_mutex;
        std::atomic<size_t> total_reads{0};
        std::atomic<size_t> reads_with_barcodes{0};
        std::atomic<size_t> filtered_reads{0};
        
        // Chunk processing function
        using ChunkFunc = std::function<void(const std::vector<read_streaming::sequence>&, const std::string&)>;
        chunk_streaming<read_streaming::sequence, ChunkFunc> streamer(chunk_size);
        
        ChunkFunc process_chunk = [&](const std::vector<read_streaming::sequence>& chunk, const std::string& /*file*/) {
            // Group records by target barcode within this chunk
            std::unordered_map<std::string, std::vector<std::string>> chunk_records;
            
            for (const auto& seq : chunk) {
                total_reads++;
                
                // Build header string
                std::string header = "@" + seq.id;
                if (!seq.comment.empty()) {
                    header += " " + seq.comment;
                }
                
                std::string raw_barcode = extract_header_bcs(header);
                if (raw_barcode.empty()) {
                    continue; // Skip reads without barcodes
                }
                
                raw_barcode = fix_messy_header_bcs(raw_barcode);
                reads_with_barcodes++;
                
                // Apply whitelist filter if enabled
                if (use_whitelist && whitelist.find(raw_barcode) == whitelist.end()) {
                    filtered_reads++;
                    continue; // Skip reads not in whitelist
                }
                
                // Determine target barcode and suffix
                std::string target_barcode = raw_barcode;
                std::string read_suffix = "";
                
                {
                    std::lock_guard<std::mutex> lock(seen_barcodes_mutex);
                    std::string revcomp_barcode = seq_utils::revcomp(raw_barcode);
                    
                    if (seen_barcodes.count(raw_barcode)) {
                        target_barcode = raw_barcode;
                    } else if (seen_barcodes.count(revcomp_barcode)) {
                        target_barcode = revcomp_barcode;
                        read_suffix = "_rc";
                    } else {
                        // First time seeing this barcode
                        seen_barcodes.insert(raw_barcode);
                        target_barcode = raw_barcode;
                    }
                }
                
                // Build FASTQ record string
                std::string fastq_record = "@" + seq.id + read_suffix;
                if (!seq.comment.empty()) {
                    fastq_record += " " + seq.comment;
                }
                fastq_record += "\n" + seq.seq + "\n+\n";
                fastq_record += (seq.is_fastq ? seq.qual : std::string(seq.seq.length(), 'I')) + "\n";
                
                chunk_records[target_barcode].push_back(fastq_record);
            }
            
            // Write all records for this chunk
            if (!chunk_records.empty()) {
                bc_writer.write_records(chunk_records);
            }
            
            // Progress reporting
            if (verbose && total_reads % 500000 == 0) {
                std::cout << "[split_fastqas] Processed " << total_reads << " reads, "
                          << seen_barcodes.size() << " unique barcodes" << std::endl;
            }
        };
        
        // Process the file in parallel chunks using your existing infrastructure
        streamer.process_chunks(input_fastq_path, process_chunk, num_threads, -1);
        
        // Explicit cleanup of writer
        bc_writer.close_all();
        
        if (verbose) {
            std::cout << "[split_fastqas] Complete! Processed " << total_reads << " total reads" << std::endl;
            std::cout << "[split_fastqas] Reads with barcodes: " << reads_with_barcodes << std::endl;
            if (use_whitelist) {
                std::cout << "[split_fastqas] Reads filtered by whitelist: " << filtered_reads << std::endl;
                std::cout << "[split_fastqas] Reads written: " << (reads_with_barcodes - filtered_reads) << std::endl;
            }
            std::cout << "[split_fastqas] Unique barcode files created: " << seen_barcodes.size() << std::endl;
        }
        
        return seen_barcodes.size();
    }

    /**
     * @brief Enhanced wrapper function with whitelist support
     */
    void split_fastqas(const std::string& output_base, bool verbose, int num_threads, int compression_level = 1) {
        std::string fastq_file = output_base + ".fq.gz";
        std::string split_dir = output_base + "_split_bcs";
        std::string whitelist_file = output_base + "_whitelist.csv";
        
        if (!boost::filesystem::exists(fastq_file)) {
            fastq_file = output_base + ".fq";
            if (!boost::filesystem::exists(fastq_file)) {
                throw std::runtime_error("Output FASTQ file not found: " + output_base + ".fq[.gz]");
            }
        }
        
        // Check for whitelist file
        if (!boost::filesystem::exists(whitelist_file)) {
            if (verbose) {
                std::cout << "[split_fastqas] No whitelist file found at: " << whitelist_file << std::endl;
                std::cout << "[split_fastqas] Processing all barcodes without filtering" << std::endl;
            }
            whitelist_file = ""; // Empty means no filtering
        }
        auto start_time = std::chrono::steady_clock::now();
        size_t num_barcodes = split_fastqa_file(fastq_file, split_dir, num_threads, whitelist_file, verbose);
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start_time);
        std::cout << "[split_fastqas] Split into " << num_barcodes 
                  << " barcode files in " << duration.count() << " seconds" << std::endl;
    }
}