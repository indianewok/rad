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
/*
template<typename T, typename function_to_run>
class chunk_streaming {
public:
    chunk_streaming(size_t chunk_size = 50000) 
        : chunk_size(chunk_size), seqs_processed(0) {}

    void process_chunks(const std::string& input_path, function_to_run chunk_func, int num_threads, int64_t max_seqs) {
        std::vector<std::vector<T>> chunks;
        {
            std::vector<T> chunk;
            chunk.reserve(chunk_size);
            file_streaming files(input_path);
            read_streaming reader(files);
            
            while (auto item = reader.next_sequence()) {
                if (max_seqs > 0 && seqs_processed >= max_seqs)
                    break;
                chunk.push_back(*item);
                seqs_processed++;  
                if (chunk.size() >= chunk_size) {
                    chunks.push_back(std::move(chunk));
                    chunk.clear();
                    chunk.reserve(chunk_size);
                }
            }
            if (!chunk.empty()) {
                chunks.push_back(std::move(chunk));
            }
        }
        
        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp for schedule(dynamic) nowait
            for (size_t i = 0; i < chunks.size(); ++i) {
                std::string current_file_path = "";
                chunk_func(chunks[i], current_file_path);
            }
        }
    }

    int64_t total_sequences() const { return seqs_processed; }

private:
    size_t chunk_size;
    int64_t seqs_processed;
};*/

//this version is a bit complicated but based off of what i understand of minimap2 (i don't understand minimap2, it's pearls before swine)
template<typename T, typename function_to_run>
class chunk_streaming {
public:
    chunk_streaming(size_t chunk_size) 
        : chunk_size(chunk_size), seqs_processed(0) {}
    
    void modify_chunk_size(size_t new_size) {
        chunk_size = new_size;
    }
    
    size_t get_chunk_size() const {
        return chunk_size;
    }

    void process_chunks(const std::string& input_path, function_to_run chunk_func, int num_threads, int64_t max_seqs = -1) {
        const int pipeline_depth = 3; // Number of batches to keep in the pipeline
        
        // Shared queues for pipeline stages
        std::queue<std::vector<T>> read_queue;    // Chunks ready for processing
        std::queue<std::vector<T>> result_queue;  // Processed chunks ready for output
        
        // Synchronization
        std::mutex read_mutex, result_mutex;
        std::condition_variable read_cv, result_cv;
        std::atomic<bool> done_reading{false};
        std::atomic<bool> done_processing{false};
        std::atomic<size_t> active_processing_threads{0};
        
        // Read thread - reads chunks from file
        std::thread reader_thread([&]() {
            file_streaming files(input_path);
            read_streaming reader(files);
            
            while (true) {
                // Create a new chunk
                std::vector<T> chunk;
                chunk.reserve(chunk_size);
                
                // Fill the chunk with sequences
                size_t local_count = 0;
                while (local_count < chunk_size) {
                    if (max_seqs > 0 && seqs_processed >= max_seqs) break;
                    
                    auto seq = reader.next_sequence();
                    if (!seq) break;
                    
                    chunk.push_back(*seq);
                    local_count++;
                    seqs_processed++;
                }
                
                if (chunk.empty()) break; // No more sequences
                
                // Add the chunk to the read queue
                {
                    std::unique_lock<std::mutex> lock(read_mutex);
                    // Wait if queue is too full (memory control)
                    read_cv.wait(lock, [&]() { 
                        return read_queue.size() < pipeline_depth;
                    });
                    
                    read_queue.push(std::move(chunk));
                    lock.unlock();
                    read_cv.notify_one();
                }
            }
            
            done_reading = true;
            read_cv.notify_all();
        });
        
        // Worker threads - process chunks
        std::vector<std::thread> worker_threads;
        for (int i = 0; i < num_threads; i++) {
            worker_threads.emplace_back([&, i]() {
                while (true) {
                    std::vector<T> chunk_to_process;
                    
                    // Get a chunk from the read queue
                    {
                        std::unique_lock<std::mutex> lock(read_mutex);
                        read_cv.wait(lock, [&]() {
                            return !read_queue.empty() || (done_reading && read_queue.empty());
                        });
                        
                        if (read_queue.empty() && done_reading) {
                            // No more work
                            if (active_processing_threads == 0) {
                                done_processing = true;
                                result_cv.notify_all();
                            }
                            break;
                        }
                        
                        if (!read_queue.empty()) {
                            chunk_to_process = std::move(read_queue.front());
                            read_queue.pop();
                            active_processing_threads++;
                            
                            lock.unlock();
                            read_cv.notify_one();
                        }
                    }
                    
                    if (!chunk_to_process.empty()) {
                        // Process the chunk using the provided function
                        std::string file_path = input_path;
                        chunk_func(chunk_to_process, file_path);
                        
                        active_processing_threads--;
                    }
                }
            });
        }
        
        // Wait for threads to finish
        reader_thread.join();
        for (auto& thread : worker_threads) {
            thread.join();
        }
    }

    int64_t total_sequences() const { return seqs_processed; }

private:
    std::atomic<size_t> chunk_size;
    std::atomic<int64_t> seqs_processed{0};
};

//working on optimizing the chunk streaming to be better
/*
template<typename T, typename function_to_run>
class chunk_streaming {
public:
    chunk_streaming(size_t chunk_size = 50000) 
        : chunk_size(chunk_size), seqs_processed(0) {}
    
    void modify_chunk_size(size_t new_size) {
        chunk_size = new_size;
    }
    
    size_t get_chunk_size() const {
        return chunk_size;
    }

    void process_chunks(const std::string& input_path, function_to_run chunk_func, int num_threads, int64_t max_seqs = -1) {
        // Create a file reader
        file_streaming files(input_path);
        read_streaming reader(files);
        
        // Create a pool of worker threads
        std::vector<std::thread> workers;
        std::atomic<bool> processing_complete{false};
        
        // Mutex and condition variable for synchronization
        std::mutex chunk_mutex;
        std::condition_variable chunk_cv;
        
        // Current chunk being processed
        std::vector<T> current_chunk;
        current_chunk.reserve(chunk_size);
        
        // Flag to indicate if chunk is ready
        bool chunk_ready = false;
        bool reading_complete = false;
        
        // Thread to read chunks
        std::thread reader_thread([&]() {
            while (true) {
                // Check if we've reached the maximum number of sequences
                if (max_seqs > 0 && seqs_processed >= max_seqs) {
                    break;
                }
                
                // Read next sequence
                auto item = reader.next_sequence();
                if (!item) {
                    // No more sequences
                    break;
                }
                
                // Add to current chunk
                {
                    std::lock_guard<std::mutex> lock(chunk_mutex);
                    current_chunk.push_back(*item);
                    seqs_processed++;
                }
                
                // If chunk is full, signal processing
                if (current_chunk.size() >= chunk_size) {
                    {
                        std::unique_lock<std::mutex> lock(chunk_mutex);
                        chunk_ready = true;
                        lock.unlock();
                        chunk_cv.notify_one();
                    }
                    
                    // Wait until chunk is processed
                    {
                        std::unique_lock<std::mutex> lock(chunk_mutex);
                        chunk_cv.wait(lock, [&]() { return !chunk_ready; });
                    }
                }
            }
            
            // Signal that the last chunk is ready (if not empty)
            {
                std::lock_guard<std::mutex> lock(chunk_mutex);
                reading_complete = true;
                if (!current_chunk.empty()) {
                    chunk_ready = true;
                }
            }
            chunk_cv.notify_one();
        });
        
        // Process chunks as they become available
        {
            std::unique_lock<std::mutex> lock(chunk_mutex);
            
            while (true) {
                // Wait for a chunk to be ready
                chunk_cv.wait(lock, [&]() { return chunk_ready || (reading_complete && !chunk_ready); });
                
                if (!chunk_ready && reading_complete) {
                    // No more chunks to process
                    break;
                }
                
                // Process the current chunk
                std::vector<T> chunk_to_process = std::move(current_chunk);
                current_chunk = std::vector<T>();
                current_chunk.reserve(chunk_size);
                chunk_ready = false;
                
                // Signal that we've taken the chunk
                lock.unlock();
                chunk_cv.notify_one();
                
                // Process the chunk in parallel
                #pragma omp parallel num_threads(num_threads)
                {
                    #pragma omp single
                    {
                        std::string current_file_path = files.current_file_path();
                        chunk_func(chunk_to_process, current_file_path);
                    }
                }
                
                // Re-acquire the lock
                lock.lock();
            }
        }
        
        // Wait for reader thread to finish
        if (reader_thread.joinable()) {
            reader_thread.join();
        }
    }

    int64_t total_sequences() const { return seqs_processed; }

private:
    size_t chunk_size;
    std::atomic<int64_t> seqs_processed{0};
};*/

class sigstring_writing {
    public:
        enum class format {
            FASTQA,
            SIGSTRING,
            CSV
        };

        std::unique_ptr<std::ostream> out;
        format fmt;
        
        // @param output_path   base filename ("myout"), or with .gz
        // @param output_format which of the three formats to emit
        // @param compress      if true, writes to `output_path + ".gz"` via ogzstream
        // @param append        if true, opens in append mode; otherwise truncates
        sigstring_writing(std::string output_path, format output_format, bool compress, bool append) : fmt(output_format) {
            std::ios::openmode mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
            if (compress) {
                // ensure .gz suffix
                if (output_path.size()<3 || output_path.substr(output_path.size()-3) != ".gz")
                    output_path += ".gz";
                auto gz = std::make_unique<ogzstream>(output_path.c_str(), mode);
                if (!gz->good())
                    throw std::runtime_error("Failed to open gz output: " + output_path);
                out = std::move(gz);
            } else {
                auto f = std::make_unique<std::ofstream>(output_path, mode);
                if (!f->is_open())
                    throw std::runtime_error("Failed to open file: " + output_path);
                out = std::move(f);
            }
        }
        // Write all items in the chunk
        template<typename T>
        void operator()(const std::vector<T>& chunk) {
            for (const auto& item : chunk) {
                switch (fmt) {
                    case format::FASTQA:
                        *out << item.to_fastqa() << "\n";
                        break;
                    case format::SIGSTRING:
                        *out << item.to_sigstring() << "\n";
                        break;
                    case format::CSV:
                        *out << item.to_csv();
                        break;
                }
            }
        }

        ~sigstring_writing() {
            if (out) {
                out->flush();
            }
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
        void write(sigstring_writing& sig_writer, 
                   sigstring_writing& csv_writer, 
                   sigstring_writing& fastqa_writer,
                   std::vector<T>& data) {
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
            job->size = data.size(); // Record size before moving
            
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
                        std::cout << "Writer: chunk " << job->chunk_id 
                                  << " (" << job->size << " entries) written in " 
                                  << write_time_ms << " ms (queued for " 
                                  << queue_time_ms << " ms)" << std::endl;
                    }
                    
                    // Release the job (shared_ptr will clean up the memory)
                    job.reset();
                }
            }
        }
    };