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

template<typename T, typename function_to_run>
class chunk_streaming {
public:
    chunk_streaming(size_t chunk_size = 50000) 
        : chunk_size(chunk_size), seqs_processed(0) {}

    void process_chunks(const std::string& input_path, function_to_run chunk_func, int num_threads = 1, int64_t max_seqs = -1) {
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
};

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