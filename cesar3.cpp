#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
#include <filesystem>
#include <system_error>
#include <cstring>
#include <csignal>
#include <algorithm>
#include <fcntl.h>
namespace fs = std::filesystem;

extern char **environ;
std::string shm_dir;

void signal_handler(int signum) {
    if (!shm_dir.empty()) {
        std::error_code ec;
        fs::remove_all(shm_dir, ec);
    }
    std::exit(signum);
}

void print_usage() {
    std::cerr << "cesar3 — CESAR 3.0 Multi-Species Parallel Runner\n"
              << "Copyright 2017 MPI-CBG/MPI-PKS Peter Schwede\n\n"
              << "Align a pair of sequences using HMM.\n"
              << "Drop-in upgrade for cesar with automatic parallelism.\n\n"
              << "Usage: cesar3 <fasta file>\n"
              << "               -c/--clade <human|mouse>\n"
              << "               -m/--matrix <matrix file>\n"
              << "               -p/--profiles <acc profile> <do profile>\n"
              << "               -i/--split_codon_emissions <acc split codon emissions> <do split codon emissions>\n"
              << "               -f/--firstexon\n"
              << "               -l/--lastexon\n"
              << "               -x/--max-memory\n"
              << "               -v/--verbosity <0 .. 2>\n"
              << "               -V/--version\n"
              << "               -s/--set name1=value1 .. nameN=valueN\n"
              << "               -S/--sanityChecks  (error-exit if sanity checks fail)\n"
              << "               -h/--help\n"
              << "               -t/--threads <N>     (parallel workers, default: 80% of cores)\n"
              << "               -o/--output <dir>    (output dir, required for >1 species)\n"
              << "               --no-shm             (use /tmp instead of /dev/shm)\n\n"
              << "Single species → stdout (same as cesar 2.0):\n"
              << "  ./cesar3 gene.fa -c human\n\n"
              << "Multi species → parallel, RAM-backed:\n"
              << "  ./cesar3 ENST00000262209.fa -o results/ -t 20\n\n"
              << "cesar binary auto-detected from same directory as cesar3.\n";
}

void get_spawn_output(const std::string& bin, const std::string& input_file, const std::vector<std::string>& pt_args, const std::string& err_file, std::string& out_str, std::string& err_str) {
    int pipefd[2];
    if (pipe(pipefd) != 0) {
        return;
    }

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_addclose(&actions, pipefd[0]);
    posix_spawn_file_actions_adddup2(&actions, pipefd[1], STDOUT_FILENO);
    posix_spawn_file_actions_addclose(&actions, pipefd[1]);
    posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, err_file.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);

    pid_t pid;
    std::vector<char*> argv_vec;
    argv_vec.push_back(strdup(bin.c_str()));
    argv_vec.push_back(strdup(input_file.c_str()));
    for (const auto& a : pt_args) {
        argv_vec.push_back(strdup(a.c_str()));
    }
    argv_vec.push_back(nullptr);
    
    if (posix_spawn(&pid, bin.c_str(), &actions, nullptr, argv_vec.data(), environ) != 0) {
        posix_spawn_file_actions_destroy(&actions);
        close(pipefd[0]);
        close(pipefd[1]);
        for (char* p : argv_vec) free(p);
        return;
    }

    posix_spawn_file_actions_destroy(&actions);
    close(pipefd[1]);
    for (char* p : argv_vec) free(p);

    char buffer[4096];
    ssize_t count;
    while ((count = read(pipefd[0], buffer, sizeof(buffer))) > 0) {
        out_str.append(buffer, count);
    }
    close(pipefd[0]);

    int status;
    waitpid(pid, &status, 0);

    std::ifstream err_in(err_file);
    if (err_in.is_open()) {
        std::string err_line;
        while (std::getline(err_in, err_line)) {
            err_str += err_line + "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    std::string arg1 = argv[1];
    if (arg1 == "-h" || arg1 == "--help") {
        print_usage();
        unsigned int hw_threads = std::thread::hardware_concurrency();
        int def = hw_threads > 0 ? std::max(1, (int)(0.8 * hw_threads)) : 4;
        std::cerr << "\nDetected default threads: " << def << "\n";
        return 0;
    }

    std::string input_fa = argv[1];
    
    fs::path exe_path = argv[0];
    fs::path bin_dir = exe_path.parent_path();
    if (bin_dir.empty()) bin_dir = ".";
    std::string cesar_bin = (bin_dir / "cesar").string();

    std::string out_dir = "";
    
    unsigned int hw_threads = std::thread::hardware_concurrency();
    int num_threads = hw_threads > 0 ? std::max(1, (int)(0.8 * hw_threads)) : 4;
    bool use_shm = true;
    int verbosity = 1;
    std::vector<std::string> passthrough_args;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                out_dir = argv[++i];
            }
        } else if (arg == "--no-shm") {
            use_shm = false;
        } else if (arg == "-v" || arg == "--verbosity") {
            if (i + 1 < argc) {
                verbosity = std::stoi(argv[i+1]);
                passthrough_args.push_back(arg);
                passthrough_args.push_back(argv[++i]);
            }
        } else if (arg == "-h" || arg == "--help") {
            print_usage();
            unsigned int hw_threads = std::thread::hardware_concurrency();
            int def = hw_threads > 0 ? std::max(1, (int)(0.8 * hw_threads)) : 4;
            std::cerr << "\nDetected default threads: " << def << "\n";
            return 0;
        } else {
            passthrough_args.push_back(arg);
        }
    }

    std::ifstream in(input_fa);
    if (!in.is_open()) {
        std::cerr << "ERROR: Cannot open input file " << input_fa << "\n";
        return 1;
    }

    std::string ref_block;
    std::string line;
    bool in_species = false;
    std::vector<std::string> species_blocks;
    std::vector<std::string> species_names;
    std::string current_species;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (!in_species) {
            ref_block += line + "\n";
            if (line == "####") {
                in_species = true;
            }
        } else {
            if (line.rfind(">", 0) == 0) {
                if (!current_species.empty()) {
                    species_blocks.push_back(current_species);
                    size_t end_pos = current_species.find_first_of(" \n\r\t");
                    species_names.push_back(current_species.substr(1, end_pos - 1));
                }
                current_species = line + "\n";
            } else {
                current_species += line + "\n";
            }
        }
    }
    if (!current_species.empty()) {
        species_blocks.push_back(current_species);
        size_t end_pos = current_species.find_first_of(" \n\r\t");
        species_names.push_back(current_species.substr(1, end_pos - 1));
    }

    int species_count = species_blocks.size();

    if (species_count == 0) {
        std::cerr << "ERROR: No query species found after #### separator.\n";
        return 1;
    }

    if (species_count == 1) {
        std::vector<char*> exec_argv;
        exec_argv.push_back(strdup(cesar_bin.c_str()));
        exec_argv.push_back(strdup(input_fa.c_str()));
        for (const auto& a : passthrough_args) {
            exec_argv.push_back(strdup(a.c_str()));
        }
        exec_argv.push_back(nullptr);
        
        execv(cesar_bin.c_str(), exec_argv.data());
        std::cerr << "ERROR: execv failed for single species mode.\n";
        for (char* p : exec_argv) free(p);
        return 1;
    }

    if (out_dir.empty()) {
        std::cerr << "ERROR: -o/--output DIR is required when fasta contains >1 species.\n";
        return 1;
    }

    std::cerr << "[INFO] " << fs::path(input_fa).filename().string() << ": " << species_count << " species, " << num_threads << " threads\n";

    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);

    std::string pid_str = std::to_string(getpid());
    
    if (use_shm && ::access("/dev/shm", W_OK) == 0) {
        shm_dir = "/dev/shm/cesar3_" + pid_str;
    } else {
        shm_dir = "/tmp/cesar3_" + pid_str;
    }

    std::error_code ec;
    fs::create_directories(shm_dir, ec);

    fs::copy_file(input_fa, shm_dir + "/input.fa", fs::copy_options::overwrite_existing, ec);

    std::mutex q_mtx;
    std::queue<int> job_queue;
    for (int i = 0; i < species_count; ++i) {
        job_queue.push(i);
    }

    std::vector<std::string> results(species_count);
    std::vector<std::string> warned_species;
    std::mutex cerr_mtx;
    int completed = 0;

    auto worker = [&]() {
        while (true) {
            int job_id;
            {
                std::lock_guard<std::mutex> lock(q_mtx);
                if (job_queue.empty()) break;
                job_id = job_queue.front();
                job_queue.pop();
            }

            std::string sp_file = shm_dir + "/sp_" + std::to_string(job_id) + ".fa";
            std::string err_file = shm_dir + "/err_" + std::to_string(job_id) + ".err";
            {
                std::ofstream out(sp_file);
                out << ref_block << species_blocks[job_id];
            }

            std::string out_str, err_str;
            get_spawn_output(cesar_bin, sp_file, passthrough_args, err_file, out_str, err_str);
            results[job_id] = std::move(out_str);

            {
                std::lock_guard<std::mutex> lock(cerr_mtx);
                completed++;
                if (verbosity > 0) {
                    std::cerr << "[INFO] species " << completed << "/" << species_count << " done\n";
                }
                if (!err_str.empty()) {
                    if (verbosity == 1) {
                        warned_species.push_back(species_names[job_id]);
                    } else if (verbosity >= 2) {
                        size_t start = 0;
                        size_t end;
                        while ((end = err_str.find('\n', start)) != std::string::npos) {
                            std::cerr << "[" << species_names[job_id] << "] " << err_str.substr(start, end - start + 1);
                            start = end + 1;
                        }
                        if (start < err_str.length()) {
                            std::cerr << "[" << species_names[job_id] << "] " << err_str.substr(start) << "\n";
                        }
                    }
                }
            }
        }
    };

    std::vector<std::thread> workers;
    for (int w = 0; w < num_threads; ++w) {
         workers.emplace_back(worker);
    }
    for (auto& w : workers) {
         w.join();
    }

    if (verbosity == 1 && !warned_species.empty()) {
        std::cerr << "[WARN] N literals found in " << warned_species.size() << " species: ";
        for (size_t i = 0; i < warned_species.size(); ++i) {
            std::cerr << warned_species[i] << (i + 1 == warned_species.size() ? "" : ", ");
        }
        std::cerr << "\n";
    }

    std::string gene_name = fs::path(input_fa).stem().string();
    std::string out_path = out_dir + "/" + gene_name + ".out";
    
    fs::create_directories(out_dir, ec);
    std::ofstream final_out(out_path);
    for (const auto& r : results) {
        final_out << r;
    }

    std::cerr << "[DONE] " << out_path << "\n";
    fs::remove_all(shm_dir, ec);

    return 0;
}
