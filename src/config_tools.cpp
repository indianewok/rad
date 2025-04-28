#include "include/rad/rad_headers.h"

void usage(const char *prog) {
    std::cerr
      << "Usage:\n"
      << "  " << prog << " layout list\n"
      << "      list all registered layout types and their .csv paths (defaults should all be in resources/read_layout) \n\n"
      << "  " << prog << " layout get [type]\n"
      << "      print the .csv path for layout type [type] \n\n"
      << "  " << prog << " layout set [type] [path]\n"
      << "      register or override the path for a layout type [type]\n\n"
      << "  " << prog << " layout rm [type]\n"
      << "      remove the layout type [type]\n\n"
      << "  " << prog << " whitelist list\n"
      << "      list all known whitelist file paths (defaults should all be in resources/wl) \n\n"
      << "  " << prog << " whitelist get [kit]\n"
      << "      print the whitelist path for [kit] \n\n"
      << "  " << prog << " whitelist set [kit] [path]\n"
      << "      set or override the global whitelist file at [path] for [kit]\n\n"
      << "  " << prog << " whitelist rm [kit]\n"
      << "      remove whitelist so it's no longer affiliated with [kit] \n\n"
      << "Flags:\n"
      << "  -h, --help    Show this message and exit\n";
  }

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    std::string domain = argv[1];
    std::string cmd = argv[2];

    try {
        if (domain == "layout") {
            if (cmd == "list") {
                for (auto &kv : config_utils::layout_files)
                    std::cout << kv.first << " => " << kv.second << "\n";
            }
            else if (cmd == "get" && argc == 4) {
                std::cout << config_utils::get_read_layout(argv[3]) << "\n";
            }
            else if (cmd == "set" && argc == 5) {
                config_utils::save_read_layout(argv[3], argv[4]);
                std::cout << "Set layout[" << argv[3] << "] = " << argv[4] << "\n";
            }
            else if (cmd == "rm" && argc == 4) {
                bool ok = config_utils::remove_read_layout(argv[3]);
                std::cout << (ok ? "Removed\n" : "Not found\n");
            }
            else { usage(argv[0]); return 1; }
        }
        else if (domain == "whitelist") {
            if (cmd == "list") {
                for (auto &kv : whitelist_utils::kit_wl_paths)
                    std::cout << kv.first << " => " << kv.second << "\n";
            }
            else if (cmd == "get" && argc == 4) {
                std::cout << whitelist_utils::get_whitelist_path(argv[3]) << "\n";
            }
            else if (cmd == "set" && argc == 5) {
                whitelist_utils::set_whitelist_path(argv[3], argv[4]);
                std::cout << "Set whitelist[" << argv[3] << "] = " << argv[4] << "\n";
            }
            else if (cmd == "rm" && argc == 4) {
                bool ok = whitelist_utils::remove_whitelist_path(argv[3]);
                std::cout << (ok ? "Removed\n" : "Not found\n");
            }
            else { usage(argv[0]); return 1; }
        }
        else {
            usage(argv[0]);
            return 1;
        }
    }
    catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}