#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "panvar/associate_command.hpp"
#include "panvar/benchmark_command.hpp"
#include "panvar/bubble_command.hpp"
#include "panvar/call_command.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/describe_command.hpp"
#include "panvar/inspect.hpp"
#include "panvar/panphorte_command.hpp"
#include "panvar/refine_command.hpp"

int main(int argc, char** argv) {
    try {
        if (argc <= 1) {
            panvar::cli::print_general_help();
            return 1;
        }

        const std::string subcommand = argv[1];

        std::vector<std::string> args;
        args.reserve(static_cast<std::size_t>(argc > 2 ? argc - 2 : 0));
        for (int i = 2; i < argc; ++i) {
            args.emplace_back(argv[i]);
        }

        if (subcommand == "-h" || subcommand == "--help") {
            panvar::cli::print_general_help();
            return 0;
        }

        if (subcommand == "bubble") {
            return panvar::run_bubble_command(args);
        }
        if (subcommand == "inspect") {
            return panvar::run_inspect_command(args);
        }
        if (subcommand == "call") {
            return panvar::run_call_command(args);
        }
        if (subcommand == "benchmark") {
            return panvar::run_benchmark_command(args);
        }
        if (subcommand == "describe") {
            return panvar::run_describe_command(args);
        }
        if (subcommand == "panphorte") {
            return panvar::run_panphorte_command(args);
        }
        if (subcommand == "refine") {
            return panvar::run_refine_command(args);
        }
        if (subcommand == "associate") {
            return panvar::run_associate_command(args);
        }

        throw std::runtime_error("Unknown subcommand: " + subcommand);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }
}
