#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "panvar/associate_command.hpp"
#include "panvar/benchmark_command.hpp"
#include "panvar/bubble_command.hpp"
#include "panvar/call_command.hpp"
#include "panvar/cli_utils.hpp"
#ifdef PANVAR_ENABLE_EXPERIMENTAL_GENOTYPE
#include "panvar/genotype_command.hpp"
#include "panvar/genotype_frag_command.hpp"
#endif
#include "panvar/describe_command.hpp"
#include "panvar/inspect.hpp"
#include "panvar/panphorte_command.hpp"
#include "panvar/rebuild_command.hpp"
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
        if (subcommand == "genotype") {
#ifdef PANVAR_ENABLE_EXPERIMENTAL_GENOTYPE
            return panvar::run_genotype_command(args);
#else
            // Named explicitly rather than falling through to "Unknown subcommand": the module exists
            // in the tree and in the history, so a user who has heard of it deserves to be told it was
            // withheld and why, not that it never existed.
            throw std::runtime_error(
                "genotype is not built in this release. It is the one module that has not completed "
                "its review pass, so it is excluded rather than shipped unreviewed. To build it for "
                "development: cmake -DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON");
#endif
        }
        if (subcommand == "genotype-frag") {
#ifdef PANVAR_ENABLE_EXPERIMENTAL_GENOTYPE
            return panvar::run_genotype_frag_command(args);
#else
            throw std::runtime_error(
                "genotype-frag is not built in this release. It is a prototype of the genotype "
                "module's fragment-level evidence model, gated with it. To build it for "
                "development: cmake -DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON");
#endif
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
        if (subcommand == "rebuild") {
            return panvar::run_rebuild_command(args);
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
