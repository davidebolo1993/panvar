#pragma once
// M15-L: the genotype command's output layout. A command run with --out P keeps only its primary reports beside P
// (.sequence_locus.tsv, .sequence_blocks.tsv, .sequence_haplotypes.tsv and the conditional .sequence_haplotypes.fa); every
// other report goes to P.diagnostics/<basename P>.<suffix>. The run works under that diagnostics prefix D, registered here
// against P, and every writer, reader or remover of a primary file asks genotype_output_path for its path. An unregistered
// prefix (the internal API without diagnostics_directory) keeps the flat layout: every path is prefix + suffix.
#include <string>

namespace panvar {

// The four primary suffixes, and their atomic-write temporaries (<suffix>.writing).
bool genotype_primary_suffix(const std::string& suffix);
// P -> P.diagnostics/<basename P>
std::string genotype_diagnostics_prefix(const std::string& primary_prefix);
// Register / forget the run's layout (one run at a time per process).
void set_genotype_output_layout(const std::string& diagnostics_prefix, const std::string& primary_prefix);
void clear_genotype_output_layout();
// True when `prefix` is the registered diagnostics prefix.
bool genotype_layout_active(const std::string& prefix);
// The path of `suffix` for a run working under `prefix`: the primary prefix for a primary suffix of a registered run,
// otherwise prefix + suffix.
std::string genotype_output_path(const std::string& prefix, const std::string& suffix);
// The user's prefix P for a registered diagnostics prefix; `prefix` itself otherwise.
std::string genotype_primary_prefix(const std::string& prefix);

// Registers the layout for one command run and forgets it when the run ends, whatever the exit path.
class GenotypeOutputLayoutScope {
public:
    GenotypeOutputLayoutScope(const std::string& diagnostics_prefix, const std::string& primary_prefix);
    ~GenotypeOutputLayoutScope();
    GenotypeOutputLayoutScope(const GenotypeOutputLayoutScope&) = delete;
    GenotypeOutputLayoutScope& operator=(const GenotypeOutputLayoutScope&) = delete;
};

}  // namespace panvar
