#include "panvar/genotype_output_layout.hpp"

#include <mutex>

namespace panvar {
namespace {
std::mutex g_layout_mutex;
std::string g_diagnostics_prefix, g_primary_prefix;   // empty: no layout registered (flat)
}  // namespace

bool genotype_primary_suffix(const std::string& suffix) {
    static const char* const primary[] = {".sequence_locus.tsv", ".sequence_blocks.tsv", ".sequence_haplotypes.tsv",
                                          ".sequence_haplotypes.fa"};
    for (const char* p : primary) {
        if (suffix == p || suffix == std::string(p) + ".writing") return true;
    }
    return false;
}

std::string genotype_diagnostics_prefix(const std::string& primary_prefix) {
    const std::size_t slash = primary_prefix.find_last_of('/');
    const std::string base = slash == std::string::npos ? primary_prefix : primary_prefix.substr(slash + 1);
    return primary_prefix + ".diagnostics/" + base;
}

void set_genotype_output_layout(const std::string& diagnostics_prefix, const std::string& primary_prefix) {
    std::lock_guard<std::mutex> lock(g_layout_mutex);
    g_diagnostics_prefix = diagnostics_prefix;
    g_primary_prefix = primary_prefix;
}

void clear_genotype_output_layout() {
    std::lock_guard<std::mutex> lock(g_layout_mutex);
    g_diagnostics_prefix.clear();
    g_primary_prefix.clear();
}

bool genotype_layout_active(const std::string& prefix) {
    std::lock_guard<std::mutex> lock(g_layout_mutex);
    return !g_diagnostics_prefix.empty() && prefix == g_diagnostics_prefix;
}

std::string genotype_output_path(const std::string& prefix, const std::string& suffix) {
    {
        std::lock_guard<std::mutex> lock(g_layout_mutex);
        if (!g_diagnostics_prefix.empty() && prefix == g_diagnostics_prefix && genotype_primary_suffix(suffix)) {
            return g_primary_prefix + suffix;
        }
    }
    return prefix + suffix;
}

std::string genotype_primary_prefix(const std::string& prefix) {
    std::lock_guard<std::mutex> lock(g_layout_mutex);
    return !g_diagnostics_prefix.empty() && prefix == g_diagnostics_prefix ? g_primary_prefix : prefix;
}

GenotypeOutputLayoutScope::GenotypeOutputLayoutScope(const std::string& diagnostics_prefix, const std::string& primary_prefix) {
    set_genotype_output_layout(diagnostics_prefix, primary_prefix);
}
GenotypeOutputLayoutScope::~GenotypeOutputLayoutScope() { clear_genotype_output_layout(); }

}  // namespace panvar
