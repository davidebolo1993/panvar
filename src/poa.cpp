#include "panvar/poa.hpp"

#include <cstdint>
#include <cstdlib>

extern "C" {
#include "abpoa.h"
}

namespace panvar {
namespace {

inline uint8_t nt_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // N / anything else
    }
}

// abPOA RC-MSA base codes: 0-3 = ACGT, 4 = N, 5 = gap (== abpt->m for the nt alphabet).
const char kDecode[] = "ACGTN-";

} // namespace

std::vector<std::string> poa_msa(const std::vector<std::string>& seqs) {
    const std::size_t n = seqs.size();
    if (n == 0) return {};
    if (n == 1) return {seqs[0]};

    // abPOA cannot take zero-length reads: align the non-empty sequences and re-insert any empty
    // ones as all-gap rows of the resulting MSA length.
    std::vector<std::size_t> nonempty;
    for (std::size_t i = 0; i < n; ++i) {
        if (!seqs[i].empty()) nonempty.push_back(i);
    }
    std::vector<std::string> rows(n);
    if (nonempty.empty()) {
        return std::vector<std::string>(n, std::string());
    }
    if (nonempty.size() == 1) {
        const std::string& s = seqs[nonempty[0]];
        for (std::size_t i = 0; i < n; ++i) {
            rows[i] = (i == nonempty[0]) ? s : std::string(s.size(), '-');
        }
        return rows;
    }

    const int m = static_cast<int>(nonempty.size());
    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->out_msa = 1;
    abpt->out_cons = 0;
    abpt->out_gfa = 0;
    abpoa_post_set_para(abpt);

    int* seq_lens = static_cast<int*>(std::malloc(sizeof(int) * m));
    uint8_t** bseqs = static_cast<uint8_t**>(std::malloc(sizeof(uint8_t*) * m));
    for (int i = 0; i < m; ++i) {
        const std::string& s = seqs[nonempty[i]];
        const int len = static_cast<int>(s.size());
        seq_lens[i] = len;
        bseqs[i] = static_cast<uint8_t*>(std::malloc(static_cast<std::size_t>(len)));
        for (int j = 0; j < len; ++j) bseqs[i][j] = nt_code(s[j]);
    }

    abpoa_msa(ab, abpt, m, nullptr, seq_lens, bseqs, nullptr, nullptr);

    const abpoa_cons_t* abc = ab->abc;
    const int msa_len = (abc != nullptr) ? abc->msa_len : 0;
    for (std::size_t i = 0; i < n; ++i) rows[i].assign(static_cast<std::size_t>(msa_len), '-');
    if (abc != nullptr && msa_len > 0 && abc->msa_base != nullptr) {
        for (int i = 0; i < m; ++i) {
            std::string& r = rows[nonempty[i]];
            for (int j = 0; j < msa_len; ++j) {
                const uint8_t b = abc->msa_base[i][j];
                r[static_cast<std::size_t>(j)] = (b <= 5) ? kDecode[b] : '-';
            }
        }
    }

    for (int i = 0; i < m; ++i) std::free(bseqs[i]);
    std::free(bseqs);
    std::free(seq_lens);
    abpoa_free(ab);
    abpoa_free_para(abpt);
    return rows;
}

} // namespace panvar
