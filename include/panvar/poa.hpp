#pragma once

#include <string>
#include <vector>

namespace panvar {

// Multiple-sequence alignment via abPOA (partial-order alignment, affine gaps). Returns gapped rows
// (equal length, '-' for gaps) in the same order as the input sequences. Non-ACGT bases map to 'N'.
// n==0 -> {}; n==1 -> the single sequence unchanged; empty input sequences become all-gap rows.
std::vector<std::string> poa_msa(const std::vector<std::string>& seqs);

} // namespace panvar
