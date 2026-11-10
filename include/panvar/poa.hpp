#pragma once

#include <string>
#include <vector>

namespace panvar {

// Multiple-sequence alignment via abPOA (partial-order alignment, affine gaps). Returns gapped rows
// (equal length, '-' for gaps) in the same order as the input sequences. Each row holds its sequence's own
// bytes (letter case and non-ACGT bases kept); abPOA aligns them as upper case with non-ACGT read as N.
// n==0 -> {}; n==1 -> the single sequence unchanged; empty input sequences become all-gap rows.
std::vector<std::string> poa_msa(const std::vector<std::string>& seqs);

} // namespace panvar
