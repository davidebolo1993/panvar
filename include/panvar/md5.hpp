#pragma once

#include <cstddef>
#include <string>

namespace panvar {

// Lowercase hex MD5 of a byte range. Present so a dumped sequence can be checked against the
// external assembly it is supposed to reproduce with plain `md5sum`, without a second tool in the
// loop -- the audit that found the LPA representation drift had to hash in Python, which made the
// caller's own numbers unverifiable from the caller's own output.
std::string md5_hex(const void* data, std::size_t len);

inline std::string md5_hex(const std::string& s) { return md5_hex(s.data(), s.size()); }

} // namespace panvar
