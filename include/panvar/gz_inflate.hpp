#pragma once

#include <cstddef>

namespace panvar {

// One gzip stream's decompressor, for GzLineReader: zlib-ng's own API when the build has it (PANVAR_ZLIB_NG; inflate
// about twice as fast on x86), else the system zlib. Its header names no zlib type, because zlib-ng's header cannot
// share a translation unit with zlib's, which other code includes.
class Inflater {
public:
    enum Status { kOk, kMemberEnd, kNeedInput, kError };
    Inflater();   // throws std::runtime_error when it cannot start
    ~Inflater();
    Inflater(const Inflater&) = delete;
    Inflater& operator=(const Inflater&) = delete;
    // Decompresses from `in` into `out`; `used` and `produced` are the bytes taken from each.
    Status run(const unsigned char* in, std::size_t in_len, char* out, std::size_t cap, std::size_t& used, std::size_t& produced);
    void next_member();          // ready for the gzip member that follows
    const char* message() const;  // the library's message after kError, if any
    static const char* library(); // "zlib-ng <version>" or "zlib <version>"

private:
    void* stream_;
};

} // namespace panvar
