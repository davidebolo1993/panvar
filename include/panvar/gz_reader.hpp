#pragma once

#include <cstddef>
#include <string>

#include <zlib.h>

namespace panvar {

// Transparent line reader for plain or gzipped text. zlib's gzread auto-detects
// the gzip magic, so the same reader handles `foo.gfa` and `foo.gfa.gz`. Lines
// are returned without the trailing '\n' (and a trailing '\r' is stripped), and
// arbitrarily long lines (e.g. GFA P-lines) are supported.
class GzLineReader {
public:
    explicit GzLineReader(const std::string& path) {
        file_ = gzopen(path.c_str(), "rb");
    }
    ~GzLineReader() {
        if (file_ != nullptr) {
            gzclose(file_);
        }
    }
    GzLineReader(const GzLineReader&) = delete;
    GzLineReader& operator=(const GzLineReader&) = delete;

    bool ok() const { return file_ != nullptr; }

    bool getline(std::string& out) {
        std::size_t from = pos_;
        while (true) {
            const std::size_t nl = buf_.find('\n', from);
            if (nl != std::string::npos) {
                out.assign(buf_, pos_, nl - pos_);
                if (!out.empty() && out.back() == '\r') {
                    out.pop_back();
                }
                pos_ = nl + 1;
                if (pos_ > (std::size_t{1} << 20)) {
                    buf_.erase(0, pos_);
                    pos_ = 0;
                }
                return true;
            }
            from = buf_.size();
            char chunk[1 << 16];
            const int n = gzread(file_, chunk, static_cast<unsigned int>(sizeof(chunk)));
            if (n > 0) {
                buf_.append(chunk, static_cast<std::size_t>(n));
                continue;
            }
            // EOF (or error): emit any trailing unterminated line, then stop.
            if (pos_ < buf_.size()) {
                out.assign(buf_, pos_, buf_.size() - pos_);
                if (!out.empty() && out.back() == '\r') {
                    out.pop_back();
                }
                pos_ = buf_.size();
                return true;
            }
            return false;
        }
    }

private:
    gzFile file_ = nullptr;
    std::string buf_;
    std::size_t pos_ = 0;
};

} // namespace panvar
