#pragma once

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>

#include <zlib.h>
#include <sys/stat.h>
#include <unistd.h>

namespace panvar {

// Transparent line reader for plain or gzipped text: the gzip magic is detected from the first bytes,
// so the same reader handles `foo.gfa` and `foo.gfa.gz`, from a file or a pipe. Lines are returned
// without the trailing '\n' (and a trailing '\r' is stripped), and arbitrarily long lines (e.g. GFA
// P-lines) are supported.
//
// A damaged compressed file throws rather than reading as a shorter valid one. GFA writers put the
// path lines last, so a truncated graph otherwise loses haplotypes and parses cleanly. Decompression
// drives zlib's inflate directly, member by member, rather than gzread: gzread takes anything after a
// complete member that is not a gzip header as trailing garbage and ends the file there, so one damaged
// block header in a BGZF file, or in any multi-member .gz, would silently drop everything after it.
class GzLineReader {
public:
    explicit GzLineReader(const std::string& path) : path_(path) {
        std::error_code ec;
        if (std::filesystem::is_directory(path, ec)) {
            throw std::runtime_error(path + ": is a directory, not a file");
        }
        fp_ = std::fopen(path.c_str(), "rb");
        if (fp_ == nullptr) return;
        refill();
        gz_ = in_len_ - in_pos_ >= 2 && in_[in_pos_] == 0x1f && in_[in_pos_ + 1] == 0x8b;
        if (gz_) {
            std::memset(&zs_, 0, sizeof(zs_));
            if (inflateInit2(&zs_, 15 + 16) != Z_OK) {
                close();
                throw std::runtime_error(path_ + ": cannot initialise decompression");
            }
            zs_open_ = true;
            if (is_bgzf_missing_eof_marker(fp_)) {
                close();
                throw std::runtime_error(
                    path + ": BGZF file has no end-of-file marker, so it was cut short (an interrupted "
                           "write ends on a block boundary, which decompresses without error). If the "
                           "file is known to be complete, decompress it and re-compress with bgzip");
            }
        }
    }
    ~GzLineReader() { close(); }
    GzLineReader(const GzLineReader&) = delete;
    GzLineReader& operator=(const GzLineReader&) = delete;

    bool ok() const { return fp_ != nullptr; }

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
            const std::size_t n = read_more(chunk, sizeof(chunk));
            if (n > 0) {
                buf_.append(chunk, n);
                continue;
            }
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
    void close() {
        if (zs_open_) {
            inflateEnd(&zs_);
            zs_open_ = false;
        }
        if (fp_ != nullptr) {
            std::fclose(fp_);
            fp_ = nullptr;
        }
    }

    // Move the unread input to the front and top the buffer up from the file.
    void refill() {
        if (in_pos_ > 0) {
            std::memmove(in_, in_ + in_pos_, in_len_ - in_pos_);
            in_len_ -= in_pos_;
            in_pos_ = 0;
        }
        while (!in_eof_ && in_len_ < sizeof(in_)) {
            const std::size_t got = std::fread(in_ + in_len_, 1, sizeof(in_) - in_len_, fp_);
            in_len_ += got;
            if (got == 0) {
                if (std::ferror(fp_)) {
                    throw std::runtime_error(path_ + ": read error (" + std::strerror(errno) + ")");
                }
                in_eof_ = true;
            }
            if (got > 0) break;
        }
    }

    // Up to `cap` decompressed (or plain) bytes; 0 only at a clean end of input.
    std::size_t read_more(char* out, std::size_t cap) {
        if (!gz_) {
            if (in_pos_ == in_len_) refill();
            const std::size_t n = std::min(cap, in_len_ - in_pos_);
            std::memcpy(out, in_ + in_pos_, n);
            in_pos_ += n;
            return n;
        }
        if (done_) return 0;
        while (true) {
            if (in_pos_ == in_len_) refill();
            zs_.next_in = in_ + in_pos_;
            zs_.avail_in = static_cast<uInt>(in_len_ - in_pos_);
            zs_.next_out = reinterpret_cast<Bytef*>(out);
            zs_.avail_out = static_cast<uInt>(cap);
            const int ret = inflate(&zs_, Z_NO_FLUSH);
            in_pos_ = in_len_ - zs_.avail_in;
            const std::size_t produced = cap - zs_.avail_out;
            if (ret == Z_STREAM_END) {
                next_member();
                if (produced > 0 || done_) return produced;
                continue;
            }
            if (ret == Z_OK || ret == Z_BUF_ERROR) {
                if (produced > 0) return produced;
                if (in_pos_ == in_len_ && in_eof_) {
                    throw std::runtime_error(path_ + ": compressed stream ends early, so the file is "
                                                     "truncated");
                }
                continue;
            }
            const std::string why = zs_.msg != nullptr ? zs_.msg : "data error";
            throw std::runtime_error(path_ + ": cannot decompress (zlib: " + why + ")");
        }
    }

    // After a member ends, what follows must be another gzip member (BGZF and concatenated .gz are
    // chains of them) or the end of the file. Zero bytes are tolerated as padding.
    void next_member() {
        while (true) {
            if (in_len_ - in_pos_ < 2) refill();
            while (in_pos_ < in_len_ && in_[in_pos_] == 0) ++in_pos_;
            if (in_pos_ < in_len_) break;
            if (in_eof_) {
                done_ = true;
                return;
            }
        }
        if (in_len_ - in_pos_ < 2) refill();
        if (in_len_ - in_pos_ >= 2 && in_[in_pos_] == 0x1f && in_[in_pos_ + 1] == 0x8b) {
            inflateReset(&zs_);
            return;
        }
        throw std::runtime_error(path_ + ": data after a complete gzip member is not another gzip "
                                         "member, so a block of the file is damaged; the rest of the "
                                         "file cannot be read");
    }

    // A BGZF file is a chain of independent gzip blocks ending in a fixed 28-byte empty block. A
    // write interrupted between blocks leaves a valid gzip stream, so only that marker shows it.
    // Checked on regular files only, through the descriptor already open, with pread, which moves no
    // file offset: reopening by name would move the offset of /dev/stdin redirected from a file (a dup
    // of fd 0 on macOS), and the reader would then decompress from the wrong place.
    static bool is_bgzf_missing_eof_marker(std::FILE* fp) {
        static const unsigned char kEof[28] = {
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
        const int fd = ::fileno(fp);
        struct stat st {};
        if (fd < 0 || ::fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size < 28) return false;
        unsigned char head[16] = {0};
        const bool bgzf = ::pread(fd, head, sizeof(head), 0) == static_cast<ssize_t>(sizeof(head)) &&
                          head[0] == 0x1f && head[1] == 0x8b && head[2] == 0x08 &&
                          (head[3] & 0x04) != 0 && head[12] == 'B' && head[13] == 'C';
        if (!bgzf) return false;
        unsigned char tail[28] = {0};
        return ::pread(fd, tail, sizeof(tail), st.st_size - 28) != static_cast<ssize_t>(sizeof(tail)) ||
               std::memcmp(tail, kEof, sizeof(kEof)) != 0;
    }

    std::string path_;
    std::FILE* fp_ = nullptr;
    bool gz_ = false;
    bool done_ = false;
    bool zs_open_ = false;
    z_stream zs_{};
    unsigned char in_[1 << 16];
    std::size_t in_pos_ = 0;
    std::size_t in_len_ = 0;
    bool in_eof_ = false;
    std::string buf_;
    std::size_t pos_ = 0;
};

} // namespace panvar
