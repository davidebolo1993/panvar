#include "panvar/gz_inflate.hpp"

#include <stdexcept>
#include <string>

#ifdef PANVAR_ZLIB_NG
#include <zlib-ng.h>
using Stream = zng_stream;
#define PV_Z(f) zng_##f
#define PV_ZLIB_NAME "zlib-ng " ZLIBNG_VERSION
#else
#include <zlib.h>
using Stream = z_stream;
#define PV_Z(f) f
#define PV_ZLIB_NAME "zlib " ZLIB_VERSION
#endif

namespace panvar {

Inflater::Inflater() : stream_(new Stream{}) {
    if (PV_Z(inflateInit2)(static_cast<Stream*>(stream_), 15 + 16) != Z_OK) {   // gzip wrapper
        delete static_cast<Stream*>(stream_);
        throw std::runtime_error("cannot initialise decompression");
    }
}

Inflater::~Inflater() {
    PV_Z(inflateEnd)(static_cast<Stream*>(stream_));
    delete static_cast<Stream*>(stream_);
}

Inflater::Status Inflater::run(const unsigned char* in, std::size_t in_len, char* out, std::size_t cap, std::size_t& used,
                               std::size_t& produced) {
    Stream* zs = static_cast<Stream*>(stream_);
    zs->next_in = const_cast<unsigned char*>(in);
    zs->avail_in = static_cast<unsigned>(in_len);
    zs->next_out = reinterpret_cast<unsigned char*>(out);
    zs->avail_out = static_cast<unsigned>(cap);
    const int ret = PV_Z(inflate)(zs, Z_NO_FLUSH);
    used = in_len - zs->avail_in;
    produced = cap - zs->avail_out;
    if (ret == Z_STREAM_END) return kMemberEnd;
    if (ret == Z_OK || ret == Z_BUF_ERROR) return produced > 0 ? kOk : kNeedInput;
    return kError;
}

void Inflater::next_member() { PV_Z(inflateReset)(static_cast<Stream*>(stream_)); }

const char* Inflater::message() const { return static_cast<const Stream*>(stream_)->msg; }

const char* Inflater::library() { return PV_ZLIB_NAME; }

} // namespace panvar
