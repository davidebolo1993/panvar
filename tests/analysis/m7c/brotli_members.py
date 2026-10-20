"""Decode a file of concatenated Brotli streams (large-window enabled) and require complete input consumption.

Adapted, after review, from Sol's independent audit (audit_streams.py): libbrotlidec through ctypes, decoder parameter 1
(BROTLI_DECODER_PARAM_LARGE_WINDOW) set to 1, one decoder instance per member, result codes 1 = success (member complete),
3 = needs more output; anything else, or input left over without progress, raises.
"""
import ctypes as C

LIB = "/Users/davide.bolognini/miniconda3/lib/libbrotlidec.1.2.0.dylib"
_lib = C.CDLL(LIB)
U8 = C.c_uint8
PU8 = C.POINTER(U8)
_lib.BrotliDecoderCreateInstance.argtypes = [C.c_void_p] * 3
_lib.BrotliDecoderCreateInstance.restype = C.c_void_p
_lib.BrotliDecoderDestroyInstance.argtypes = [C.c_void_p]
_lib.BrotliDecoderSetParameter.argtypes = [C.c_void_p, C.c_int, C.c_uint32]
_lib.BrotliDecoderSetParameter.restype = C.c_int
_lib.BrotliDecoderDecompressStream.argtypes = [C.c_void_p, C.POINTER(C.c_size_t), C.POINTER(PU8), C.POINTER(C.c_size_t), C.POINTER(PU8), C.POINTER(C.c_size_t)]
_lib.BrotliDecoderDecompressStream.restype = C.c_int
_lib.BrotliDecoderGetErrorCode.argtypes = [C.c_void_p]
_lib.BrotliDecoderGetErrorCode.restype = C.c_int


def decode_members(data):
    """Returns (text, [(consumed_bytes, produced_bytes) per member]); raises on any invalid or truncated member."""
    source = (U8 * len(data)).from_buffer_copy(data)
    remain = C.c_size_t(len(data))
    ip = C.cast(source, PU8)
    chunks, members = [], []
    while remain.value:
        state = _lib.BrotliDecoderCreateInstance(None, None, None)
        assert state
        assert _lib.BrotliDecoderSetParameter(state, 1, 1)
        before = remain.value
        out = bytearray()
        try:
            while True:
                buf = (U8 * 65536)()
                avail = C.c_size_t(len(buf))
                op = C.cast(buf, PU8)
                total = C.c_size_t()
                status = _lib.BrotliDecoderDecompressStream(state, C.byref(remain), C.byref(ip), C.byref(avail), C.byref(op), C.byref(total))
                out.extend(bytes(buf[:len(buf) - avail.value]))
                if status == 1:
                    break
                if status != 3:
                    raise ValueError(f"invalid or truncated member {len(members)}: status {status}, error {_lib.BrotliDecoderGetErrorCode(state)}, "
                                     f"remaining {remain.value} bytes")
            if remain.value >= before:
                raise ValueError("no progress")
        finally:
            _lib.BrotliDecoderDestroyInstance(state)
        members.append((before - remain.value, len(out)))
        chunks.append(out.decode("utf8"))
    assert sum(n for n, _ in members) == len(data), "input not completely consumed"
    return "".join(chunks), members
