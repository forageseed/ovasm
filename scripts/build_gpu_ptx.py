"""Compile src/gpu_delta.cu to src/gpu_delta.ptx with NVRTC (no nvcc needed).

usage: python3 scripts/build_gpu_ptx.py PATH/TO/libnvrtc.so.12
(e.g. from the nvidia-cuda-nvrtc-cu12 wheel). The PTX targets compute_75 and is compiled for the
actual GPU by the driver when ovasm loads it, so building ovasm itself needs no CUDA.
"""
import ctypes
import pathlib
import sys

src = pathlib.Path(__file__).resolve().parent.parent / "src"
nvrtc = ctypes.CDLL(sys.argv[1])
code = (src / "gpu_delta.cu").read_bytes()
prog = ctypes.c_void_p()
assert nvrtc.nvrtcCreateProgram(ctypes.byref(prog), code, b"gpu_delta.cu", 0, None, None) == 0
opts = [b"--gpu-architecture=compute_75", b"--std=c++14"]
arr = (ctypes.c_char_p * len(opts))(*opts)
status = nvrtc.nvrtcCompileProgram(prog, len(opts), arr)
size = ctypes.c_size_t()
nvrtc.nvrtcGetProgramLogSize(prog, ctypes.byref(size))
log = ctypes.create_string_buffer(size.value)
nvrtc.nvrtcGetProgramLog(prog, log)
if status != 0:
    sys.exit(log.value.decode())
nvrtc.nvrtcGetPTXSize(prog, ctypes.byref(size))
ptx = ctypes.create_string_buffer(size.value)
nvrtc.nvrtcGetPTX(prog, ptx)
(src / "gpu_delta.ptx").write_bytes(ptx.value)
print(f"wrote {src / 'gpu_delta.ptx'} ({len(ptx.value)} bytes)")
