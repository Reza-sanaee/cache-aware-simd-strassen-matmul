# Build and execution

## Contents

- `src/matrix.h`, `src/matrix.c`: stride-aware multiplication library, SIMD leaf and Strassen recursion.
- `src/main.c`: reproducible command-line benchmark.
- `tests/test_matrix.c`: numerical and output-padding checks.
- `bin/`: executables built on Windows x86-64 with MinGW GCC.
- `results/benchmark.csv`: all measured repetitions; `results/ENVIRONMENT.md` gives the measurement settings.

## Windows

Install MinGW GCC and put `gcc` on PATH. From this directory:

```powershell
./build.ps1
./bin/matrix_benchmark.exe 1024 512 42
```

The arguments are matrix dimension, Strassen leaf cutoff and nonzero seed. Defaults are 128,512,42. Output is CSV and exit code zero means numerical agreement. No interactive keyboard input is required. Rebuild on a different target platform instead of using the Windows binaries.

## POSIX

```sh
make
make test
./bin/matrix_benchmark 1024 512 42
```

To exercise scalar-only dispatch, add `-DMATRIX_DISABLE_AVX` to the compiler flags. No global `-mavx` flag is required. Library callers must supply separate input/output storage, valid leading dimensions and a positive cutoff. Callers own the main matrices; internal temporary buffers are allocated and freed by the library.
