# Measurement environment

- Date: 2026-10-12.
- Processor: 13th Gen Intel Core i7-13650HX.
- Operating system: Windows, x86-64.
- Compiler: MSYS2 MinGW GCC 13.2.0 (Rev4).
- Benchmark flags: -O3 -std=c11 -Wall -Wextra -Werror.
- SIMD: AVX selected at runtime; separate scalar-only tests also pass.
- Workload: dimensions 64,128,256,512,1024; seed 42; leaf cutoff 512; three runs per dimension.
- Timer: QueryPerformanceCounter wall-clock time.
- Timing includes internal transpose, temporary allocation, recursion and arithmetic, and excludes initial main-matrix allocation and filling.
- Benchmark tolerance: absolute 1e-8 plus relative 1e-10 times the larger result magnitude.
- Test tolerance: absolute 1e-10 plus relative 1e-10; 52 patterned, padded-stride cases per backend, plus zero/nonfinite comparison checks.
