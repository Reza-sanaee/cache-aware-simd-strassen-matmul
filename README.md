# Cache-Aware SIMD and Strassen Matrix Multiplication

A double-precision C library combines transposed storage, a 4-by-4 AVX leaf, Strassen recursion, scalar fallback, and explicit submatrix strides.

## Technical highlights

- Implemented tail-safe AVX and scalar kernels with fixed-seed benchmarks and numerical validation.
- Validated 52 patterned and padded-stride cases per backend; recorded an 18.6x median speedup at n=1024 versus a direct i-j-k reference on an Intel Core i7-13650HX across three runs. This is not an optimized BLAS comparison.

## Explore the project

- [Source and execution guide](code/RUNNING.md)
- [Technical report (PDF)](report/Report_EN.pdf)
- [Technical report (Markdown)](report/Report_EN.md)
- [Assignment/specification (PDF)](assignment/Assignment_EN.pdf)

## Reproduction and attribution

Read the execution guide and report before running the code. Dependencies, input data, licensed tools, and hardware requirements vary by project. Recorded results are distinguished from rerun results in the accompanying documentation. Team and third-party source attribution in the reports and source headers is retained. This portfolio does not claim independent authorship of supplied frameworks or dependencies.
