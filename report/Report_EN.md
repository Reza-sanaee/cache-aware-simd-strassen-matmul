# Cache-Aware SIMD and Strassen Matrix Multiplication

## Objective and computational problem

This project computes the product of two square matrices of double-precision values and compares a direct reference implementation with a hybrid implementation that combines memory-layout transformation, a four-by-four SIMD microkernel, and Strassen recursion. The mathematical contract is C[i,j] = sum over k of A[i,k] times B[k,j], for i,j,k from zero to n-1. The inputs are read-only and the output occupies separate storage. Each matrix is row-major, and an explicit leading dimension specifies the distance between successive rows. A logical n-by-n submatrix can therefore be processed within a larger parent allocation without copying its rows into a separate matrix.

The direct path preserves the simple pointer-based multiplication used as a numerical reference. The optimized path first makes the inner products contiguous, then evaluates multiple output elements together. Strassen reduces the number of recursive products from eight to seven when a matrix is split into four equal quadrants. The crossover between arithmetic savings and preparation overhead is controlled by a leaf threshold. The default threshold is 512, matching the intended hybrid structure. Dimensions at or below that threshold use the transposed leaf kernel directly. Odd dimensions also use that kernel because an equal-quadrant split is not defined for them.

## Row-major layout, strides, and the reference algorithm

For a row-major matrix A with leading dimension s_A, element A[i,k] is stored at offset i*s_A+k. In the reference i-j-k loop, A's row is read sequentially, while B[k,j] advances by s_B doubles at each inner-loop iteration. C[i,j] is written only after its dot product has been accumulated. The algorithm performs n cubed scalar multiplications and approximately n cubed additions. Its arithmetic complexity is cubic even when pointer arithmetic replaces repeated array indexing.

The strided traversal of B can reduce spatial locality: consecutive values used by one dot product may reside in different cache lines, and a large working set can repeatedly displace useful data. The source-level use of register variables does not guarantee a particular hardware-register allocation. The compiler determines that allocation. This project consequently identifies the memory access pattern and arithmetic scheduling explicitly instead of treating a language keyword as a measured optimization.

The public reference function accepts separate s_A, s_B and s_C values. This matters for the validation experiments, which deliberately use different padded row lengths. A correct result must be stored only in the first n entries of each output row; the padding after those entries is left intact. The same contract applies to the recursive and SIMD paths.

## Transposing the second operand

The leaf routine allocates a packed n-by-n buffer B_T and stores B_T[j*n+k] = B[k*s_B+j]. An inner product then reads both A[i*s_A+k] and B_T[j*n+k] sequentially. Transposition requires quadratic preparation work and an additional n squared doubles. The buffer is released after the leaf product, including the scalar fallback path.

The allocation has 64-byte alignment. Row starts and submatrix offsets are not assumed to retain that alignment. Unaligned AVX loads and stores are therefore used for vector operations, which makes the kernel valid for padded strides and submatrix views. Alignment of the allocation and alignment of every accessed vector are separate properties.

The timing of the optimized function includes this transpose, its allocation and deallocation, arithmetic, and any recursive temporary storage. It is an end-to-end multiplication time rather than a measurement of arithmetic instructions alone. For small matrices the preparation cost can consume much of the arithmetic benefit. No hardware cache-miss counters are collected, so locality is explained from the access pattern and not presented as a measured cache-miss reduction.

## Four-by-four AVX leaf microkernel

An output tile contains up to four rows and four columns. For each group of four k values, the kernel loads one four-double vector from each participating A row and one from each participating B_T row. Every A vector is reused against the B_T vectors for that tile, and every B_T vector is reused against the A vectors. Sixteen independent vector accumulators represent a complete four-by-four tile. The operations are vector multiplication followed by vector addition; fused multiply-add is not required.

After the vector loop, each accumulator's four lanes are summed to produce one output scalar. Remaining k values are accumulated scalarly. Partial tiles at the last rows or columns are also handled explicitly. Thus n need not be divisible by four. A one-by-one matrix, a three-by-three matrix and a matrix with an odd leading dimension obey the same indexing contract as a large matrix.

The AVX routine is isolated in a target-specific function when compiled with GCC on x86. Runtime CPU dispatch selects it only when AVX is available. Otherwise the routine uses the same transposed storage with scalar dot products. Defining MATRIX_DISABLE_AVX builds the scalar path alone. The reference, allocation, recursion, comparison and command-line code do not require a global -mavx compiler flag. The distributed test suite exercises both the AVX-capable build and the scalar-only build.

## Strassen partitioning and seven recursive products

For an even dimension n above the leaf threshold, m=n/2. A and B are partitioned into A11, A12, A21, A22 and B11, B12, B21, B22. Each quadrant is a view into its parent allocation, with its parent's leading dimension. The seven products follow the standard Strassen identities shown below. Temporary sums and differences are packed m-by-m matrices. The output products M1 through M7 are also packed, while the final C quadrants use the caller's output stride.

The products are evaluated sequentially. Two temporary matrices hold sums and differences, and their contents are reused once the corresponding recursive multiplication has completed. Seven product matrices retain the intermediate results needed for final recombination. Additions and subtractions use both source strides explicitly, so a parent view and a packed temporary need not share a row length.

Every recursive call receives its logical dimension and all three leading dimensions separately. The logical dimension is not inferred from an output stride. A nonzero cutoff is required. Even dimensions greater than the cutoff are split recursively; a small or odd subproblem is computed by the complete leaf kernel. This avoids discarding an unmatched last row or column. It also means that the familiar Strassen asymptotic bound applies to the repeatedly splittable case, while an odd problem may use cubic leaf work.

## Recombination, storage, and complexity

The output quadrants are obtained from the seven products using the recombination equations below. Each element is written at the correct offset in the caller's C allocation. The implementation does not overwrite A or B. The input and output buffers must not overlap because the computation reads input values throughout the multiplication.

For repeated even splits, the recurrence is T(n)=7*T(n/2)+O(n squared), which gives O(n to the power log2(7)), approximately O(n to the power 2.8074), rather than cubic arithmetic. This is an arithmetic-count statement, not a guarantee of lower elapsed time for every n. Transposition, matrix addition, allocation and memory traffic influence the crossover in practice.

Nine m-by-m temporary arrays are held at one recursive level: seven products and two work matrices. A descending sequential call adds further temporary levels, and its leaf adds one transpose buffer. For n=1024 with cutoff 512, the top-level temporaries occupy 18 MiB and a 512-by-512 transpose occupies 2 MiB. The four main 1024-by-1024 double arrays occupy 32 MiB, giving approximately 52 MiB of matrix storage at the peak, excluding allocator and runtime overhead. No measured resident-memory trace is claimed.

An allocation failure propagates to the caller. Allocations already obtained by that call are released before returning failure. Matrix dimensions and leading dimensions are checked before entering the kernel. Allocation-size calculations use size_t and guard multiplication overflow. The benchmark accepts positive dimensions and cutoffs up to 8192 and reports allocation failure instead of continuing with a null pointer.

## Floating-point agreement and numerical validation

Floating-point products can differ because SIMD lane reduction and Strassen recombination change the order of arithmetic. Exact bitwise equality is therefore not the numerical contract. Each result element must satisfy an absolute-plus-relative error bound: absolute error no greater than atol + rtol*max(abs(reference),abs(candidate)). Nonfinite values cause comparison failure. The rule remains defined when the reference value is zero, unlike division by the reference value.

The deterministic test suite uses dimensions 1,2,3,4,5,7,8,12,16,31,32,33 and 64. Each dimension is tested with four patterns: a zero left operand, an identity left operand, an identity right operand, and signed nontrivial operands. A cutoff of four exercises multiple recursive levels. Leading dimensions n+3, n+5 and n+7 test independent row strides. Output padding is initialized with a sentinel and checked after every multiplication. This yields 52 matrix cases per backend. Both AVX-capable and scalar-only builds pass all 52 cases, as well as explicit zero and nonfinite comparison checks. The tests use atol=rtol=1e-10 for their bounded signed inputs.

The measured benchmark uses larger positive rational inputs and compares every element against the direct reference with atol=1e-8 and rtol=1e-10. All fifteen benchmark executions pass that criterion. The maximum absolute error is also recorded; its magnitude alone is not the pass criterion because the result values vary in scale.

## Benchmark procedure and measured results

The benchmark is single-threaded. It was compiled with GCC at -O3 and run on Windows on an Intel Core i7-13650HX with AVX enabled. Inputs are generated deterministically from seed 42. The allocation and initialization of the four main matrices occur before timing; the optimized function's internal allocations and transposition occur inside its timed interval. Windows QueryPerformanceCounter provides elapsed wall-clock time. Each dimension is run three times, always with cutoff 512. The baseline is timed before the optimized path in each process.

The table reports the median reference time, median optimized time, minimum-to-maximum optimized time and median of the three per-run speedup ratios. The speedup column is not calculated by dividing the two separately chosen medians. The underlying CSV retains every repetition and its numerical check. Dimensions 64 through 512 use the AVX/transposed leaf directly; n=1024 uses one Strassen split into seven 512-dimensional leaf products.

These timings describe this machine and executable. The three repetitions are a descriptive sample, not a confidence interval or a general performance guarantee. The large spread at n=1024 and the first-run difference at n=256 show the influence of execution conditions. The comparison does not isolate the independent contribution of transposition, AVX tiling and Strassen, nor compare against a tuned BLAS library. Such ablations and additional machines would be needed for broader architectural conclusions.

## Build, execution, and interpretation

The code package separates the matrix library, command-line benchmark, deterministic tests and measured CSV. On Windows, build.ps1 compiles the benchmark and tests with MinGW GCC and executes the tests. On a POSIX system, make and make test provide the corresponding build. The benchmark accepts dimension, leaf cutoff and seed as positional arguments. A run such as matrix_benchmark 1024 512 42 emits a CSV header and one measurement row. Every measurement includes the AVX-dispatch state, both elapsed times, speedup, maximum absolute error and verification flag.

The resulting implementation demonstrates memory-layout transformation, tiled vector arithmetic, recursive matrix algorithms and numerical validation in one program. Its conclusions are limited to measured end-to-end timing and verified numerical agreement for the stated inputs. It does not measure power, prove a universal cache policy improvement or establish optimal performance across processors.
