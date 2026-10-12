$ErrorActionPreference = 'Stop'
$projectRoot = $PSScriptRoot
$outputDir = Join-Path $projectRoot 'bin'
New-Item -ItemType Directory -Path $outputDir -Force | Out-Null
& gcc -O3 -std=c11 -Wall -Wextra -Werror -I (Join-Path $projectRoot 'src') (Join-Path $projectRoot 'src/matrix.c') (Join-Path $projectRoot 'src/main.c') -o (Join-Path $outputDir 'matrix_benchmark.exe') -lm
if ($LASTEXITCODE -ne 0) { throw 'Benchmark compilation failed' }
& gcc -O2 -std=c11 -Wall -Wextra -Werror -I (Join-Path $projectRoot 'src') (Join-Path $projectRoot 'src/matrix.c') (Join-Path $projectRoot 'tests/test_matrix.c') -o (Join-Path $outputDir 'test_matrix.exe') -lm
if ($LASTEXITCODE -ne 0) { throw 'Test compilation failed' }
& (Join-Path $outputDir 'test_matrix.exe')
if ($LASTEXITCODE -ne 0) { throw 'Matrix tests failed' }
& gcc -O2 -std=c11 -Wall -Wextra -Werror -DMATRIX_DISABLE_AVX -I (Join-Path $projectRoot 'src') (Join-Path $projectRoot 'src/matrix.c') (Join-Path $projectRoot 'tests/test_matrix.c') -o (Join-Path $outputDir 'test_matrix_scalar.exe') -lm
if ($LASTEXITCODE -ne 0) { throw 'Scalar test compilation failed' }
& (Join-Path $outputDir 'test_matrix_scalar.exe')
if ($LASTEXITCODE -ne 0) { throw 'Scalar matrix tests failed' }
