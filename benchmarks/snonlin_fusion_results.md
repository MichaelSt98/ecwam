# SNONLIN CPU interaction-loop fusion

## Change

Fuse interaction calculation and ordered updates in all three CPU spectral
branches. Keep the original split loops for `WAM_GPU` / `_OPENACC`. Four
horizontal scratch arrays become scalar intermediates on the CPU. Every update
retains its original arithmetic and destination order, including coincident
destinations. Public interfaces and floating-point compiler flags are unchanged.

## Matched O1280 measurement

Slurm job **30299213**, 2026-10-20, completed 45 model runs in 1h38m13s. This was
a nine-condition physics attribution matrix, with one warm-up per condition and
four interleaved, rotated/reversed timed rounds. The fusion contrast differs
only in `snonlin.F90`; both libraries retain the original SINPUT and PROPAGS2.

- AMD EPYC 7742, four nodes, 64 MPI ranks, eight OpenMP threads per rank.
- Intel classic 2021.4, single precision, AVX2, `-O2`, `-fp-model precise`,
  `-fp-speculation=safe`, `-fast-transcendentals`, `-fpe0`, `-no-fma`.
- O1280, NPROMA 32, 36 directions / 29 frequencies, ETOPO1 and subgrid.
- Six-hour integration, 2022-12-31 12:00–18:00; 900 s physics, 450 s propagation,
  225 s fast propagation through frequency 7. Initial/final SWH GRIB output;
  automatic end-date binary restart suppressed equally.
- DR_HOOK profiling identical; no fine-grained propagation phase hooks.
- `OMP_PLACES=cores`, `OMP_PROC_BIND=close`, `OMP_DYNAMIC=FALSE`; identical
  Slurm core binding, no explicit NUMA memory binding in either condition.

| Metric: median of four runs | Split baseline | CPU fusion |
|---|---:|---:|
| Mean-thread accumulated SNONLIN self time | 2.703952 s | 2.193033 s |
| Maximum-thread accumulated SNONLIN self time | 3.236500 s | 2.518500 s |
| Model wall time, including initialisation and configured I/O | 126 s | 122 s |

Kernel mean-time reduction by ratio of medians: **18.90%**. The four within-round
reductions were **17.74%, 18.98%, 19.58%, 19.26%**. A separate matched contrast
with vector-math compilation also improved in every round (16.47–18.13%).

Model times include substantial propagation/exchange variability: the 126/122 s
medians are observations, not an attributed forecast-wide fusion speedup.
Accumulated thread times are not synchronized critical-path measurements.

Artifacts:

```text
/ec/res4/scratch/nams/ecwam-benchmarks/physics-series-20261020/
  timing-summary.json
  repeat-{1,2,3,4}-{baseline,fusion-strict}/logs/model/
```

Frozen sources, compiler caches, library hashes and reports are in
`tmp/performance-series-20261020/variants/{baseline,fusion-strict}`. The
experiment starts from tracked revision `703abb40`, with identical IFS-like Atos
toolchain alignment and compiler-report settings in both conditions. Those
common settings are not part of the measured fusion effect.

## Numerical and regression checks

- All eight O1280 baseline/fusion statistics files are byte-identical, SHA-256
  `3387befece6efb7d0bc64f2a723a6d79591aa66fd8756577336d7b5e36b25b38`.
- O48 with depth/current refraction: all **3,381,264** unquantised restart
  spectrum values at hour 12 are bit-identical. All **68** requested fields
  captured before GRIB packing at hours 12 and 18 are bit-identical.
- The same full-spectrum comparison passes bitwise against a matched
  vector-math baseline. This is distinct from comparing vector math with strict
  scalar math, which can change results.
- Baseline: all 30 model/setup CTests passed. Fusion: all 18 model tests passed,
  with 12 existing setup fixtures skipped; includes O96 GRIB restart and O48
  serial, MPI4, OMP4 and source-disabled coverage.
- The dedicated `ecwam_test_snonlin_fusion_sp` compares complete source and
  derivative arrays in 15 cases: short/non-unit horizontal ranges, interior and
  boundary frequencies, coincident angular destinations and coincident spectral
  destinations. It compiles the CPU and retained split formulations with
  matched precise/no-FMA settings and bounds checking. The test is enabled for
  Intel classic Bit CPU SP builds without Loki.
- The historical synthetic fixture had an out-of-range input neighbor at its
  highest synthetic frequency. The new fixture bounds every input read while
  retaining guarded out-of-domain output destinations; the corrected cases pass
  under both strict and vector-math compilation. This replaces that part of the
  historical fixture evidence.
- fortls reports unavailable-module hints for the fixture; actual compilation
  and execution pass. The new fixture passes the Fortran quick checker.

Run the registered fixture after sourcing the bundle environment:

```sh
ctest --output-on-failure -R ecwam_test_snonlin_fusion_sp
```

For separate strict and vector-math comparisons using an existing Intel classic
SP bundle build (the output directory must not exist):

```sh
bash helpers/ecwam-check-snonlin-fusion.sh BUILD_DIRECTORY NEW_OUTPUT_DIRECTORY
```

The retained split path is exercised as CPU reference code; these checks do not
execute on a GPU. Bitwise evidence above concerns the tested compiler, precision,
flags and inputs.
