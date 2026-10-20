# Validated performance commit series

## Selected changes

| Change | Evidence | Numerical status / scope |
|---|---|---|
| SNONLIN CPU interaction fusion | 18.90% lower matched kernel median; all four rounds faster | Full O1280 18-hour spectrum and raw fields byte-identical |
| Rank-local NUMA wrapper | 11.75% lower matched complete-propagation median; all four rounds faster | Full O1280 spectrum/raw fields byte-identical; configuration-dependent placement gain |
| SINPUT vector math, opt-in | 45.51% lower matched kernel median; all four rounds faster | Default OFF; localized O1280 differences reach 1.5 cm SWH and 2.19% scaled spectral-bin error |
| Compact coefficients, opt-in | 28.57% lower coefficient allocation; observed live RSS 7.36→5.64 GiB/rank | Default OFF; higher preparation peak 7.98→9.40 GiB/rank; localized raw SWH differences up to 1.02 cm |

The compact representation's cleaned confirmation has about 11.9% lower kernel
median but no repeatable complete-propagation gain. It is a live-memory
optimisation, not a peak-memory or model-speedup claim. Individual percentages
are not added to describe a combined package.

The strict/default implementation changes no FP flags or public entry-point
signatures. The two changed-rounding options have been validated separately;
their simultaneous use is not the measured package.

## Conservative normal-use configuration

Keep numerical-relaxation options OFF for the bitwise-validated source mode:

```text
ENABLE_SINPUT_ARD_VECTOR_MATH=OFF
ENABLE_PROPAG_FACTORED_PAIRS=OFF
```

SNONLIN CPU fusion is active without a new option. The retained `WAM_GPU` /
`_OPENACC` split path remains available. Invoke `helpers/ecwam-numa-local.py`
inside an affinity-setting launcher when choosing the tested local memory
policy; see [placement usage and evidence](numa_binding_results.md).

On existing bundle builds use `--reconfigure` when changing options, and inspect
the saved cache. One preflight demonstrated that ecbundle's `--when-needed`
path can otherwise reuse an old option value despite changed `--cmake` inputs.

## Accuracy assessment

The reference and strict fused source produce the same final binary O1280
spectrum: **6,106,530,384** values in 1,296 records. All 68 raw output fields at
two dates are byte-identical. NUMA placement preserves those outputs in both
strict and vector modes. O48 binary restart continuation at hour 12 reproduces
uninterrupted hour-18 spectra bitwise in reference, strict and vector modes.

The opt-in modes do not pass a universal pointwise `1e-6 + 1e-6*abs(reference)`
screen. Their measured global RMS errors are small, but rare local changes and
discrete-product boundaries matter. Evaluate the actual product requirements
before selecting them; the detailed results give maxima, RMS values and failure
counts rather than claiming blanket numerical equivalence.

- [SNONLIN implementation and validation](snonlin_fusion_results.md)
- [SINPUT vector-math accuracy/performance](sinput_vector_math_results.md)
- [NUMA placement](numa_binding_results.md)
- [Compact coefficients and memory trade-off](compact_pairs_results.md)

## Tested environment and runs

Intel classic 2021.4, Bit, SP, default Atos architecture, AVX2, no-FMA strict
baseline; fckit, Ninja, four Zen 2 nodes, 64 ranks, eight threads/rank, O1280,
NPROMA 32. The configured 29 frequencies are the propagated/reduced frequencies;
physics and binary restart spectra use 36 frequencies and 36 directions.

| Job | Purpose | Status |
|---|---|---|
| 30299213 | Nine-condition physics attribution, 45 model runs | Completed |
| 30312134 | Initial propagation matrix | Cancelled after missing inclusive diagnostics were detected; secondary data retained |
| 30352747 | Verified fourteen-condition propagation/placement matrix, 70 runs | Completed |
| 30372070 | Clean accepted-tree 18-hour full-spectrum/raw-field validation | Completed |
| 30375110 | Known-phase all-rank live/peak memory comparison | Completed |
| 30398792 | Cleaned compact-coefficient four-pair confirmation | Completed |
| 30399910 | Cleaned compact-coefficient 18-hour numerical validation | Completed |

The accepted-tree validation isolates committed sources from the root worktree's
preserved experiments. No result relies on compiler cost estimates alone.
All 33 ecWAM CTests pass for the cleaned coefficient build; dedicated fixtures
cover coincident SNONLIN destinations, partial/empty stencil ranges, paired
weights and extra fast steps. fortls gives unavailable-module hints without
syntax errors; the focused Fortran patches pass the quick checker.

## Timing/data-quality safeguards

Complete propagation includes caller-side automatic allocation and cleanup.
The verified matrix uses an actual outer hook and cross-checks its scope against
SYSTEM_CLOCK around the whole call. Every case must retain diagnostics from all
64 ranks, the expected kernel call count and final-date norms. Accumulated PE
times are not a synchronized global critical path.

A pre-existing GSTATS fallback-unit bug closed non-master logs before final
timers. The separate RUNWAM fix keeps caller-owned units open and has an MPI4
regression test. Missing hooks were also detected before acceptance; fresh
source mtimes and per-variant compile commands prevent the old timestamp issue.

Numerical comparisons use unquantised spectra and fields before GRIB packing.
The raw-capture/memory diagnostic libraries were checked to preserve spectra
and are excluded from the repeated timing matrices. The 18-hour diagnostic
runs include dense raw output and full binary checkpoints, so their single
model-wall observations are not a clean production speedup estimate.

## Reproduction and retained artifacts

Scratch base: `/ec/res4/scratch/nams/ecwam-benchmarks/`. Each named run directory
retains profiles, all-rank logs, effective namelists, loading/placement records
and summaries. `tmp/performance-series-20261020/` retains source exports, build
caches, library/source hashes, reports, standalone fixtures and original failed
preflight artifacts. Investigation scripts named `ecwam-performance-series-*`
are session artifacts; normal-use code is the focused committed implementation.

Build with the repository bundle workflow, default Atos architecture, Ninja and
`--with-fckit`. Source the build environment before running CTest inside its
`ecwam` directory. To stream SP spectrum comparisons with bounded memory:

```sh
cc -std=c11 -O2 helpers/ecwam-compare-sp.c -lm -o compare-sp
./compare-sp REFERENCE_BLS ACTUAL_BLS
python3 tests/ecwam_compare_sp_check.py ./compare-sp NEW_TEST_OUTPUT_DIRECTORY
```

The tool accepts big-endian Fortran sequential SP spectra, including continuation
subrecords, or `--raw-little` for the diagnostic field dumps. It checks finite
values, counts beyond 32 bits, exact bit patterns and the documented numerical
screen. It has signed-zero, odd-subrecord/buffer-boundary, truncation and NaN
checks and agrees with the independent Python O48 reader.

## Remaining questions

- Product-specific acceptance of localized changed-rounding deviations.
- Coupled IFS, other decompositions/hardware, DP and GPU performance are not
  established by this standalone SP study.
- Reduce compact preparation peak by producing bases without expanded/new
  overlap, preserving CFL checks and original activity flags.
- Source-only spectral fusion/unrolling and aggressive propagation FP flags did
  not establish a sufficiently robust independent benefit for this series.
- Exchange overlap and persistent communication buffers remain unimplemented
  leads; allocation pooling, native output and tiling have not earned inclusion.
