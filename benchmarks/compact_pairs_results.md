# Opt-in compact propagation coefficients

`ENABLE_PROPAG_FACTORED_PAIRS=ON` enables CPU SP latitude/corner coefficient
pairs. It defaults OFF and excludes OpenACC/Loki. The SP target expression keeps
the active DP propagation implementation unchanged in a mixed build. No public
Fortran entry-point signature is changed.

The default preprocessed PROPAGS2 body is retained. The opt-in refraction path
uses the tested paired traversal and ordered spectral terms, with expanded
fallback for direct callers without a prepared private plan. Tiling, native
output, persistent workspaces, ownership changes and relaxed propagation FP
flags are not part of this option.

## Storage and lifecycle

Expanded weights remain authoritative through generation, CFL retries,
obstruction application and activity-flag construction. Only successful final
weights are converted to private latitude/corner bases. Their expanded arrays
are then deallocated. Existing geometry reconstructs the two incoming legs;
`SUMWN` remains independent outgoing loss and is not reconstructed.

Private bases are invalidated before rebuilding weights and released at
WVDEALLOC. Weight/dot invalidation ensures correct reinitialisation. This uses
the repository's legacy YOW-module private-storage ownership pattern; no private
state is added to the public API argument lists.

Allocated coefficient storage falls from 21 to 15 field equivalents, including
unused zero-index planes in both counts: **28.57%**, approximately
**6.456 → 4.612 GB/rank**. This is allocation size, not measured DRAM traffic.

## Observed live memory and peak trade-off

Job **30375110** compares expanded/paired libraries with identical local binding
on O1280, four Zen 2 nodes, 64 ranks, eight threads/rank. A validated diagnostic
wrapper samples `/proc/self/status` after known completed propagation calls,
before final deallocation. All 64 ranks contribute.

| Mean process memory | Expanded | Paired |
|---|---:|---:|
| RSS after propagation call 2 | 7.361267 GiB | 5.642888 GiB |
| RSS after call 48 | 7.353611 GiB | 5.637191 GiB |
| Whole-process peak RSS | 7.981904 GiB | 9.404983 GiB |

Live RSS falls about **1.72 GiB/rank / 23.3%**, consistent with the coefficient
storage saving. **Preparation raises mean peak RSS by 1.42 GiB/rank / 17.8%**
because old expanded weights and new bases temporarily overlap. Maximum process
peaks are 8.15/9.57 GiB. This option reduces the live working set, not the memory
capacity needed to survive weight preparation.

## Timing evidence and scope

The broad verified matrix, job **30352747**, finds that strict binding plus
pairs improves complete propagation against its strict baseline in all four
rounds. Pair increments over binding alone are +0.63%, +2.53%, +10.47%, +2.38%.
Those prototypes used a common refactored baseline.

Job **30398792** rechecks the cleaned production extraction, with original
default propagation, committed SNONLIN fusion, strict FP and local binding held
fixed. It runs two warm-ups and four alternating matched pairs:

| Median of four runs | Pairs OFF | Pairs ON |
|---|---:|---:|
| Complete propagation, mean accumulated PE time | 30.712500 s | 30.632031 s |
| PROPAGS2 mean-thread self time | 13.681597 s | 12.054179 s |
| Model wall time | 111 s | 110 s |

The kernel median is about **11.9% lower**, but complete propagation does not
improve repeatably: paired reductions are +16.79%, −0.75%, −1.29%, −15.62%.
The cleaned option is therefore selected as a **live-memory/working-set
optimisation**, without a repeatable complete-propagation speedup claim. User
acceptance permits a useful memory saving with neutral overall performance.

## Numerical checks

Coefficient reconstruction and latitude/corner accumulation order change
rounding. The option is not bit-identical:

- The dedicated full-output fixture checks 4,431,840 actual paired outputs over
  eight cases and three extra fast steps, including flags/fallbacks, single,
  empty and partial ranges, geometry endpoints, obstruction, sentinel reads,
  unchanged output regions and frozen slow input. Its pointwise tolerance check
  passes. Reference dump files contain expanded outputs; actual paired outputs
  are compared in memory.
- All **33 CTests pass** with the cleaned option ON, including serial/MPI/OMP,
  source-disabled, O96 GRIB restart, rank-log retention and the paired fixture.
  The final fixture also passes after replacing implicit reference-array
  allocation with explicit allocation/bounded copies.
- O48 hour-12 full spectrum: max absolute error 3.0441284e-3, scaled
  3.3734661e-4, RMS 2.9816173e-6. OFF is bitwise with the accepted strict baseline.
- O1280 18-hour job **30399910** compares the cleaned bound option against the
  bound strict accepted sources, with SINPUT vector math OFF. All 6,106,530,384
  final spectrum values are checked: max absolute error **0.01880646**, scaled
  **0.00399670**, RMS **1.4946893e-6**. The initial 1e-6 absolute-plus-relative
  screen is exceeded by 3,130,583 values (about 0.051%).
- Raw SWH max differences are **0.01019189 m at hour 12** and **0.00405896 m at
  hour 18**; output-grid RMS errors are **3.0769e-6 / 1.5342e-6 m**, including
  unchanged missing cells. All 136 field/date dumps are compared before GRIB
  packing. These localized differences require product-specific accuracy
  assessment; they are not an all-fields tolerance pass.
- fortls reports expected unavailable-module hints without syntax errors; the
  focused patch passes the Fortran quick checker.

Enable only the desired tested mode on a clean source tree, reconfiguring an
existing bundle build:

```sh
--reconfigure --cmake="ENABLE_PROPAG_FACTORED_PAIRS=ON ENABLE_SINPUT_ARD_VECTOR_MATH=OFF"
```

The two numerical-relaxation options have been validated separately, not as an
automatically combined package. Peak reduction through a streamed producer
without expanded/new overlap remains an unimplemented direction.

Artifacts under `/ec/res4/scratch/nams/ecwam-benchmarks/`:

```text
pairs-package-20261020/timing-summary.json
pairs-memory-20261020/memory-summary.json
pairs-long-20261020/repeat-1-pairs-bound/full-numerical-comparison/summary.json
```
