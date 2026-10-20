# Opt-in SINPUT_ARD vector math

## Configuration

`ENABLE_SINPUT_ARD_VECTOR_MATH=ON` appends `-fpe3` to **SINPUT_ARD in the
single-precision library only**, for Intel classic Bit CPU builds without Loki. The option
defaults OFF. The source, precise FP model, safe speculation and no-FMA settings
are retained; SNONLIN and other files receive no new FP overrides. Mixed SP/DP
generated commands were checked: SP has the override and DP does not.

Example addition to a default Atos, SP, fckit, Ninja bundle build:

```sh
--reconfigure --cmake="ENABLE_SINPUT_ARD_VECTOR_MATH=ON"
```

This changes compilation-time floating-point exception assumptions and enables
SVML vector math; rounding can change. It does not explicitly disable global
runtime traps. Caller/main-program and DR_HOOK exception settings still apply.
The tested path has `TAUWSHELTER=0` and `ISNONLIN=0`.
For an existing bundle build, `--reconfigure` is necessary to apply a changed
option rather than reuse its cached value; verify the saved CMake cache.

## Measured performance

Job **30299213**, 2026-10-20, compares original-source SINPUT_ARD with/without
its per-file `-fpe3`. The other physics and propagation sources/flags are fixed.
Hardware, inputs, scopes and run design are described in
[`snonlin_fusion_results.md`](snonlin_fusion_results.md): O1280, four Zen 2 nodes,
64 ranks, eight threads/rank, one warm-up/condition, four interleaved timed
rounds, one six-hour interval.

| Metric: median of four runs | Original `-fpe0` compilation | SINPUT-only `-fpe3` |
|---|---:|---:|
| Mean-thread SINPUT_ARD self time | 3.816859 s | 2.079940 s |
| Maximum-thread SINPUT_ARD self time | 4.375000 s | 2.383000 s |
| Model wall time | 126 s | 125.5 s |

Mean-thread kernel reduction by ratio of medians: **45.51%**. Within-round
reductions: **44.96%, 45.76%, 44.90%, 45.25%**. The model result is small and
subject to propagation/arrival variability; the kernel reduction is not a
forecast-wide speedup claim.

The neutral-operand source rewrite was also measured independently. Under
strict flags its kernel median is **56.09% slower** than the original source;
with vector math it is **36.20% slower** than original-source vector math.
The opt-in implementation therefore uses the original guarded source.

## Numerical differences

O48 integrates 18 hours, with 12 directions, 25 propagated/reduced frequencies
and 36 physics frequencies. The comparison reads all 3,381,264 unquantised
restart values at hour 12, and all 68 requested raw fields before GRIB packing
at hours 12 and 18.

| Quantity | Observed difference vs strict original |
|---|---:|
| Spectrum maximum absolute error | 1.21974945e-3 |
| Spectrum maximum scaled error, `abs(error)/max(1,abs(ref))` | 2.09913426e-4 |
| Spectrum RMS error | 1.34466786e-6 |
| SWH maximum absolute error | 1.27005577e-3 m |
| SWH RMS error across the two raw output grids, including unchanged missing cells | 6.49390784e-6 m |
| Largest scaled error among continuous raw fields | 8.69508124e-4 |

These are small deviations in the tested case, **not bitwise equality**. The
initial pointwise `1e-6 + 1e-6*abs(reference)` screen is exceeded by 1,154 spectrum
values and 66 SWH values. Some derived products amplify the last-bit changes;
one discrete event count changes by one. Field-specific accuracy requirements
must therefore be assessed from the actual differences rather than a final
global norm alone. The largest absolute error across fields of different units
is not a useful combined accuracy metric.

Decoded GRIB comparisons can exaggerate differences because changed packing
parameters shift many decoded values. Raw fields were captured immediately
before packing, and the diagnostic capture was verified not to alter any saved
restart spectrum. The numerical results above use raw values.

## Implementation validation and artifacts

- Production option OFF reproduces the strict fused-SNONLIN baseline restart
  bitwise. ON reproduces the isolated, measured SINPUT-vector candidate restart
  bitwise; the CPU fusion has no numerical contribution to those deviations.
- ON passes nine selected model tests: O48 normal/source-disabled and O96 GRIB
  restart, each serial/MPI4/OMP4. Six existing setup fixtures skip. The SNONLIN
  bitwise fixture also passes.
- Mixed SP/DP configuration succeeds and generated SINPUT commands confirm the
  SP-only flag scope. This is a configuration check, not DP numerical testing.
- Tests-disabled configuration succeeds.
- Compiler reports show the appended `-fpe3` only for ON; vector math is also
  supported by the preserved assembly/SVML evidence from the original-source
  experiment.

```text
/ec/res4/scratch/nams/ecwam-benchmarks/physics-series-20261020/timing-summary.json
tmp/performance-series-20261020/numerics-o48/spectrum-comparison.json
tmp/performance-series-20261020/numerics-o48-capture/raw-field-comparison.json
tmp/performance-series-20261020/variants/sinput-option-final-{off,on}/
tmp/performance-series-20261020/mixed-config/sinput-precision-audit.json
```

## Larger-domain accuracy qualification

The clean accepted-tree O1280 18-hour run, job **30372070**, gives a stronger
accuracy check than the O48 screen. It compares all 6,106,530,384 final spectrum
values and raw fields at hours 12 and 18. The SINPUT-only vector mode has:

- Spectrum max absolute error **0.12002325**, max scaled error **0.02186057**,
  RMS error **4.85055483e-6**; 3,002,553 values exceed the initial pointwise
  absolute-plus-relative screen, about 0.049% of the spectrum.
- Raw SWH max difference **0.01158351 m at hour 12** and **0.01492977 m at hour
  18**. RMS differences over the output grids are **3.4768e-6 / 4.5156e-6 m**,
  including unchanged missing-grid cells.
- 124 of 136 field/date dumps differ. Bitwise binding comparisons confirm that
  these changes come from the vector-math mode, not NUMA placement.

Thus the option is **not an all-fields tolerance pass**: global RMS differences
are small, but localized spectral-bin changes reach 2.19% under the stated
scaled metric, and SWH changes reach 1.5 cm. It remains default-off and requires
case/product-specific accuracy assessment. The bitwise-reproducible selected
configuration leaves it OFF.

Both strict and vector modes reproduce their own uninterrupted O48 hour-18
spectrum bitwise after a binary restart at hour 12. The vector-vs-strict hour-18
O48 spectrum max scaled difference is 2.0522e-3; the earlier table explicitly
describes the hour-12 spectrum.

Full-size results: `accepted-series-long-20261020/large-validation-summary.json`
under the scratch benchmark base. These diagnostic runs include large raw field
output and a complete binary checkpoint; their single wall-time observations
are not a repeated production speedup estimate.
