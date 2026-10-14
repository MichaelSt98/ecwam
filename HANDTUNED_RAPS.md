# RAPS ecWAM handtuned SP integration (preparation)

This checkout starts from RAPS ecWAM 1.5.13 at
`6ef5609d74094fc806d0aaf06d5e8d3ae5063253`. The checked-in GPU
implementation in `handtuned/` was copied from Lukas's reproduced revision
`283d1e6273ed0de5481ed5928ac57bd158f0eeed` at
`implementations/handtuned/source/`. The latter's standalone CMakeLists.txt
is included for provenance but is not used by this build.

`ECWAM_HANDTUNED_SP=ON` selects the checked-in handtuned Fortran sources
instead of the usual Loki-transformed ecWAM sources. The RAPS project still
generates the `ecwam_intfb` interfaces and exports `ecwam_sp` and the usual
ecWAM package name, program targets and coupling dependencies. UNWAM sources
come from the RAPS checkout. FCM scans `handtuned/` and the RAPS `unwam/`
sources, avoiding duplicate routine interfaces with the baseline sources.
The switch is opt-in; the existing source selection is the default.

## Build contract

- Single precision, NVIDIA OpenACC, separate memory with CUDA; no Loki
  transformation of ecWAM (`ECWAM_ENABLE_LOKI=OFF`). Other IFS components may
  still use Loki.
- NVHPC 26.9 was used for the reproduced handtuned result and the preparation
  build. NVHPC 25.9 rejects the checked-in `minctapersm` directives. The
  CMake switch requires at least 26.5; that lower bound is a syntax check,
  not a validated performance configuration.
- Pass a numeric `CMAKE_CUDA_ARCHITECTURES` (90 for the reproduced SM90 case).
  Keep `FIELD_API_ENABLE_ACC=ON`, `FIELD_API_ENABLE_CUDA=ON`, and a matching
  compiler/CUDA/MPI dependency build. Do not mix 25.9 Fortran modules with a
  26.9 compilation or use the unified-memory field_api prefix.
- `ECWAM_ENABLE_UNWAM=ON` is supported in the standalone preparation build.
  IFS ocean coupling, ecflow and the full `wamassi` path still need coupled
  build and runtime validation.

For example, in a **fresh** IFS bundle configuration based on `50r3_opt.yml`,
set the ecwam project to `dir: /path/to/ecwam-raps-handtuned-sp` instead of
`git`/`version`, and arrange for these options at bundle configure time:

```
ECWAM_HANDTUNED_SP=ON
ECWAM_ENABLE_SINGLE_PRECISION=ON
ECWAM_ENABLE_DOUBLE_PRECISION=OFF
ECWAM_ENABLE_LOKI=OFF
ECWAM_ENABLE_ACC=ON
ECWAM_ENABLE_CUDA=ON
CMAKE_CUDA_ARCHITECTURES=90
```

RAPS's `--with-ecwam-gpu` and `--with-gpu` options currently set
`ECWAM_ENABLE_LOKI=ON`; override it explicitly for the handtuned variant.
Use a new IFS build directory and dependencies rebuilt with the selected
compiler. Start with forecast-only SP, then check the full `wamassi` and NEMO
configuration, a short coupled forecast and numerical comparison against the
unmodified RAPS build before performance measurements.

## Preparation build result (14 October 2026)

Using the NVHPC 26.9 separate-memory GPU dependency prefix from the private
reproduction, an isolated Release configure with SP/ACC/CUDA/UNWAM ON and
LOKI/OCEAN_COUPLING/ECFLOW OFF built both `libecwam_sp.so` and
`ecwam-chief-sp`. This is a compile/link smoke check; it does not establish
IFS compatibility, full-feature compilation, numerical equivalence or the
standalone reproduced performance. The build is in the ignored
`build-handtuned-26.9/` directory.
