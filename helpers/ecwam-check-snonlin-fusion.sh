#!/usr/bin/env bash
# Intel classic SP check of the CPU fusion against the retained split-loop path.
set -eo pipefail
root=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
build=$(realpath "${1:?Pass an Intel classic single-precision bundle build}")
out=${2:?Pass a new output directory with an existing parent}
[[ -d $(dirname "${out}") && ! -e ${out} ]]
source "${build}/env.sh"
set -u
shopt -s inherit_errexit
mkdir "${out}"
out=$(realpath "${out}")
flags=(-O2 -g -traceback -fpp -DWAM_HAVE_SINGLE_PRECISION
  '-D__FILENAME__="snonlin.F90"' -fp-model precise -fp-speculation=safe
  -fast-transcendentals -march=core-avx2 -no-fma -qopenmp -heap-arrays 32
  -check bounds
  "-I${build}/module/ecwam_sp" "-I${build}/module/parkind_sp"
  "-I${build}/module/fiat" "-I${build}/ecwam/src/ecwam/intfb/interfaces/include")
for mode in strict vector; do
  fp=-fpe0
  [[ ${mode} != vector ]] || fp=-fpe3
  ifort "${flags[@]}" "${fp}" -DWAM_GPU -DSNONLIN=SNONLIN_REFERENCE \
    -c "${root}/src/ecwam/snonlin.F90" -o "${out}/reference-${mode}.o"
  ifort "${flags[@]}" "${fp}" \
    -c "${root}/src/ecwam/snonlin.F90" -o "${out}/fused-${mode}.o"
  ifort "${flags[@]}" -fpe0 "${root}/tests/snonlin_fusion_check.F90" \
    "${out}/reference-${mode}.o" "${out}/fused-${mode}.o" \
    "-L${build}/lib64" -lecwam_sp -lfiat -lparkind_sp \
    "-Wl,-rpath,${build}/lib64" -o "${out}/check-${mode}"
  "${out}/check-${mode}" | tee "${out}/${mode}.log"
done
