# Rank-local NUMA memory placement

`helpers/ecwam-numa-local.py` runs inside the MPI/Slurm launcher, after rank CPU
binding. It derives the rank's NUMA node from its actual allowed CPU set and
executes `numactl --membind=NODE` before the model starts. It requires Linux
sysfs, numactl and CPU affinity contained in one NUMA node. Each rank derives its
own node; there is no fixed node-0 policy.

Example prefix for the tested allocation, after sourcing the build environment:

```sh
export OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE
export LAUNCH_PARALLEL="srun --nodes=4 --ntasks=64 --ntasks-per-node=16 --cpus-per-task=8 --hint=nomultithread --distribution=block:block --cpu-bind=cores python3 /absolute/ecwam/helpers/ecwam-numa-local.py"
ecwam-run-model --config=CONFIG --run-dir=NEW_RUN --prec=sp
```

The model driver appends the executable. The wrapper also accepts a rank
executable and its arguments directly. This policy affects the entire process,
including initialization; it is an explicit launcher choice.

## Matched performance evidence

Job **30352747**, 2026-10-20, completed 70 O1280 cases in 2h23m45s on four AMD
EPYC 7742 nodes, 64 MPI ranks and eight OpenMP threads/rank. Its fourteen
conditions separate source changes, compiler flags, coefficients and binding.
Each condition has one warm-up and four interleaved timed rounds.

The strict binding contrast uses the **identical frozen baseline library** in
both conditions; only the NUMA policy differs. All ranks record actual affinity
and selected policy. Primary is inclusive complete propagation, measured outside
`CALL PROPAG_WAM`, including allocation and cleanup, and cross-checked against
the retained per-rank built-in timers.

| Median of four runs | Default memory policy | Rank-local binding |
|---|---:|---:|
| Mean accumulated PE complete propagation | 37.158234 s | 32.793633 s |
| Mean accumulated thread PROPAGS2 self time | 12.097423 s | 13.721139 s |
| Model wall time, including initialization/I/O | 115 s | 111.5 s |

Complete propagation reduction by ratio of medians: **11.75%**. Within-round
reductions: **11.91%, 16.51%, 1.45%, 34.68%**, all positive. The mean kernel time
increases, demonstrating why kernel-only selection would miss this placement
benefit. Accumulated PE times are not a synchronized global critical-path trace.
Arrival waiting and placement imbalance contribute to the enclosing result;
these measurements do not attribute every exchange-region second to the network.

The contrast against the original strict PROPAGS2 library is also favorable in
all four rounds (0.74%, 24.64%, 12.24%, 11.12%). Its source differs from the
baseline used for the identical-library contrast, so use the latter to isolate
the placement effect.

## Numerical and longer-run checks

- O48 MPI4 full restart spectrum: all 3,381,264 values are bit-identical with the
  normal-use wrapper enabled/disabled under matched CPU placement.
- O1280 18-hour integration, job **30372070**: all **6,106,530,384** final
  spectrum values and all **136 raw field/date dumps** are byte-identical for
  bound/unbound strict accepted sources. The same binding comparison is bitwise
  in the optional SINPUT-vector mode.
- The longer diagnostic run observes complete propagation 128.207188 s unbound
  versus 126.304844 s bound. Model times are 409 s versus 426 s; those include
  dense raw field output and a full binary checkpoint. This single diagnostic
  comparison is not a repeatable production model speedup estimate.
- A local test with affinity spanning two NUMA nodes was correctly rejected;
  the subsequent aligned MPI test passed. SMT sibling CPU IDs do not by
  themselves imply oversubscription.

These are standalone O1280 / Zen 2 results for the stated decomposition. Memory
placement gains depend on process footprint, affinity and NUMA pressure; coupled
IFS performance is not established by these measurements.

Artifacts:

```text
/ec/res4/scratch/nams/ecwam-benchmarks/propagation-series-verified-20261020/timing-summary.json
/ec/res4/scratch/nams/ecwam-benchmarks/accepted-series-long-20261020/large-validation-summary.json
tmp/performance-series-20261020/numerics-binding-o48-numa-aligned/spectrum-comparison.json
```
