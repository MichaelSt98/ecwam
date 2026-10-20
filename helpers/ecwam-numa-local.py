#!/usr/bin/env python3
"""Execute one MPI rank with memory bound to its CPU-affinity NUMA node.

Invoke inside the MPI/Slurm launcher, after rank CPU binding has been applied.
Requires Linux sysfs and numactl. CPU affinity must be contained in one NUMA node.
"""
import argparse
import os
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Rank executable and arguments")
    args = parser.parse_args()
    if not args.command:
        parser.error("Pass the rank executable and its arguments")
    nodes = {int(path.name[4:]) for cpu in os.sched_getaffinity(0)
             for path in Path("/sys/devices/system/cpu/cpu%d" % cpu).glob("node[0-9]*")}
    if len(nodes) != 1:
        parser.error("Rank CPU affinity must identify exactly one NUMA node; found %s" % sorted(nodes))
    numactl = shutil.which("numactl")
    if not numactl:
        parser.error("numactl is required")
    command = [numactl, "--membind=" + str(next(iter(nodes)))] + args.command
    os.execv(numactl, command)


if __name__ == "__main__":
    main()
