#!/usr/bin/env python3
"""Check streaming comparisons across odd subrecords and buffer boundaries."""
import json
from pathlib import Path
import struct
import subprocess
import sys

tool = str(Path(sys.argv[1]).resolve())
output = Path(sys.argv[2]).resolve()
output.mkdir(exist_ok=False)


def record(data, split=None):
    if split is None:
        return struct.pack(">i", len(data)) + data + struct.pack(">i", len(data))
    rest = len(data) - split
    return (struct.pack(">i", -split) + data[:split] + struct.pack(">i", split) +
            struct.pack(">i", rest) + data[split:] + struct.pack(">i", -rest))


def run(left, right, raw=False):
    command = [tool] + (["--raw-little"] if raw else []) + [str(left), str(right)]
    return json.loads(subprocess.check_output(command, universal_newlines=True))


reference = [1.0, 0.0, 2.0]
actual = [1.0 + 2.0 ** -23, -0.0, 2.5]
left = output / "reference.bin"
right = output / "actual.bin"
left.write_bytes(record(struct.pack(">3f", *reference), 5) + record(struct.pack(">f", 3.0)))
right.write_bytes(record(struct.pack(">3f", *actual), 5) + record(struct.pack(">f", 3.0)))
result = run(left, right)
assert result["values"] == 4 and result["records"] == 2 and result["subrecords"] == 3
assert not result["bitwise_equal"] and result["differing_bit_patterns"] == 3
assert result["differing_values"] == 2 and result["tolerance_failures"] == 1
assert result["max_absolute_error"] == 0.5 and result["max_scaled_error"] == 0.25
assert result["worst_flat_index"] == 2
assert run(left, left)["bitwise_equal"]

left = output / "reference-raw.bin"
right = output / "actual-raw.bin"
left.write_bytes(struct.pack("<3f", *reference))
right.write_bytes(struct.pack("<3f", *actual))
result = run(left, right, raw=True)
assert result["values"] == 3 and result["tolerance_failures"] == 1
assert result["differing_bit_patterns"] == 3 and result["differing_values"] == 2

# The second subrecord exceeds the C reader's block size and starts mid-float.
count = 262149
left = output / "reference-block.bin"
right = output / "actual-block.bin"
left.write_bytes(record(struct.pack(">f", 1.0) * count, 7))
right.write_bytes(record(struct.pack(">f", 1.0) * (count - 1) + struct.pack(">f", 1.25), 7))
result = run(left, right)
assert result["values"] == count and result["worst_flat_index"] == count - 1
assert result["tolerance_failures"] == 1 and result["max_absolute_error"] == 0.25

bad = output / "truncated.bin"
bad.write_bytes(struct.pack(">i", 4) + b"xx")
assert subprocess.call([tool, str(bad), str(bad)]) != 0
bad.write_bytes(struct.pack("<f", float("nan")))
assert subprocess.call([tool, "--raw-little", str(bad), str(bad)]) != 0
print("Streaming comparison checks passed.")
