#!/usr/bin/env python3
"""Check that GSTATS leaves every rank's model log open for final timers."""
import math
from pathlib import Path
import re
import sys

run = Path(sys.argv[1])
ranks = int(sys.argv[2])
assert (run / "model.completed").is_file(), "Model fixture did not complete"
for rank in range(1, ranks + 1):
    path = run / ("logs/model/wam.log.%d" % rank)
    text = path.read_text()
    for label in ("WAVE PROPAGATION TIME", "SOURCE TERM TIME"):
        matches = re.findall(label + r"[^\n]*\n\s*\+\s*([-+0-9.Ee]+)\s*\+", text)
        assert len(matches) == 1, "Missing or duplicate timer in %s: %s" % (path, label)
        value = float(matches[0])
        assert math.isfinite(value) and value >= 0, (path, label, value)
    owners = re.findall(r"ON PE\s*:\s*(\d+)", text)
    assert owners == [str(rank)], (path, owners)
print("Verified retained timer footers on all %d ranks." % ranks)
