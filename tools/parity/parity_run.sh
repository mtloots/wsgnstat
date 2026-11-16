#!/bin/sh
# Byte-identical parity of the R and Python fronts over every entry point of the shared C back-end.
# Usage: sh parity_run.sh [R library path] [python executable]
# The guard runs first: a parity verdict against a stale copy of the C source means nothing.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
sh "$ROOT/parity_sync_guard.sh" "$ROOT/wsgnstat/src" "$ROOT/wsgnstat-py"
RLIB="${1:-}"
PY="${2:-python3}"
TMP=${TMPDIR:-/tmp}
Rscript "$HERE/parity_R.R" $RLIB > "$TMP/wsgn_parity_R.txt"
"$PY" "$HERE/parity_py.py" > "$TMP/wsgn_parity_py.txt"
NR=$(wc -l < "$TMP/wsgn_parity_R.txt" | tr -d ' ')
NP=$(wc -l < "$TMP/wsgn_parity_py.txt" | tr -d ' ')
if cmp -s "$TMP/wsgn_parity_R.txt" "$TMP/wsgn_parity_py.txt"; then
  echo "VERDICT: BYTE-IDENTICAL -- $NR values (%.15e) agree between R and Python"
else
  ND=$(diff "$TMP/wsgn_parity_R.txt" "$TMP/wsgn_parity_py.txt" | grep -c '^<' || true)
  echo "VERDICT: DIFFER -- R $NR lines, Python $NP lines, $ND differing; first differences:"
  diff "$TMP/wsgn_parity_R.txt" "$TMP/wsgn_parity_py.txt" | head -20
  exit 1
fi
