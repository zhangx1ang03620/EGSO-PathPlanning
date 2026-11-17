#!/bin/bash
# Launch the full equal-NFE experiment: 5 maps x {15000, 5000} MaxFEs x 20 runs,
# split into 5-run chunks, longest jobs first, NPROC single-threaded MATLAB
# processes in parallel. Results -> runs/equal_nfe/, per-job logs -> runs/equal_nfe/logs/.
# Afterwards run Fixup_EarlyFinish and then Aggregate_EqualNFE('<repo>/runs/equal_nfe').
#
# Usage:  NPROC=8 MATLAB=/path/to/matlab ./run_all_jobs.sh
set -u
cd "$(dirname "$0")"
OUT=../../runs/equal_nfe
mkdir -p "$OUT/logs"
NPROC=${NPROC:-8}
MATLAB=${MATLAB:-matlab}
jobs=()
for fe in 15000 5000; do
  for m in 5 4 3 1 2; do
    for a in 1 6 11 16; do
      b=$((a+4))
      jobs+=("$m $fe $a $b")
    done
  done
done
export OUT MATLAB
printf '%s\n' "${jobs[@]}" | xargs -P "$NPROC" -L 1 bash -c '
  m=$0; fe=$1; a=$2; b=$3
  tag=$(printf "map%d_fe%d_r%02d-%02d" $m $fe $a $b)
  if [ -f "$OUT/EqNFE_${tag}.mat" ]; then echo "SKIP $tag"; exit 0; fi
  "$MATLAB" -singleCompThread -batch "Run_EqualNFE_Job($m, $fe, $a, $b)" > "$OUT/logs/${tag}.log" 2>&1
  if [ -f "$OUT/EqNFE_${tag}.mat" ]; then echo "DONE $tag"; else echo "FAIL $tag"; fi
'
echo "ALL JOBS FINISHED"
