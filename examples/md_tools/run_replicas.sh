#!/bin/bash
# Automates the run of all replicas of a pyMDMix project with a single equilibration.
#
# Usage (with the pymdmix3 environment and Amber loaded):
#   cd MY_PROJECT/MD
#   bash /path/to/pymdmix/examples/md_tools/run_replicas.sh
#
# 1. Runs the complete COMMANDS.sh (minimization, eq1 to eq5 and md) in the first replica.
# 2. Copies its eq/eq5.rst7 to eq/ of the other replicas (same starting system).
# 3. In the other replicas runs only the md part of COMMANDS.sh (no min or eq).
#    Replicas 2..N run sequentially, or in parallel distributed over several GPUs.
#
# Variables:
#   GPUS="0 1"    available GPUs; with several, replicas 2..N run in parallel in batches (default "0")
#   FIRST=ETAO_1  folder of the replica that equilibrates (default: the first in alphabetical order)
#   ONLY_FIRST=1  run only step 1 and stop
#   RUN_FIRST=0   skip step 1 (the first replica already ran)

set -e
cd "${1:-.}"
GPUS=${GPUS:-0}
REPS=$(ls -d */ | tr -d '/' | while read d; do [ -f "$d/COMMANDS.sh" ] && echo "$d"; done)
[ -n "$REPS" ] || { echo "ERROR: no folders with COMMANDS.sh; run from the project MD folder"; exit 1; }
FIRST=${FIRST:-$(echo "$REPS" | head -1)}
OTHERS=$(echo "$REPS" | grep -vx "$FIRST" || true)
echo "Equilibrating replica: $FIRST | others: $(echo $OTHERS) | GPUs: $GPUS"

if [ "${RUN_FIRST:-1}" = "1" ]; then
  echo "== $FIRST: minimization, equilibration and production"
  (cd "$FIRST" && CUDA_VISIBLE_DEVICES=$(echo $GPUS | awk '{print $1}') sh COMMANDS.sh) || echo "WARNING: COMMANDS.sh of $FIRST exited with a non-zero code"
fi
[ "${ONLY_FIRST:-0}" = "1" ] && exit 0
[ -s "$FIRST/eq/eq5.rst7" ] || { echo "ERROR: $FIRST/eq/eq5.rst7 is missing (the equilibration did not finish)"; exit 1; }

# md part of COMMANDS.sh: from the 'cd ../md' line onwards, run from the replica folder
for r in $OTHERS; do
  cp "$FIRST/eq/eq5.rst7" "$r/eq/"
  sed -n '/^cd \.\.\/md$/,$p' "$r/COMMANDS.sh" | sed '1s#.*#cd md#' > "$r/COMMANDS_md.sh"
done

NG=$(echo $GPUS | wc -w)
i=0
for r in $OTHERS; do
  g=$(echo $GPUS | awk -v k=$(( i % NG + 1 )) '{print $k}')
  echo "== $r: production on GPU $g"
  ( cd "$r" && CUDA_VISIBLE_DEVICES=$g sh COMMANDS_md.sh || echo "WARNING: $r exited with a non-zero code" ) &
  i=$(( i + 1 ))
  # full batch: wait for it to finish before launching more
  [ $(( i % NG )) -eq 0 ] && wait
done
wait
echo "DONE: all replicas finished."
