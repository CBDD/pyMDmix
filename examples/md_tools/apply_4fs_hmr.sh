#!/bin/bash
# Optional: switch a pyMDMix project to 4 fs time steps with hydrogen mass repartitioning (HMR).
#
# Usage (with the pymdmix3 environment active, which includes parmed):
#   cd MY_PROJECT/MD
#   bash /path/to/pymdmix/examples/md_tools/apply_4fs_hmr.sh
#
# 1. Edits every */*/*.in file (min, eq and md of each replica):
#      dt 0.002 -> 0.004;  ntwx and ntpr 500 -> 250;  nstlim 500000 -> 250000 (md.in and eq5.in);
#      nstlim 100000 -> 50000 (eq1 to eq4);  iwrap 1 -> 0 in eq1 to eq5 and md.
#    Every ps of simulation yields the same number of frames and each md file still covers 1 ns.
# 2. Repartitions the hydrogen mass (parmed HMassRepartition, 3.024 Da by default) in the topology of
#    the first replica and copies the result to the others, using each replica's own topology name.
#    pyMDMix replicas start from the same solvated system (same .rst7), so one topology is enough.
#    A copy *.parm7.no_hmr of every original topology is kept.
# It is idempotent: running it again changes nothing (marker file HMR_4FS_DONE).
# Variables: HMASS=3.024 (hydrogen mass after repartitioning).

set -e
cd "${1:-.}"
HMASS=${HMASS:-3.024}
command -v python > /dev/null || { echo "ERROR: python is not in the PATH (conda activate pymdmix3)"; exit 1; }
python -c "import parmed" 2> /dev/null || { echo "ERROR: parmed is missing (conda activate pymdmix3)"; exit 1; }
ls */*/*.in > /dev/null 2>&1 || { echo "ERROR: no */*/*.in files here; run the script from the project MD folder"; exit 1; }

echo "== Input files: dt = 4 fs"
for f in */*/*.in; do
  sed -i -e 's/dt=0.002/dt=0.004/' -e 's/ntwx=500/ntwx=250/' -e 's/ntpr=500/ntpr=250/' \
         -e 's/nstlim= 100000/nstlim= 50000/' "$f"
done
for f in */*/md.in */*/eq5.in; do
  [ -f "$f" ] && sed -i 's/nstlim=500000/nstlim=250000/' "$f"
done
for f in */*/eq[1-5].in */*/md.in; do
  [ -f "$f" ] && sed -i 's/iwrap=1/iwrap=0/' "$f"
done
grep -h "dt=" */*/*.in | sort | uniq -c

echo "== Hydrogen mass repartitioning (${HMASS} Da)"
if [ -f HMR_4FS_DONE ]; then
  echo "   already applied (HMR_4FS_DONE); not repeated"
else
  REPS=$(ls -d */ | tr -d '/' | while read d; do ls "$d"/*.parm7 > /dev/null 2>&1 && echo "$d"; done)
  FIRST=$(echo "$REPS" | head -1)
  TOP1=$(ls "$FIRST"/*.parm7 | head -1)
  python - "$TOP1" "$HMASS" << 'PYEOF'
import sys, shutil, parmed
from parmed.tools import HMassRepartition
top, mass = sys.argv[1], float(sys.argv[2])
shutil.copy(top, top + '.no_hmr')
p = parmed.load_file(top)
HMassRepartition(p, str(mass)).execute()
p.save(top, overwrite=True)
hs = [a.mass for a in p.atoms if a.atomic_number == 1 and a.residue.name not in ('WAT', 'HOH')]
print('   %s: minimum H mass of solute/cosolvent = %.3f, maximum = %.3f' % (top, min(hs), max(hs)))
PYEOF
  for r in $REPS; do
    [ "$r" = "$FIRST" ] && continue
    t=$(ls "$r"/*.parm7 | head -1)
    cp "$t" "$t.no_hmr"
    cp "$TOP1" "$t"
    echo "   $t <- $TOP1"
  done
  touch HMR_4FS_DONE
fi
echo "DONE. Check an input file (for example: head -20 */md/md.in) and launch with run_replicas.sh."
