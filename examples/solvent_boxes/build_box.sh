#!/bin/bash
# =============================================================================
# Building COSOLVENT / OPC WATER boxes for pyMDMix
#
# Usage (from an EMPTY working folder per solvent, with pymdmix3 and AMBER active):
#   mkdir ANT && cd ANT
#   bash /path/to/build_box.sh ANT        # runs all stages
#   bash /path/to/build_box.sh ANT 5      # resumes from stage 5
#
# Cosolvents with a built-in profile: ETA ANT MAM ISO ISO5 MOH PYR1 ION
# (composition and side taken from the TIP3P boxes in pyMDMix/data/solventlib).
# For another cosolvent (with parameters already available in an .off) define
# COSOLV, SRC_OFF, COMP, NW0, SIDE_EQ, SOLV and UNIT as environment variables;
# boxcalc.py computes COMP and NW0 for a given % v/v.
#
# Stages:
#   1  Extract the cosolvents (coordinates and parameters) from the TIP3P .off
#   2  Get one OPC water molecule (with the EPW virtual site)
#   3  Pack with packmol and generate parm7/rst7 with tleap
#   4  Minimize and heat
#   5  Equilibrate in NPT (GPU)
#   6  Analyze density and box, choose the representative frame
#   7  Build the final .off with tleap
#   8  Write the solvent .config and the aggregation RDF
# =============================================================================
set -euo pipefail

COSOLV=${1:-}
START=${2:-1}
[ -n "$COSOLV" ] || { echo "Uso: bash build_box.sh COSOLV [etapa_inicial]   (ETA ANT MAM ISO ISO5 MOH PYR1 ION)"; exit 1; }

# ------------------------------ profiles --------------------------------------
# SRC: source TIP3P .off | CFG: original .config | COMP: cosolvent residues:number (original box)
# NW0: waters in the original box | SIDE_EQ: equilibrated side (A) of the original box
# PARMS: PARMSET objects (cosolvent-specific parameters) stored in the .off | SOLV: new solvent name | UNIT: new box unit | RDFA/RDFB: masks for the cosolvent-cosolvent RDF
case "$COSOLV" in
  ETA)  P_SRC=ETAWAT20.off;    P_CFG=eta.config;  P_COMP="ETA:17";           P_NW0=209; P_SIDE=19.99; P_SOLV=ETAO;  P_UNIT=ETAOPC20; P_RDFA=":ETA@O1"; P_RDFB=":ETA@O1";;
  ANT)  P_PARMS=ANTparms; P_SRC=ANTWAT20.off;    P_CFG=ant.config;  P_COMP="ANT:17";           P_NW0=209; P_SIDE=19.86; P_SOLV=ANTO;  P_UNIT=ANTOPC20; P_RDFA=":ANT@N1"; P_RDFB=":ANT@N1";;
  MAM)  P_SRC=MAMWAT20.off;    P_CFG=mam.config;  P_COMP="MAM:18";           P_NW0=209; P_SIDE=19.85; P_SOLV=MAMO;  P_UNIT=MAMOPC20; P_RDFA=":MAM@O1"; P_RDFB=":MAM@O1";;
  ISO)  P_SRC=ISOWAT20.off;    P_CFG=iso.config;  P_COMP="ISO:14";           P_NW0=209; P_SIDE=20.06; P_SOLV=ISOO;  P_UNIT=ISOOPC20; P_RDFA=":ISO@O1"; P_RDFB=":ISO@O1";;
  ISO5) P_SRC=ISOWAT5_ok.off;  P_CFG=iso5.config; P_COMP="ISO:4";            P_NW0=310; P_SIDE=21.41; P_SOLV=ISO5O; P_UNIT=ISOOPC5;  P_RDFA=":ISO@O1"; P_RDFB=":ISO@O1";;
  MOH)  P_SRC=MOHWAT20.off;    P_CFG=moh.config;  P_COMP="MOH:56";           P_NW0=470; P_SIDE=26.23; P_SOLV=MOHO;  P_UNIT=MOHOPC20; P_RDFA=":MOH@O1"; P_RDFB=":MOH@O1";;
  PYR1) P_PARMS=PYRparm;  P_SRC=PYRWAT1.off;     P_CFG=pyr1.config; P_COMP="PYR:1";            P_NW0=438; P_SIDE=23.65; P_SOLV=PYR1O; P_UNIT=PYROPC1;  P_RDFA=":PYR@N1"; P_RDFB=":PYR@N1";;
  ION)  P_SRC=IONWAT20.off;    P_CFG=ion.config;  P_COMP="CN3:12 COO:12";   P_NW0=315; P_SIDE=22.17; P_SOLV=IONO;  P_UNIT=IONOPC20; P_RDFA=":CN3@N1"; P_RDFB=":COO@O1";;
  *)    P_SRC=; P_CFG=; P_COMP=; P_NW0=; P_SIDE=; P_SOLV=; P_UNIT=; P_RDFA=; P_RDFB=;;
esac
P_PARMS=${P_PARMS:-}

# ---------------------------- parameters -------------------------------------
# Any variable can be overridden from the environment (for example REP=3 NPT_NS=20).
SRC_DIR=${SRC_DIR:-$(python -c "import os, pyMDMix; print(os.path.join(os.path.dirname(pyMDMix.__file__), 'data', 'solventlib'))" 2>/dev/null || echo "")}
SRC_OFF=${SRC_OFF:-$SRC_DIR/$P_SRC}
CFG_SRC=${CFG_SRC:-$SRC_DIR/$P_CFG}
COMP=${COMP:-$P_COMP}              # "RES:n RES:n" for the ORIGINAL (unreplicated) box
NW0=${NW0:-$P_NW0}
SIDE_EQ=${SIDE_EQ:-$P_SIDE}        # side of the original equilibrated box (A)
SOLV=${SOLV:-$P_SOLV}              # solvent name in pyMDMix
UNIT=${UNIT:-$P_UNIT}              # box unit name in the final .off
PARMS=${PARMS:-$P_PARMS}          # PARMSET objects in the .off (e.g. ANTparms); empty if none
RDFA=${RDFA:-$P_RDFA}
RDFB=${RDFB:-$P_RDFB}
for v in SRC_OFF COMP NW0 SIDE_EQ SOLV UNIT; do
  [ -n "${!v}" ] || { echo "ERROR: perfil '$COSOLV' desconocido o falta $v. Define $v en el entorno."; exit 1; }
done
# A minimum side of ~36 A for the 9 A cutoff (pmemd requires sphere radius > cutoff + 2 A and
# pmemd.cuda rejects small boxes). The original box is replicated REP times per side; the
# composition does not change.
REP=${REP:-$(awk "BEGIN{r=int(36.0/$SIDE_EQ); if (r*$SIDE_EQ<36.0) r++; if (r<1) r=1; print r}")}
R3=$(( REP * REP * REP ))
NWAT=$(( NW0 * R3 ))
SIDE0=${SIDE0:-$(awk "BEGIN{printf \"%.2f\", $SIDE_EQ*1.05*$REP}")}   # slightly looser initial box
TEMP=${TEMP:-300.0}               # K
NPT_NS=${NPT_NS:-10}              # duration of the NPT phase (ns)
DISCARD_NS=${DISCARD_NS:-2}       # ns of NPT discarded in the analysis
# Solute force fields: ff19SB (designed to be used with OPC) and gaff2 (organic cosolvents)
FFS=${FFS:-"leaprc.protein.ff19SB leaprc.gaff2"}
PMEMD_CPU=${PMEMD_CPU:-pmemd}
# Heating and NPT on GPU. If the box is small (REP=1) use CPU: NPT_CMD=pmemd HEAT_CMD=pmemd
PMEMD_GPU=${PMEMD_GPU:-pmemd.cuda}
HEAT_CMD=${HEAT_CMD:-$PMEMD_GPU}
NPT_CMD=${NPT_CMD:-$PMEMD_GPU}

NSTEPS=$(( NPT_NS * 500000 ))     # dt = 2 fs
NFRAMES=$(( NPT_NS * 100 ))       # one frame every 10 ps
DISCARD_FR=$(( DISCARD_NS * 100 ))
RESLIST=$(for c in $COMP; do echo -n "${c%%:*} "; done)
SRC_LINES=$(for f in $FFS; do echo "source $f"; done)
OBJLIST="$RESLIST $PARMS"
LOAD_LIBS=$(for r in $OBJLIST; do echo "loadOff $r.lib"; done)

# ---------------------------- checks -----------------------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: no se encuentra '$1' en el PATH"; exit 1; }; }
need tleap; need cpptraj; need packmol; need python
command -v "$PMEMD_CPU" >/dev/null 2>&1 || PMEMD_CPU=sander
need "$PMEMD_CPU"; need "${HEAT_CMD%% *}"; need "${NPT_CMD%% *}"
[ -f "$SRC_OFF" ] || { echo "ERROR: no existe SRC_OFF=$SRC_OFF"; exit 1; }
python -c "import netCDF4, numpy" || { echo "ERROR: faltan netCDF4/numpy en este entorno"; exit 1; }
[ -n "${AMBERHOME:-}" ] || { echo "ERROR: AMBERHOME no esta definido (source \$AMBERHOME/amber.sh)"; exit 1; }
for f in $FFS leaprc.water.opc; do
  [ -f "$AMBERHOME/dat/leap/cmd/$f" ] || { echo "ERROR: no existe $AMBERHOME/dat/leap/cmd/$f (ajusta FFS=... o revisa tu instalacion)"; echo "Disponibles:"; ls "$AMBERHOME/dat/leap/cmd" | grep -E "^leaprc\.(protein|ff|gaff|water)" ; exit 1; }
done
COMPTXT=$(for c in $COMP; do echo -n "$(( ${c##*:} * R3 )) ${c%%:*} + "; done)
echo ">> $COSOLV | Etapa inicial: $START | REP=$REP | ${COMPTXT}$NWAT WAT | caja inicial $SIDE0 A | NPT ${NPT_NS} ns | $SOLV / $UNIT"

run_stage() { [ "$START" -le "$1" ]; }

# ------------------------------ STAGE 1 --------------------------------------
if run_stage 1; then
echo "== Etapa 1: extraer cosolvente(s) [$RESLIST] y parametros [${PARMS:-ninguno}] de $SRC_OFF"
{
  echo "$SRC_LINES"
  echo "loadOff $SRC_OFF"
  for r in $RESLIST; do
    echo "savePdb $r $r.pdb"
  done
  for r in $OBJLIST; do
    echo "saveOff $r $r.lib"
  done
  echo "quit"
} > 1_cosolv.leap
tleap -f 1_cosolv.leap > 1_cosolv.log || { echo "ERROR: tleap fallo con 1_cosolv.leap, revisa 1_cosolv.log"; tail -5 1_cosolv.log; exit 1; }
for r in $OBJLIST; do
  [ -s $r.lib ] && { [ -s $r.pdb ] || [[ " $PARMS " == *" $r "* ]]; } || { echo "ERROR etapa 1: falta $r.pdb o $r.lib (el residuo $r no esta en $SRC_OFF?). Revisa 1_cosolv.log"; exit 1; }
done
fi

# ------------------------------ STAGE 2 --------------------------------------
if run_stage 2; then
echo "== Etapa 2: molecula de agua OPC"
cat > 2_wat.leap << EOF
$SRC_LINES
source leaprc.water.opc
savePdb WAT wat_opc.pdb
quit
EOF
tleap -f 2_wat.leap > 2_wat.log || { echo "ERROR: tleap fallo con 2_wat.leap, revisa 2_wat.log"; tail -5 2_wat.log; exit 1; }
[ -s wat_opc.pdb ] || { echo "ERROR etapa 2: revisa 2_wat.log"; exit 1; }
WRES=$(awk '/^ATOM|^HETATM/{print substr($0,18,3); exit}' wat_opc.pdb | tr -d ' ')
WATOMS=$(awk '/^ATOM|^HETATM/{printf "%s ", substr($0,13,4)}' wat_opc.pdb)
echo "   residuo del agua: '$WRES' | atomos: $WATOMS"
[ "$WRES" = "WAT" ] || { echo "ERROR: el residuo de agua OPC se llama '$WRES', no 'WAT'. Avisame antes de seguir."; exit 1; }
echo "$WATOMS" | grep -q "EPW" || { echo "ERROR: no aparece el sitio virtual EPW en wat_opc.pdb"; exit 1; }
fi

# ------------------------------ STAGE 3 --------------------------------------
if run_stage 3; then
echo "== Etapa 3: empaquetar y construir topologia"
# 1 A margin on each face: packmol does not know about periodicity and, without a margin, two molecules
# placed on opposite faces end up almost overlapping once the periodic box is applied.
PLO=1.0
PHI=$(awk "BEGIN{printf \"%.2f\", $SIDE0-1.0}")
PACK_BLOCKS=$(for c in $COMP; do
  printf 'structure %s.pdb\n  number %d\n  inside box %s %s %s %s %s %s\nend structure\n\n' "${c%%:*}" $(( ${c##*:} * R3 )) $PLO $PLO $PLO $PHI $PHI $PHI
done)
cat > mix.inp << EOF
tolerance 2.0
filetype pdb
output mix_packed.pdb
add_amber_ter
seed 12345

$PACK_BLOCKS
structure wat_opc.pdb
  number $NWAT
  inside box $PLO $PLO $PLO $PHI $PHI $PHI
end structure
EOF
packmol < mix.inp > 3_packmol.log
grep -q "Success" 3_packmol.log || { echo "ERROR etapa 3: packmol no convergio (ver 3_packmol.log)"; exit 1; }
python - mix_packed.pdb $SIDE0 << 'PYEOF'
import sys
import numpy as np
from scipy.spatial import cKDTree
pdb, side = sys.argv[1], float(sys.argv[2])
xyz, res = [], []
for l in open(pdb):
    if l.startswith(('ATOM', 'HETATM')):
        xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        res.append(l[22:27])            # residue number + insertion code
xyz = np.array(xyz) % side
res = np.array([hash(r) for r in res])
pairs = cKDTree(xyz, boxsize=side).query_pairs(1.5, output_type='ndarray')
inter = pairs[res[pairs[:, 0]] != res[pairs[:, 1]]] if len(pairs) else pairs
if len(inter):
    dd = xyz[inter[:, 0]] - xyz[inter[:, 1]]
    dd -= side * np.round(dd / side)
    d = np.linalg.norm(dd, axis=1)
    print('   AVISO: %d pares de atomos de moleculas distintas a menos de 1.5 A (minimo %.2f A)' % (len(inter), d.min()))
    if d.min() < 0.8:
        sys.exit('ERROR etapa 3: hay solapamientos graves en mix_packed.pdb')
else:
    print('   empaquetado sin contactos intermoleculares < 1.5 A (con periodicidad)')
PYEOF

cat > 3_build.leap << EOF
$SRC_LINES
source leaprc.water.opc
$LOAD_LIBS
mix = loadPdb mix_packed.pdb
set mix box { $SIDE0 $SIDE0 $SIDE0 }
check mix
charge mix
saveAmberParm mix mix.parm7 mix.rst7
quit
EOF
tleap -f 3_build.leap > 3_build.log || { echo "ERROR: tleap fallo con 3_build.leap, revisa 3_build.log"; tail -5 3_build.log; exit 1; }
[ -s mix.parm7 ] && [ -s mix.rst7 ] || { echo "ERROR etapa 3: revisa 3_build.log"; exit 1; }
echo "   carga total reportada por tleap:"; grep -i "total unperturbed charge\|Total perturbed charge" 3_build.log || true
grep -E "Errors = [1-9]|Fatal Error|FATAL" 3_build.log && { echo "AVISO: tleap reporto errores, revisa 3_build.log"; } || true
# Old torsion parameters (PARMSET without SCEE/SCNB) leave SCEE/SCNB uninitialized (garbage or 0) in the
# parm7 and give 1-4 VDW = -Infinity. Out-of-range values are detected and set to the Amber values (1.2 / 2.0).
python - << 'PYEOF'
import parmed
p = parmed.load_file('mix.parm7', 'mix.rst7')
n = 0
for t in p.dihedral_types:
    if not (1.0 <= t.scnb <= 3.0 and 1.0 <= t.scee <= 2.0):
        t.scee, t.scnb = 1.2, 2.0
        n += 1
if n:
    p.save('mix.parm7', overwrite=True)
    print('   AVISO: %d tipos de torsion con SCEE/SCNB fuera de rango; fijados a 1.2/2.0 en mix.parm7' % n)
PYEOF
fi

# ------------------------------ STAGE 4 --------------------------------------
if run_stage 4; then
echo "== Etapa 4: minimizacion y calentamiento"
cat > min.in << EOF
minimizacion
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=9.0, ntpr=500,
 /
EOF
cat > heat.in << EOF
calentamiento NVT
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=25000, dt=0.002,
  ntc=2, ntf=2, cut=9.0,
  ntb=1,
  ntt=3, gamma_ln=2.0, ig=-1,
  tempi=100.0, temp0=$TEMP,
  ntpr=1000, ntwx=0, ntwr=25000,
 /
EOF
$PMEMD_CPU -O -i min.in  -o min.out  -p mix.parm7 -c mix.rst7  -r min.rst7
$HEAT_CMD -O -i heat.in -o heat.out -p mix.parm7 -c min.rst7  -r heat.rst7
fi

# ------------------------------ STAGE 5 --------------------------------------
if run_stage 5; then
echo "== Etapa 5: NPT de ${NPT_NS} ns con: $NPT_CMD"
cat > npt.in << EOF
NPT
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=$NSTEPS, dt=0.002,
  ntc=2, ntf=2, cut=9.0,
  ntb=2, ntp=1, pres0=1.0, barostat=2,
  ntt=3, gamma_ln=2.0, temp0=$TEMP, ig=-1,
  iwrap=1,
  ntpr=5000, ntwx=5000, ntwr=500000, ioutfm=1,
 /
EOF
$NPT_CMD -O -i npt.in -o npt.out -p mix.parm7 -c heat.rst7 -r npt.rst7 -x npt.nc
fi

# ------------------------------ STAGE 6 --------------------------------------
if run_stage 6; then
echo "== Etapa 6: analisis de densidad y caja"
cat > analyze_density.py << 'EOF'
import sys
import numpy as np
from netCDF4 import Dataset

parm, nc, discard = sys.argv[1], sys.argv[2], int(sys.argv[3])

# total mass (amu) from the MASS section of the parm7 (EPW has zero mass)
lines = open(parm).read().split('\n')
i = next(k for k, l in enumerate(lines) if l.startswith('%FLAG MASS'))
masses = []
for l in lines[i + 2:]:
    if l.startswith('%FLAG'):
        break
    masses += [float(x) for x in l.split()]
mass = sum(masses)

ds = Dataset(nc)
cell = np.array(ds.variables['cell_lengths'][:], dtype=float)
angles = np.array(ds.variables['cell_angles'][:], dtype=float)
nf = len(cell)
if nf <= discard + 10:
    sys.exit('ERROR: muy pocos cuadros (%d) para descartar %d' % (nf, discard))
if np.abs(angles - 90.0).max() > 1e-3:
    sys.exit('ERROR: la caja no es ortogonal, este script asume caja ortogonal')

vol = cell.prod(axis=1)
rho = mass * 1.66053907 / vol          # g/cm3
sl = slice(discard, nf)
v, r, c = vol[sl], rho[sl], cell[sl]
half = len(r) // 2
print('cuadros totales: %d | usados: %d' % (nf, len(r)))
print('masa total (uma): %.3f' % mass)
print('densidad media: %.4f +/- %.4f g/cm3' % (r.mean(), r.std()))
print('  primera mitad: %.4f | segunda mitad: %.4f (deriva: %.4f)' %
      (r[:half].mean(), r[half:].mean(), r[half:].mean() - r[:half].mean()))
print('lado medio de la caja (A): %.3f %.3f %.3f' % tuple(c.mean(axis=0)))
print('volumen medio (A^3): %.1f' % v.mean())
best = int(np.argmin(np.abs(v - v.mean())))
idx = discard + best
print('cuadro representativo (1-based): %d | caja: %.4f %.4f %.4f | densidad %.4f' %
      (idx + 1, cell[idx][0], cell[idx][1], cell[idx][2], rho[idx]))
open('best_frame.txt', 'w').write('%d %.4f %.4f %.4f\n' % (idx + 1, cell[idx][0], cell[idx][1], cell[idx][2]))
np.savetxt('density_vs_frame.dat', np.c_[np.arange(1, nf + 1), rho, cell], fmt='%d %.5f %.4f %.4f %.4f',
           header='frame density(g/cm3) a b c')
EOF
python analyze_density.py mix.parm7 npt.nc "$DISCARD_FR" | tee 6_analysis.log
fi

# ------------------------------ STAGE 7 --------------------------------------
if run_stage 7; then
echo "== Etapa 7: construir el .off final ($UNIT)"
read FRAME BA BB BC < best_frame.txt
cpptraj -p mix.parm7 > 7_cpptraj.log << EOF
trajin npt.nc $FRAME $FRAME
trajout best.pdb pdb include_ep
EOF
[ -s best.pdb ] || { echo "ERROR etapa 7: revisa 7_cpptraj.log"; exit 1; }
NEP=$(grep -c " EPW" best.pdb || true)
echo "   sitios virtuales EPW en best.pdb: $NEP (esperados: $NWAT)"
[ "$NEP" -eq "$NWAT" ] || { echo "ERROR: el numero de EPW no coincide"; exit 1; }

cat > 7_off.leap << EOF
$SRC_LINES
source leaprc.water.opc
$LOAD_LIBS
$UNIT = loadPdb best.pdb
set $UNIT box { $BA $BB $BC }
check $UNIT
saveOff { $OBJLIST WAT $UNIT } $UNIT.off
quit
EOF
tleap -f 7_off.leap > 7_off.log || { echo "ERROR: tleap fallo con 7_off.leap, revisa 7_off.log"; tail -5 7_off.log; exit 1; }
[ -s $UNIT.off ] || { echo "ERROR etapa 7: tleap no escribio $UNIT.off (revisa 7_off.log; si rechazo la lista de unidades, avisame)"; exit 1; }
echo "   unidades guardadas en $UNIT.off:"
sed -n '2,8p' $UNIT.off | grep -v '^!'
echo
echo "LISTO: $UNIT.off. Caja: $BA x $BB x $BC A."
fi

# ------------------------------ STAGE 8 --------------------------------------
if run_stage 8; then
echo "== Etapa 8: archivo .config del solvente y RDF de agregacion"
if [ -f "$CFG_SRC" ]; then
python - "$CFG_SRC" "$SOLV" "$UNIT" "$UNIT.off" > "$SOLV.config" << 'PYEOF'
import re, sys
src, solv, unit, off = sys.argv[1:5]
lines = open(src).read().split('\n')
out, seen_wm, last_box = [], False, None
for l in lines:
    if re.match(r'\s*name\s*=', l):         l = 'name = %s' % solv
    elif re.match(r'\s*info\s*=', l):       l = re.sub(r'\s+$', '', l) + ' in OPC water'
    elif re.match(r'\s*objectfile\s*=', l): l = 'objectfile = %s' % off
    elif re.match(r'\s*boxunit\s*=', l):    l = 'boxunit = %s' % unit; last_box = len(out)
    elif re.match(r'\s*watermodel\s*=', l): l = 'watermodel = OPC'; seen_wm = True
    out.append(l)
if not seen_wm and last_box is not None:
    out.insert(last_box + 1, 'watermodel = OPC')
print('\n'.join(out))
PYEOF
echo "   escrito $SOLV.config (revisa name, objectfile, boxunit, watermodel y PROBES)"
else
echo "   AVISO: no se encontro $CFG_SRC; escribe el .config a mano (ver TUTORIAL_CAJAS.md)"
fi
if [ -f npt.nc ]; then
cat > 8_rdf.in << EOF
parm mix.parm7
trajin npt.nc $(( DISCARD_FR + 1 )) last
radial out rdf_cosolv_cosolv.dat 0.1 12.0 volume $RDFA $RDFB
radial out rdf_cosolv_wat.dat 0.1 12.0 volume :WAT@O $RDFA
EOF
cpptraj -i 8_rdf.in > 8_rdf.log 2>&1 || { echo "ERROR etapa 8: cpptraj fallo, revisa 8_rdf.log"; tail -5 8_rdf.log; exit 1; }
python - << 'PYEOF'
import numpy as np
for f, lab in (('rdf_cosolv_cosolv.dat', 'cosolvente-cosolvente'), ('rdf_cosolv_wat.dat', 'cosolvente-agua (O)')):
    d = np.loadtxt(f)
    r, g = d[:, 0], d[:, 1]
    m = (r > 2.0) & (r < 6.0)
    k = np.argmax(np.where(m, g, -1))
    tail = g[r > 9.0].mean()
    print('   RDF %-22s primer maximo: g=%.2f en r=%.1f A | media g(9-12 A)=%.2f (esperado ~1)' % (lab, g[k], r[k], tail))
PYEOF
echo "   Referencia: un g(r) cosolvente-cosolvente con maximo alto (> 3) en r < 5 A sugiere agregacion."
fi
echo
echo "SIGUIENTE: mdmix create solvents -f $SOLV.config   (con $UNIT.off en la misma carpeta)"
fi
