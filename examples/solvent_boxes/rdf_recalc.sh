#!/bin/bash
# Recomputes the aggregation RDFs of an already built box.
# Usage (from the box folder): bash rdf_recalc.sh ":ANT@N1" ":ANT@N1" [first_frame]
A=${1:?mascara A}; B=${2:-$1}; F=${3:-201}
cat > rdf_recalc.in << IN
parm mix.parm7
trajin npt.nc $F last
radial out rdf_cosolv_cosolv.dat 0.1 12.0 volume $A $B
radial out rdf_cosolv_wat.dat 0.1 12.0 volume :WAT@O $A
IN
cpptraj -i rdf_recalc.in > rdf_recalc.log 2>&1 || { tail -5 rdf_recalc.log; exit 1; }
python - << 'PY'
import numpy as np
for f, lab in (('rdf_cosolv_cosolv.dat','cosolvente-cosolvente'),('rdf_cosolv_wat.dat','cosolvente-agua')):
    d=np.loadtxt(f); r,g=d[:,0],d[:,1]
    m=(r>2)&(r<6); k=np.argmax(np.where(m,g,-1))
    print('%-22s max g=%.2f en r=%.1f | media g(9-12)=%.2f'%(lab,g[k],r[k],g[r>9].mean()))
PY
