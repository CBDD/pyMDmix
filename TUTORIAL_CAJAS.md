# Building mixed-solvent boxes with OPC water, step by step

This document describes how a mixed-solvent box (organic cosolvent in OPC water) is built, validated and registered in pyMDMix. It applies to the cosolvents of the library (ETA, ANT, MAM, ISO, ISO5, MOH, PYR1 and ION) and to any new cosolvent whose parameters are already stored in an `.off` file. The files are in `examples/solvent_boxes/`.

> Validation status. The ethanol box (ETAO) was built and validated on a real installation: density 0.9678 ± 0.0052 g/cm³ (approximate experimental value for 17 % w/w at 300 K: 0.97 g/cm³), no drift between the halves of the trajectory and 1672 EPW virtual sites. The other cosolvent profiles follow the same procedure and are validated one by one. Every box must pass the validation of section 5 before it is used.

## 1. Rationale

1. **Cosolvents are not reparametrized.** The RESP charges and atom types of the cosolvent were derived in the gas phase or with GAFF and do not depend on the water model. What changes is the solvent-water balance, which is why the box is rebuilt, re-equilibrated and validated (density and aggregation).
2. **OPC is a four-point model.** The water molecule includes a virtual site, `EPW`. Therefore `leaprc.water.opc` must be loaded, PDB files must be written with `include_ep`, and the number of `EPW` sites must be checked against the number of waters.
3. **ff19SB and OPC were designed together.** The script loads `leaprc.protein.ff19SB` and `leaprc.gaff2` (organic cosolvents use GAFF2 atom types).
4. **Box size.** The original MDMix boxes measure about 20 Å. With a 9 Å cutoff, `pmemd` requires the radius of the largest sphere that fits in the cell to be larger than cutoff + 2 Å, and `pmemd.cuda` rejects small boxes ("Small box detected"). The script replicates the original box `REP` times per side (by default, the smallest integer that gives at least 36 Å). The mixture fraction does not change.

## 2. Requirements

```bash
conda activate pymdmix3
source $AMBERHOME/amber.sh
which tleap cpptraj packmol pmemd pmemd.cuda
ls $AMBERHOME/dat/leap/cmd | grep -E "^leaprc\.(protein.ff19SB|gaff2|water.opc)$"
```

## 3. Choosing the composition with `boxcalc.py`

`boxcalc.py` replaces the `molecule_calculator.ods` spreadsheet. It computes the number of molecules from the volumes of the pure components (V = M / (ρ N_A)). It does not consider the excess volume of the mixture, which is why the final density is obtained from the NPT run.

```bash
cd examples/solvent_boxes
python boxcalc.py table                                   # available cosolvents
python boxcalc.py library                                 # actual concentration of the bundled TIP3P boxes
python boxcalc.py mix ETA --vv 20 --side 20               # molecules for 20 % v/v in 20 Å
python boxcalc.py check ETA --nco 136 --nwat 1672 --side 39.67
python boxcalc.py pure ETA --side 23                      # pure cosolvent box
python boxcalc.py molar UREA --molarity 7 --volume 0.1    # solution of a given molarity
```

Actual concentration of the TIP3P boxes of the library, computed with `boxcalc.py library`:

| Box | Cosolvent + water | % v/v | % w/w | Mole fraction | Molarity |
|---|---|---|---|---|---|
| ETAWAT20 | 17 + 209 | 20.9 | 17.2 | 0.0752 | 3.5 M |
| ANTWAT20 | 17 + 209 | 19.1 | 15.6 | 0.0752 | 3.6 M |
| MAMWAT20 | 18 + 209 | 19.6 | 22.0 | 0.0793 | 3.8 M |
| ISOWAT20 | 14 + 209 | 22.1 | 18.3 | 0.0628 | 2.9 M |
| ISOWAT5 | 4 + 310 | 5.2 | 4.1 | 0.0127 | 0.7 M |
| MOHWAT20 | 56 + 470 | 21.1 | 17.5 | 0.1065 | 5.2 M |
| PYRWAT1 | 1 + 438 | 1.0 | 1.0 | 0.0023 | 0.08 M (see note) |

Remarks on the original spreadsheet:

- In the `ethanol_20%` sheet, the molecular volume of ethanol (`ETHOH_vol`, cell B16) is computed with the molar mass of urea (60.06 g/mol, cell B6) instead of that of ethanol (46.07 g/mol). The correct result is 97.0 Å³ per molecule, not 126.4 Å³. With the correct value, 17 ETA + 209 WAT corresponds to 20.9 % v/v, consistent with the name "20".
- The densities of ANT, MAM, ISO, MOH and PYR used by `boxcalc.py` were not in the spreadsheet. They are literature values (marked `lit` in `boxcalc.py table`); confirm them with your supplier before publishing.
- The PYRWAT1 box stored in the library has 438 waters in a side of about 28 Å (starting density close to 0.60 g/mL), that is, it contains a large void. `build_box.sh` does not use that side: for PYR1 it starts from 23.65 Å, which corresponds to the actual composition.

## 4. Building the box

One empty working folder per solvent:

```bash
mkdir -p ~/boxes/ANT && cd ~/boxes/ANT
nohup bash /path/to/pyMDmix/examples/solvent_boxes/build_box.sh ANT > build.log 2>&1 &
tail -f build.log
```

Bundled profiles: `ETA`, `ANT`, `MAM`, `ISO`, `ISO5`, `MOH`, `PYR1`, `ION`. Each profile takes the TIP3P `.off` from `pyMDMix/data/solventlib/`, the number of molecules and the side of the original box.

Stages:

| Stage | What it does | Main output |
|---|---|---|
| 1 | Extracts the cosolvent units from the TIP3P `.off` | `RES.pdb`, `RES.lib` |
| 2 | Obtains an OPC water molecule and checks the `EPW` atom | `wat_opc.pdb` |
| 3 | Packs with packmol and builds the topology with tleap | `mix.parm7`, `mix.rst7` |
| 4 | Minimizes and heats (50 ps) | `heat.rst7` |
| 5 | 10 ns NPT on GPU | `npt.nc`, `npt.out` |
| 6 | Computes density and box size; selects the representative frame | `6_analysis.log`, `best_frame.txt` |
| 7 | Writes `best.pdb` (with `include_ep`) and the final `.off` | `UNIT.off` |
| 8 | Generates the solvent `.config` and the aggregation RDFs | `SOLV.config`, `rdf_*.dat` |

To resume from a given stage: `bash build_box.sh ANT 5`. Variables that can be overridden from the environment: `REP`, `NPT_NS`, `DISCARD_NS`, `TEMP`, `FFS`, `HEAT_CMD`, `NPT_CMD`, `SRC_OFF`, `SOLV`, `UNIT`, `IONS`.

Examples:

```bash
REP=3 bash build_box.sh MOH            # larger box
NPT_NS=20 bash build_box.sh ISO5       # longer NPT
REP=1 HEAT_CMD=pmemd NPT_CMD=pmemd bash build_box.sh PYR1   # CPU, only if the box is already large
```

Note: the messages printed by the scripts at run time are still in Spanish (for example `LISTO`, `AVISO`, `FALLO`, meaning done, warning and failure).

### 4.1 New cosolvent

If you already have an `.off` with the parametrized cosolvent unit, define the profile through variables:

```bash
python boxcalc.py mix TFE --vv 20 --side 21 --rep 2    # gives NCO and NWAT
COSOLV=TFE SRC_OFF=/path/TFE.off COMP="TFE:NN" NW0=NNN SIDE_EQ=21 \
SOLV=TFEO UNIT=TFEOPC20 CFG_SRC=/path/tfe.config RDFA=":TFE@O1" RDFB=":TFE@O1" \
bash build_box.sh TFE
```

Here `COMP` and `NW0` describe the original, unreplicated box (the script applies `REP`). The `.off` must contain one unit named after each residue in `COMP`. If the cosolvent is not parametrized, charges and parameters must be derived first (antechamber, RESP), which this procedure does not do.

## 5. Validating the box

Minimum conditions before registering:

1. `6_analysis.log`: the density drift between the first and the second half must be close to zero (|drift| smaller than the standard deviation, typically < 0.002 g/cm³), and the mean density must be within a few percent of the experimental value of the mixture at that weight fraction and temperature.
2. Stage 7 must end with "sitios virtuales EPW en best.pdb" (EPW virtual sites in best.pdb) equal to the expected value, and the `.off` must list the cosolvent units, `WAT` and the box unit.
3. Stage 8 prints the height of the first peak and the tail value of g(r) between 9 and 12 Å (it must be close to 1). A cosolvent-cosolvent g(r) with a high first maximum (> 3) below 5 Å indicates aggregation. Compare with the TIP3P box if in doubt. Boxes with few cosolvent molecules (ISO5, PYR1) give noisy RDFs.
4. Inspect the last trajectory visually (no molecules crossing the box in an anomalous way).

To obtain the experimental reference density of an ethanol-water mixture, for example, use the density tables by % w/w and temperature (see the "% w/w" column of `boxcalc.py check`).

## 6. Registering the solvent in pyMDMix

```bash
mkdir ~/mdmix_boxes && cd ~/mdmix_boxes
cp ~/boxes/ANT/ANTOPC20.off ~/boxes/ANT/ANTO.config .
grep -n "name\|objectfile\|boxunit\|watermodel" ANTO.config    # check
mdmix create solvents -f ANTO.config
mdmix info solvents | grep ANTO
```

The `.config` generated by stage 8 reuses the `[PROBES]` and `[TYPES]` sections of the original `.config`. If the command fails with `InvalidPath`, the message states which file was not found (normally `objectfile`).

Solvent names used by this construction: ETAO, ANTO, MAMO, ISOO, ISO5O, MOHO, PYR1O and IONO. They keep the TIP3P boxes (ETA, ANT, ...) separate from the OPC ones. The validated boxes are stored in `examples/solvent_boxes/opc_boxes/` and can be registered all at once with `register_opc_boxes.sh`.

## 7. Testing the box with a system

```ini
[SYSTEM]
NAME = PEP
OFF = pep.off
# With OPC water its force field must be loaded explicitly (or rely on the package default)
EXTRAFF = leaprc.water.opc

[MDSETTINGS]
SOLVENTS = ETAO
NREPL = 1
NANOS = 1
RESTR = FREE
```

```bash
mdmix create project -n PEP_OPC -f pep_opc.cfg
```

In the created replica, `tleap` must load `leaprc.water.opc` and the topology must keep the `EPW` sites.

## 8. Common problems

| Symptom | Cause and solution |
|---|---|
| `max pairlist cutoff must be less than unit cell max sphere radius` or `Small box detected` | The box is too small. Increase `REP`. |
| tleap cannot find an atom type or parameter of the cosolvent | The force field of that type is not loaded (for example `leaprc.gaff2`) or an `frcmod` is missing. Adjust `FFS`. |
| `saveOff` rejects the unit list | Check `7_off.log`. The `saveOff { ... }` syntax was verified with the ethanol box (ETAOPC20.off contains ETA, WAT and the box unit). |
| The number of `EPW` sites does not match | The PDB was written without `include_ep`, or the water residue is not `WAT`. |
| `tleap` warning with "Errors = N" different from 0 | Check `3_build.log` before continuing. |
| The density does not converge | Lengthen the NPT with `NPT_NS=20` and increase `DISCARD_NS`. |
