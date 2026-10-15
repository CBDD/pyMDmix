# Tutorial rápido de pyMDMix 3 con cajas de solvente OPC

Este tutorial recorre el flujo completo de la versión Python 3 de pyMDMix: preparar el entorno, construir una caja de solvente con agua OPC, crear un proyecto, ajustar el paso de integración, analizar las réplicas y obtener hot spots. Los comandos se escriben para una terminal con AmberTools disponible.

> Estado de validación. El instalador, la interfaz `mdmix`, la lectura de la base de solventes y el análisis de trayectorias se verificaron. Las rutas que dependen de AMBER (tleap, pmemd, cpptraj) y de agua OPC con sitio virtual se están validando en una instalación real; los puntos marcados con **Verificar** indican qué comprobar.

## 1. Entorno

```bash
git clone https://github.com/CBDD/pyMDmix.git
cd pyMDmix
conda env create -f environment.yml        # crea el entorno pymdmix3
conda activate pymdmix3
python -m pip install -e .
source $AMBERHOME/amber.sh                 # tleap, cpptraj, pmemd en el PATH
```

El archivo `environment.yml` fija `numpy<2` y `biopython<1.86` para que no choquen con las herramientas de Python de AmberTools. Además instala `packmol`, que solo se necesita para construir cajas de solvente.

Comprobación:

```bash
which mdmix
python -c "import pyMDMix; print(pyMDMix.__file__)"
mdmix -h
mdmix info solvents
```

Si ya tienes AmberTools instalado fuera de conda (por ejemplo `/opt/ambertools26`), elimina la línea `- ambertools` de `environment.yml` antes de crear el entorno y carga tu instalación con `source /ruta/amber.sh`. En ese caso comprueba que sus herramientas de Python sean compatibles con `numpy<2`.

pmemd y pmemd.cuda no vienen con el paquete `ambertools` de conda. Si no los tienes, `sander` sirve para pruebas pequeñas: cambia `AMBER_PROD_EXE` y `AMBER_MIN_EXE` en `settings.cfg`.

Registro de las cajas OPC incluidas (etanol, acetonitrilo, metanol, isopropanol al 20 % y 5 %, acetamida y la mezcla iónica acetato/metilamonio):

```bash
bash examples/solvent_boxes/register_opc_boxes.sh
mdmix info solvents
```

Prueba mínima de un proyecto (si tienes tleap y cpptraj): `mdmix create template`, ajusta `1YER.cfg` o el PDB que uses con `SOLVENTS = ETAO` y ejecuta `mdmix create project -n PRUEBA -f proyecto.cfg`. En el `.parm7` resultante deben aparecer tantos `EPW` como moléculas `WAT`.

Si `which mdmix` apunta a un script antiguo en `~/.local/bin`, activa primero el entorno `pymdmix3` o llama al módulo directamente: `python -m pyMDMix -h`.

Nombres de campos de fuerza. En AmberTools recientes ya no existe `leaprc.ff14SB`. Esta versión carga por defecto `leaprc.protein.ff19SB`, `leaprc.gaff2` y `leaprc.water.opc` (opción `FF` de `md-settings.cfg`). El modelo de agua por defecto es OPC (`DEF_AMBER_WATBOX`), el buffer de solvatación es de 10 A y la temperatura por defecto es 310 K. Las cajas TIP3P antiguas (`ETAWAT20` y las demás) siguen siendo válidas si su archivo `.config` declara `watermodel = TIP3P`; en ese caso añade `leaprc.water.tip3p` en `EXTRAFF`. Comprueba que existan en tu instalación:

```bash
ls $AMBERHOME/dat/leap/cmd | grep -E "^leaprc\.(protein.ff19SB|gaff|water.opc)"
```

## Pasos de 4 fs y corrida automática de réplicas (opcional)

Después de `mdmix create project`, dentro de la carpeta `MD` del proyecto:

```bash
bash /ruta/a/pymdmix/examples/md_tools/apply_4fs_hmr.sh    # dt = 4 fs, iwrap = 0 y masa de H = 3.024 Da en todas las réplicas
bash /ruta/a/pymdmix/examples/md_tools/run_replicas.sh     # una equilibración y producción de todas las réplicas
GPUS="0 1" bash /ruta/a/pymdmix/examples/md_tools/run_replicas.sh   # réplicas 2..N en paralelo en dos GPU
```

`apply_4fs_hmr.sh` edita los `*.in` (dt 0.004, `ntwx` y `ntpr` 250, `nstlim` 250000 en `md` y `eq5`, 50000 en `eq1` a `eq4`) y reparte la masa del hidrógeno en la topología de la primera réplica con parmed. Esa topología se copia a las demás, que parten del mismo sistema solvatado. Guarda copias `*.parm7.no_hmr` y es idempotente.

`run_replicas.sh` corre `COMMANDS.sh` completo en la primera réplica, copia su `eq5.rst7` a las demás y en ellas ejecuta solo la parte `md`. La minimización usa `mpirun -np 10 pmemd.MPI` (`AMBER_MIN_EXE` en `settings.cfg`; ajusta el número de procesos a tu equipo).

## 2. Caja de solvente con agua OPC

Las cajas incluidas (`ETAWAT20` y las demás) están equilibradas con TIP3P. Para usar OPC hay que reconstruirlas y reequilibrarlas con ese modelo de agua. No hace falta reparametrizar los cosolventes: las cargas RESP y los parámetros GAFF del cosolvente no dependen del agua. Sí conviene validar densidad y agregación (sección 2.3).

### 2.1 Construcción

El script `examples/solvent_boxes/build_box.sh` construye la caja de cualquiera de los cosolventes de la biblioteca (ETA, ANT, MAM, ISO, ISO5, MOH, PYR1, ION). El procedimiento completo, la calculadora de composiciones (`boxcalc.py`), la validación y los problemas frecuentes están en [TUTORIAL_CAJAS.md](TUTORIAL_CAJAS.md). Para etanol:

```bash
mkdir -p ~/cajas/ETA && cd ~/cajas/ETA
nohup bash /ruta/a/pyMDmix/examples/solvent_boxes/build_box.sh ETA > build.log 2>&1 &
tail -f build.log
```

Con `REP=2` (el valor que el script elige para etanol) son 136 ETA y 1672 aguas en unos 40 Å, y los 10 ns de NPT toman pocos minutos en GPU. **Verificar**: que `6_analysis.log` muestre una deriva de densidad cercana a cero entre la primera y la segunda mitad, y que `best.pdb` contenga tantos átomos `EPW` como aguas (1672 con `REP=2`).

### 2.2 Registrar el solvente

```bash
cp ETAOPC20.off /ruta/de/trabajo/ && cd /ruta/de/trabajo
cp /ruta/a/pyMDmix/pyMDMix/data/solventlib/eta.config eta_opc.config
```

Edita `eta_opc.config`:

```ini
[GENERAL]
name = ETAO
info = Ethanol 20%% mixture in OPC water
objectfile = ETAOPC20.off
watermodel = OPC
boxunit = ETAOPC20

[PROBES]
WAT=WAT@O
CT=ETA@C1
OH=ETA@O1

[TYPES]
OH=Don,Acc
CT=Hyd
WAT=Wat
```

Luego:

```bash
mdmix create solvents -f eta_opc.config
mdmix info solvents | grep ETAO
```

Usa un nombre distinto de `ETA` para no sobrescribir la mezcla TIP3P. El archivo `examples/solvent_boxes/opc.config` es el equivalente para agua OPC pura.

### 2.3 Validación recomendada

Antes de usar la caja en producción, compara la densidad media de `6_analysis.log` con el valor experimental de la mezcla y revisa que no haya agregación anormal del cosolvente (función de distribución radial cosolvente-cosolvente). El equilibrio entre interacciones solvente-agua y solvente-solvente es lo que más puede cambiar al pasar de TIP3P a OPC.

## 3. Crear un proyecto

Genera una plantilla de configuración y edítala:

```bash
mdmix create template -f proyecto.cfg
```

La plantilla completa y comentada está en `pyMDMix/data/templates/complete.cfg`. Los puntos esenciales:

```ini
[SYSTEM]
NAME = EF1
OFF = ef1.off
EXTRARES = ATP,ADP
# Con agua OPC es obligatorio cargar su campo de fuerza
EXTRAFF = leaprc.water.opc

[MDSETTINGS]
SOLVENTS = ETAO
NREPL = 5
NANOS = 20
TEMP = 310
RESTR = HA
FORCE = 0.1
```

Crea el proyecto:

```bash
mdmix create project -n EF1_S2 -f proyecto.cfg
cd EF1_S2/MD
```

Las carpetas de réplicas quedan dentro de `MD`, con el nombre `SOLVENTE_N` (por ejemplo `ETAO_1`). Los archivos de topología y coordenadas se guardan ahora como `*.parm7` y `*.rst7`.

## 4. Ajustar el paso de integración y la frecuencia

No es necesario editar los archivos `.in` con `sed`. Todas esas modificaciones se declaran en la sección `[MDSETTINGS]` del `cfg` antes de crear el proyecto:

| Cambio en los `.in` | Opción equivalente en `[MDSETTINGS]` |
|---|---|
| `dt=0.002` a `dt=0.004` | `MD_TIMESTEP = 4` (en fs) |
| `ntwx` y `ntpr` de 500 a 250 | `TRAJFREQUENCY = 250` |
| `nstlim=500000` a `250000` en `md.in` | `PROD_STEPS = 250000` |
| `nstlim=500000` a `250000` en `eq5.in` | `NPT_EQ_STEPS = 250000` |
| `nstlim= 100000` a `50000` (calentamiento) | `HEATING_STEPS = 50000` |
| `iwrap=1` a `iwrap=0` | `IWRAP = 0` |

Ejemplo para pasos de 4 fs que conservan 1 ns por archivo y 1 ps entre cuadros:

```ini
MD_TIMESTEP = 4
TRAJFREQUENCY = 250
HEATING_STEPS = 50000
NPT_EQ_STEPS = 250000
PROD_STEPS = 250000
IWRAP = 0
```

Un paso de 4 fs solo es estable si la topología tiene reparto de masa de hidrógenos (HMR). pyMDMix no lo aplica; hay que hacerlo antes sobre el `parm7` (por ejemplo con ParmEd, `HMassRepartition`). **Verificar** con una corrida corta que la energía total sea estable.

Si ya creaste el proyecto y prefieres editar los archivos a mano, los `.in` están en `*/min`, `*/eq` y `*/md` de cada réplica.

## 5. Ejecutar

Genera los archivos para tu sistema de colas (o revisa `mdmix queue list`):

```bash
mdmix queue list
mdmix queue write -n NOMBRE_COLA
```

Cada réplica tiene `min/`, `eq/` y `md/`. Los pasos de equilibración van de `eq1` a `eq5`, y la producción se divide en archivos `md1`, `md2`, etcétera.

## 6. Análisis

Se ejecuta desde la carpeta `MD` del proyecto, cuando las simulaciones han avanzado. La máscara es la lista de residuos del soluto sobre cuyo esqueleto se alinea.

```bash
MASK=331,332,335,336,337,343,344,371,372,374,375,378,404,405

for i in 1 2 3 4 5; do
    mdmix analyze align byname -s ETAO_$i -C4 --mask $MASK
done

mdmix plot rmsd bysolvent -s ETAO -o etao_rmsd.png
mdmix analyze density bysolvent -s ETAO -C4
mdmix analyze energy bysolvent -s ETAO
```

El alineamiento procesa los pasos de producción disponibles. Si lo lanzas mientras corren las simulaciones, se puede limitar con `-N 1:5` (del primer al quinto paso).

Los hot spots se crean a partir de los mapas de energía libre promediados:

```bash
mdmix analyze hotspots create \
    -i PROBE_AVG/ETAO_CT_DG0.dx PROBE_AVG/ETAO_OH_DG0.dx PROBE_AVG/ETAO_WAT_DG.dx \
    -o ETAO_allprobes_hotspots
```

Salen un PDB con los hot spots y un archivo `.hset` (pickle) con el conjunto.

## 7. Problemas frecuentes

| Síntoma | Causa y solución |
|---|---|
| `Could not open file leaprc.ff14SB` | Ese archivo ya no existe en AmberTools recientes. Usa `leaprc.protein.ff19SB` (valor por defecto en esta versión). |
| `Small box detected` o `max pairlist cutoff must be less than unit cell max sphere radius` | La caja es demasiado pequeña para el corte de 9 Å. Usa una caja mayor (`REP=2` en el script de cajas). |
| `mdmix` ejecuta un script de Python 2 | Hay otro `mdmix` en `~/.local/bin`. Activa el entorno `pymdmix3` o usa `python -m pyMDMix`. |
| `AMBERHOME env variable not defined!` | Ejecuta `source $AMBERHOME/amber.sh` antes de usar `mdmix`. |
| Conflictos de `numpy` con herramientas de AmberTools | Fija `numpy<2` y `biopython<1.86` (ya incluido en `environment.yml`). |
| Faltan los átomos `EPW` en un PDB escrito por cpptraj | Añade `include_ep` al comando `trajout ... pdb` (ya incluido en `check_com`). |
