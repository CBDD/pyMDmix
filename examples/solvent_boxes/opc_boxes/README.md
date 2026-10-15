# Cajas de solvente OPC validadas

Cada caja se construyó con `build_box.sh`, se equilibró 10 ns en NPT a 300 K y se validó (densidad, deriva y RDF de agregación). Ver `TUTORIAL_CAJAS.md`.

| Solvente | Mezcla | Composición | Densidad (g/cm3) | Lado (A) |
|---|---|---|---|---|
| ETAO | etanol | 136 ETA + 1672 WAT | 0.9678 | 39.67 |
| ANTO | acetonitrilo | 136 ANT + 1672 WAT | 0.9565 | 39.58 |
| MOHO | metanol | 448 MOH + 3760 WAT | 0.9682 | 52.02 |
| ISOO | isopropanol 20 % | 112 ISO + 1672 WAT | 0.9696 | 39.81 |
| ISO5O | isopropanol 5 % | 32 ISO + 2480 WAT | 0.9890 | 42.77 |
| MAMO | acetamida | 144 MAM + 1672 WAT | 1.0185 | 39.79 |
| IONO | acetato + metilamonio | 96 CN3 + 96 COO + 2520 WAT | 1.0472 | 44.12 |
| PYR1O | piridina (~0.13 M) | 8 PYR + 3504 WAT | 0.9957 | 47.38 |

Registro en la base de datos de solventes (entorno `pymdmix3` activo):

```bash
bash register_opc_boxes.sh
```

Cada `.config` usa `watermodel = OPC` y rutas relativas al `.off` de esta carpeta. Cargan los tipos GAFF de los cosolventes originales de pyMDMix; con `leaprc.gaff2` los parametros de enlace pueden diferir ligeramente de los usados al equilibrar la caja.
