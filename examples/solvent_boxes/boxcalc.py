#!/usr/bin/env python3
"""
boxcalc.py: composition calculator for mixed-solvent boxes.

Replaces the molecule_calculator.ods spreadsheet. All calculations use the
molecular volumes of the pure components (V = M / (rho * NA)), as the
spreadsheet did. The real volume of the mixture differs by a few percent (excess
volume); that is why the box is later equilibrated in NPT and the final density
is validated against experiment.

Subcommands
-----------
  table                              list of cosolvents (density, molar mass)
  pure   COSOLV --side 23            molecules in a cubic box of the pure cosolvent
  mix    COSOLV --vv 20 --side 20    cosolvent and water molecules for a given % v/v
  check  COSOLV --nco 17 --nwat 209 --side 20
                                     actual concentration of a given composition
  molar  COSOLV --molarity 7 --volume 0.1
                                     grams and milliliters for a solution of a given molarity
  dilute COSOLV --molar 0.25 --side 40
                                     molecules for a low molar concentration (apolar probes)
  library                            actual concentration of the boxes included in pyMDMix

Examples
--------
  python boxcalc.py mix ETA --vv 20 --side 20
  python boxcalc.py mix ETA --vv 20 --side 21 --rep 2
  python boxcalc.py check ETA --nco 17 --nwat 209 --side 20
  python boxcalc.py molar UREA --molarity 7 --volume 0.1
"""
import argparse
import math
import sys

NA = 6.02214076e23          # 1/mol
RHO_WAT = 1.0               # g/mL (the spreadsheet uses 1.0 and 18.02 g/mol)
MW_WAT = 18.02              # g/mol
V_WAT = MW_WAT / (RHO_WAT * NA) * 1e24   # A^3 per water molecule (29.92)

# short code: (name, density g/mL, molar mass g/mol, source)
# "hoja" = values from molecule_calculator.ods. "lit" = literature values
# added for the cosolvents of the included boxes; confirm with your
# supplier before using them in published work.
COSOLVENTS = {
    'ETA':  ('Ethanol',                        0.789, 46.07,   'hoja'),
    'TFE':  ('Trifluoroethanol',               1.39,  100.04,  'hoja'),
    'DMSO': ('Dimethyl sulfoxide',             1.1,   78.13,   'hoja'),
    'IMID': ('Imidazole',                      1.23,  68.077,  'hoja'),
    'PYRM': ('Pyrimidine',                     1.02,  80.088,  'hoja'),
    'PHEN': ('Phenol',                         1.07,  94.113,  'hoja'),
    'BENZ': ('Benzene',                        0.876, 78.11,   'hoja'),
    'ACET': ('Acetone',                        0.784, 58.08,   'hoja'),
    'ACOH': ('Acetic acid',                    1.05,  60.052,  'hoja'),
    'PROP': ('Propionic acid',                 0.99,  74.08,   'hoja'),
    'MPD':  ('2-Methyl-1,3-propanediol',       1.015, 90.12,   'hoja'),
    'DMAC': ('Dimethylacetylcarbinol',         1.1,   102.132, 'hoja'),
    'HQ':   ('Hydroquinone',                   1.33,  110.111, 'hoja'),
    'HBA':  ('4-hydroxybenzoic acid',          1.46,  138.12,  'hoja'),
    'MPP':  ('2-(2-Methoxypropoxy)-1-propanol', 1.0,  148.2,   'hoja'),
    'HMP':  ('4-hydroxy-6-methyl-2-pyrone',    1.348, 126.12,  'hoja'),
    'BUOH': ('Butanol',                        0.81,  73.123,  'hoja'),
    'PYRN': ('Pyran-4-one',                    1.2,   96.08,   'hoja'),
    'MALT': ('Maltol',                         1.348, 126.111, 'hoja'),
    'EMAL': ('Ethyl maltol',                   1.1624, 140.1366, 'hoja'),
    'ABUT': ('3-amino-1-butanol',              0.927, 89.136,  'hoja'),
    'UREA': ('Urea',                           1.32,  60.06,   'hoja'),
    # cosolvents of the boxes included in pyMDMix (literature values)
    'MOH':  ('Methanol',                       0.792, 32.04,   'lit'),
    'ISO':  ('Isopropanol',                    0.786, 60.10,   'lit'),
    'ANT':  ('Acetonitrile',                   0.786, 41.05,   'lit'),
    'MAM':  ('Acetamide',                      1.159, 59.07,   'lit'),
    'PYR':  ('Pyridine',                       0.982, 79.10,   'lit'),
    'IBN':  ('Isobutane (liquido, ~0 C)',      0.551, 58.12,   'lit'),
    'ACD':  ('Acetaldehyde',                   0.784, 44.05,   'lit'),
}

# Boxes included in pyMDMix/data/solventlib: (cosolvent, n cosolvent, n water, side A)
# read from the .off files. ION (ammonium acetate, charged) is not treated as a neutral cosolvent.
LIBRARY = {
    'ETAWAT20':  ('ETA', 17, 209, 19.990),
    'ANTWAT20':  ('ANT', 17, 209, 19.861),
    'MAMWAT20':  ('MAM', 18, 209, 19.847),
    'ISOWAT20':  ('ISO', 14, 209, 20.059),
    'ISOWAT5':   ('ISO', 4,  310, 21.407),
    'MOHWAT20':  ('MOH', 56, 470, 26.226),
    'PYRWAT1':   ('PYR', 1,  438, 27.99),
}


def get(code):
    key = code.upper()
    if key not in COSOLVENTS:
        sys.exit("Cosolvente '%s' no esta en la tabla. Usa 'table' para ver las claves." % code)
    return COSOLVENTS[key]


def v_mol(rho, mw):
    """Molecular volume (A^3) of a pure liquid."""
    return mw / (rho * NA) * 1e24


def composition(nco, nwat, side, rho, mw):
    """Concentrations of a box with nco cosolvents and nwat waters."""
    vbox = side ** 3
    mass_co = nco * mw / NA
    mass_w = nwat * MW_WAT / NA
    vco = nco * v_mol(rho, mw)
    vw = nwat * V_WAT
    return {
        'molar': nco / (NA * vbox * 1e-27),                 # mol/L of box
        'x_co': nco / float(nco + nwat),                    # mole fraction
        'vv_ideal': 100.0 * vco / (vco + vw),               # % v/v, pure volumes
        'ww': 100.0 * mass_co / (mass_co + mass_w),         # % p/p
        'dens_ini': (mass_co + mass_w) * 1e24 / vbox,       # g/mL in the initial box
        'vbox': vbox,
    }


def show(nco, nwat, side, rho, mw, label=''):
    c = composition(nco, nwat, side, rho, mw)
    print('%s%d cosolvente + %d agua | caja %.2f A (%.0f A^3)' % (label, nco, nwat, side, c['vbox']))
    print('   %% v/v (volumenes puros): %.1f | %% p/p: %.1f | fraccion molar: %.4f | %.2f M | densidad de la caja inicial: %.3f g/mL'
          % (c['vv_ideal'], c['ww'], c['x_co'], c['molar'], c['dens_ini']))


def cmd_table(a):
    print('%-5s %-34s %8s %10s  %s' % ('clave', 'nombre', 'g/mL', 'g/mol', 'origen'))
    for k, (n, r, m, o) in COSOLVENTS.items():
        print('%-5s %-34s %8.3f %10.3f  %s' % (k, n, r, m, o))


def cmd_pure(a):
    _, rho, mw, _ = get(a.cosolv)
    vbox = a.side ** 3
    n = rho * vbox / 1e24 / mw * NA
    print('Caja pura de %s: %.1f A -> %.1f moleculas (%d redondeado)' % (a.cosolv.upper(), a.side, n, round(n)))


def cmd_mix(a):
    _, rho, mw, _ = get(a.cosolv)
    side = a.side * a.rep
    phi = a.vv / 100.0
    vbox = side ** 3
    nco = vbox * phi / v_mol(rho, mw)
    nwat = vbox * (1.0 - phi) / V_WAT
    print('Objetivo: %.1f %% v/v de %s en agua, caja %.2f A' % (a.vv, a.cosolv.upper(), side))
    print('   cosolvente: %.2f -> %d | agua: %.2f -> %d' % (nco, round(nco), nwat, round(nwat)))
    show(round(nco), round(nwat), side, rho, mw, label='Composicion redondeada: ')
    need = math.ceil(36.0 / a.side)
    if a.rep < need:
        print('AVISO: con corte de 9 A conviene una caja de al menos unos 36 A de lado. '
              'Para el lado %.1f A usa --rep %d (o replica con REP=%d en build_box.sh).' % (a.side, need, need))
    print('Para build_box.sh: NCO=%d NWAT=%d SIDE0=%.2f' % (round(nco), round(nwat), side))


def cmd_check(a):
    _, rho, mw, _ = get(a.cosolv)
    show(a.nco, a.nwat, a.side, rho, mw)


def cmd_molar(a):
    _, rho, mw, _ = get(a.cosolv) if a.cosolv else ('', a.density, a.mw, '')
    if a.density:
        rho = a.density
    if a.mw:
        mw = a.mw
    moles = a.molarity * a.volume
    grams = moles * mw
    ml = grams / rho
    print('%.3f M en %.3f L: %.4f mol = %.3f g' % (a.molarity, a.volume, moles, grams))
    print('Volumen del cosolvente puro: %.2f mL (%.1f %% v/v si fuera aditivo)' % (ml, 100.0 * ml / (a.volume * 1000.0)))


def cmd_dilute(a):
    """Box with one cosolvent at low molar concentration: n = M * V; the rest is water."""
    _, rho, mw, _ = get(a.cosolv)
    side = a.side * a.rep
    vbox = side ** 3
    nco = a.molar * NA * vbox * 1e-27
    vco = nco * v_mol(rho, mw)
    nwat = (vbox - vco) / V_WAT
    nco_i, nwat_i = max(1, round(nco)), round(nwat)
    print('Objetivo: %.3f M de %s en una caja de %.2f A (%.0f A^3)' % (a.molar, a.cosolv.upper(), side, vbox))
    print('   cosolvente: %.2f -> %d | agua: %.1f -> %d' % (nco, nco_i, nwat, nwat_i))
    show(nco_i, nwat_i, side, rho, mw, label='Composicion redondeada: ')
    print('Para build_box.sh: COMP="NOMBRE:%d" NW0=%d SIDE_EQ=%.2f REP=1' % (nco_i, nwat_i, side))


def cmd_library(a):
    print('Concentracion real de las cajas incluidas (volumenes de componentes puros)')
    for name, (code, nco, nwat, side) in LIBRARY.items():
        _, rho, mw, _ = COSOLVENTS[code]
        show(nco, nwat, side, rho, mw, label='%-9s ' % name)
    print('ION (acetato de amonio, 12 pares + 315 aguas, 22.17 A): sal cargada, no aplica %v/v.')
    print('Densidad de la caja inicial = valor antes de equilibrar; la final sale del NPT.')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    sub.add_parser('table').set_defaults(f=cmd_table)

    s = sub.add_parser('pure'); s.add_argument('cosolv'); s.add_argument('--side', type=float, required=True)
    s.set_defaults(f=cmd_pure)

    s = sub.add_parser('mix'); s.add_argument('cosolv')
    s.add_argument('--vv', type=float, required=True, help='cosolvent % v/v')
    s.add_argument('--side', type=float, required=True, help='side of the base box (A)')
    s.add_argument('--rep', type=int, default=1, help='replicas per side (final side = side*rep)')
    s.set_defaults(f=cmd_mix)

    s = sub.add_parser('check'); s.add_argument('cosolv')
    s.add_argument('--nco', type=int, required=True); s.add_argument('--nwat', type=int, required=True)
    s.add_argument('--side', type=float, required=True)
    s.set_defaults(f=cmd_check)

    s = sub.add_parser('molar'); s.add_argument('cosolv', nargs='?')
    s.add_argument('--molarity', type=float, required=True); s.add_argument('--volume', type=float, required=True, help='liters')
    s.add_argument('--density', type=float); s.add_argument('--mw', type=float)
    s.set_defaults(f=cmd_molar)

    s = sub.add_parser('dilute'); s.add_argument('cosolv')
    s.add_argument('--molar', type=float, required=True, help='desired concentration (mol/L)')
    s.add_argument('--side', type=float, required=True, help='box side (A)')
    s.add_argument('--rep', type=int, default=1)
    s.set_defaults(f=cmd_dilute)

    sub.add_parser('library').set_defaults(f=cmd_library)

    a = p.parse_args()
    a.f(a)


if __name__ == '__main__':
    main()
