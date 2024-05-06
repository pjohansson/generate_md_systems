import numpy as np

from argparse import ArgumentParser
from typing import Sequence
from sys import exit

from generate_md_systems.gmx_conf_utils import (
    Atom,
    Gromos87,
    Vec3,
    read_gromos87,
    write_gromos87,
    write_gromos87_stdout,
)


# DEFAULT_RESIDUE_LENGTHS = {
#     'SOL': 3,  # 3-point water
# }


def cut_atoms(
    atoms: Sequence["Atom"],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    zlim: tuple[float, float],
) -> list["Atom"]:
    def parse_lims(lims: tuple[float, float]) -> tuple[float, float]:
        v0, v1 = lims

        if v0 < 0.:
            v0 = -np.inf
        if v1 < 0.:
            v1 = +np.inf

        return v0, v1

    def atom_in_lims(atom: "Atom"):
        x, y, z = atom.position

        return (
            (x >= x0) and (x <= x1)
            and (y >= y0) and (y <= y1)
            and (z >= z0) and (z <= z1)
        )

    x0, x1 = parse_lims(xlim)
    y0, y1 = parse_lims(ylim)
    z0, z1 = parse_lims(zlim)

    # num_left_in_res = None

    atoms_inside = []

    for atom in atoms:

        if atom_in_lims(atom):
            atoms_inside.append(atom)

    return atoms_inside


def shift_atoms_to_origin(atoms: Sequence["Atom"]) -> list["Atom"]:
    if len(atoms) == 0:
        return []

    xmin, ymin, zmin = atoms[0].position

    for atom in atoms[1:]:
        x, y, z = atom.position

        if x < xmin:
            xmin = x
        if y < ymin:
            ymin = y
        if z < zmin:
            zmin = z

    atoms_shifted = []

    for atom in atoms:
        x, y, z = atom.position

        position = Vec3(x - xmin, y - ymin, z - zmin)

        atoms_shifted.append(
            Atom(
                position=position,
                velocity=atom.velocity,
                name=atom.name,
                residue=atom.residue,
            )
        )

    return atoms_shifted


def calc_box_size(atoms: Sequence["Atom"]) -> Vec3:
    if len(atoms) == 0:
        raise ValueError("no atoms in list: cannot calculate box size")

    xmax, ymax, zmax = atoms[0].position

    for atom in atoms[1:]:
        x, y, z = atom.position

        if x > xmax:
            xmax = x
        if y > ymax:
            ymax = y
        if z > zmax:
            zmax = z

    return Vec3(xmax, ymax, zmax)


def main():
    parser = ArgumentParser()

    parser.add_argument(
        'conf',
        help='.gro configuration to cut from',
    )

    parser.add_argument(
        '-o', '--output',
        default=None, type=str, metavar='PATH',
        help='write cut configuration to path (default: stdout)',
    )
    parser.add_argument(
        '-s', '--shift-to-origin',
        action='store_true',
        help='shift remaining atoms to origin and resize box',
    )

    parser_cut = parser.add_argument_group('cut options')
    parser_cut.add_argument(
        '-x', '--xlim',
        nargs=2, type=float,
        default=(-1., -1.), metavar=('XMIN', 'XMAX'),
        help='limits along the x axis (<0 = no limit)',
    )
    parser_cut.add_argument(
        '-y', '--ylim',
        nargs=2, type=float,
        default=(-1., -1.), metavar=('YMIN', 'YMAX'),
        help='limits along the y axis (<0 = no limit)',
    )
    parser_cut.add_argument(
        '-z', '--zlim',
        nargs=2, type=float,
        default=(-1., -1.), metavar=('ZMIN', 'ZMAX'),
        help='limits along the z axis (<0 = no limit)',
    )

    args = parser.parse_args()

    try:
        conf = read_gromos87(args.conf)
    except Exception as exc:
        print(f'error: could not read configuration ({exc})')
        exit(1)

    atoms = cut_atoms(conf.atoms, args.xlim, args.ylim, args.zlim)

    if args.shift_to_origin:
        atoms = shift_atoms_to_origin(atoms)
        box_size = calc_box_size(atoms)
    else:
        box_size = conf.box_size

    conf = Gromos87(title=conf.title, atoms=atoms, box_size=box_size)

    if args.output is None:
        write_gromos87_stdout(conf)
    else:
        write_gromos87(args.output, conf)
