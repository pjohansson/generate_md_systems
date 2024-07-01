#!/usr/bin/env python3

import datetime
import numpy as np
import os
import subprocess

from generate_md_systems.gmx_conf_utils import *

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory
from typing import Sequence


def scale_atom_coords(
    atoms: Sequence[Atom],
    scale_x: float,
    scale_y: float,
    scale_z: float,
    translate: Vec3 | None,
    residue_length: int,
) -> list[Atom]:
    if residue_length != 1:
        raise ValueError(
            "error: residue_length != 1 not implemented (scale_atom_coords)"
        )

    if translate is None:
        dx = 0.
        dy = 0.
        dz = 0.
    else:
        dx, dy, dz = translate

    buf = []

    for atom in atoms:
        x0, y0, z0 = atom.position

        x1 = dx + (scale_x * x0)
        y1 = dy + (scale_y * y0)
        z1 = dz + (scale_z * z0)

        buf.append(Atom(
            position=Vec3(x1, y1, z1),
            velocity=atom.velocity,
            name=atom.name,
            residue=atom.residue,
        ))

    return buf


def create_gromos87_conf_with_density(
    conf: Gromos87,
    box_x: float,
    box_y: float,
    box_z: float,
    density: float,
    residue_length: int = 1,
    translate: Vec3 | None = None,
) -> Gromos87:
    num_final = int(np.round(density * box_x * box_y * box_z))
    num_conf = len(conf.atoms)

    # We start with a box with the current density but cut it
    # to a size which contains slightly more than the desired
    # number of atoms. We then remove a random selection of atoms
    # from it until we reach the correct number, and finally
    # scale the coordinates to the correct density.
    scale_vol = (num_final / num_conf)**(1./3.)
    box_x_in, box_y_in, box_z_in = conf.box_size

    while True:
        scale_vol *= 1.01
        conf_scaled = create_gromos87_conf_with_size(
            conf,
            scale_vol * box_x_in,
            scale_vol * box_y_in,
            scale_vol * box_z_in,
            residue_length=residue_length,
        )

        if len(conf_scaled.atoms) >= num_final:
            break

    np.random.shuffle(conf_scaled.atoms)
    atoms_unscaled = conf_scaled.atoms[:num_final]

    scale_x = box_x / conf_scaled.box_size[0]
    scale_y = box_y / conf_scaled.box_size[1]
    scale_z = box_z / conf_scaled.box_size[2]

    atoms = scale_atom_coords(
        atoms_unscaled,
        scale_x,
        scale_y,
        scale_z,
        translate=translate,
        residue_length=residue_length,
    )

    return Gromos87(
        title=conf.title,
        atoms=atoms,
        box_size=Vec3(box_x, box_y, box_z),
    )


def main():
    parser = ArgumentParser(
        description="""
            Create a system with a fixed size from an input .gro file.
            """
    )

    parser.add_argument(
        'input',
        type=str, metavar='PATH',
        help="path to original configuration")
    parser.add_argument(
        'x',
        type=float,
        help="size of system along x")
    parser.add_argument(
        'y',
        type=float,
        help="size of system along y")
    parser.add_argument(
        'z',
        type=float,
        help="size of system along z")

    parser.add_argument(
        '-l', '--residue-length',
        default=1, type=int, metavar='N',
        help="number of atoms per molecule (residue) (default: %(default)s)",
    )

    parser.add_argument(
        '-d', '--density',
        default=None, type=float, metavar='RHO',
        help="number density of output configuration (default: none)",
    )
    parser.add_argument(
        '--translate',
        default=None, nargs=3, metavar=('DX', 'DY', 'DZ'), type=float,
        help="translate final configuration")

    parser_output = parser.add_argument_group('output options')
    parser_output.add_argument(
        '-o', '--output',
        type=str, default=None, metavar='PATH',
        help="output path for final system configuration (default: stdout)")
    parser_output.add_argument(
        '-t', '--topology',
        type=str, metavar='PATH', default=None,
        help="optionally write topology `[molecules]` directive to this path")
    parser_output.add_argument(
        '--title',
        type=str, default="Gromacs configuration",
        help="set title for output configuration")

    args = parser.parse_args()

    conf = read_gromos87(args.input)

    if args.density is None:
        conf_final = create_gromos87_conf_with_size(
            conf,
            args.x, args.y, args.z,
            residue_length=args.residue_length,
            translate=args.translate,
        )
    else:
        conf_final = create_gromos87_conf_with_density(
            conf,
            args.x, args.y, args.z,
            args.density,
            residue_length=args.residue_length,
            translate=args.translate,
        )

    if args.output is None:
        write_gromos87_stdout(conf_final)
    else:
        write_gromos87(args.output, conf_final)
