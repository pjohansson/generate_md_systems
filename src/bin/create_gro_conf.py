#!/usr/bin/env python3

import datetime
import numpy as np
import os
import subprocess

from generate_md_systems.gmx_conf_utils import *

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory


def stack_conf_to_minimum_size(conf, to_size_x, to_size_y, to_size_z):
    """Stack multiples of the configuration to a minimum final size.
    
    Note that this only stacks the boxes, it does not cut the 
    configuration to the correct size. It only guarantuees that
    the output size is at least this large.
    
    """

    def move_atom(atom, ix, iy, iz, dx, dy, dz):
        x = atom.position.x + ix * dx
        y = atom.position.y + iy * dy
        z = atom.position.z + iz * dz

        return Atom(
                name=atom.name,
                residue=atom.residue,
                position=Vec3(x, y, z),
                velocity=atom.velocity,
                )

    size_x, size_y, size_z = conf.box_size

    nx = int(np.ceil(to_size_x / size_x))
    ny = int(np.ceil(to_size_y / size_y))
    nz = int(np.ceil(to_size_z / size_z))

    atoms = conf.atoms.copy()

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if ix > 0 or iy > 0 or iz > 0:
                    atoms += [
                        move_atom(atom, ix, iy, iz, size_x, size_y, size_z) 
                        for atom in conf.atoms]

    box_size = size_x * nx, size_y * ny, size_z * nz

    return Gromos87(
            title=conf.title,
            box_size=box_size,
            atoms=atoms,
            )


def cut_conf_to_size(conf, xmax, ymax, zmax, residue_length=1):
    """Cut a configuration to the set size.
    
    Molecules are kept if any of their atoms lie inside the box.

    """

    def check_atom(atom, xmax, ymax, zmax):
        return (atom.position.x <= xmax 
                and atom.position.y <= ymax
                and atom.position.z <= zmax)

    # Keep molecule if any atom is inside, else discard
    def check_molecule(atoms, xmax, ymax, zmax):
        return np.any(
            [check_atom(atom, xmax, ymax, zmax) for atom in atoms])

    num_mols = len(conf.atoms) // residue_length
    kept_atoms = []

    for i in range(num_mols):
        i0 = i * residue_length
        i1 = (i + 1) * residue_length

        atoms = conf.atoms[i0:i1]

        if check_molecule(atoms, xmax, ymax, zmax):
            kept_atoms += atoms

    box_size = xmax, ymax, zmax

    return Gromos87(
            title=conf.title,
            box_size=box_size,
            atoms=kept_atoms
            )


def create_conf_with_size(conf, size_x, size_y, size_z, residue_length=1, translate=None):
    """Expand or cut the given configuration to the set size."""

    conf_stacked = stack_conf_to_minimum_size(conf, size_x, size_y, size_z)
    conf_resized = cut_conf_to_size(conf_stacked, size_x, size_y, size_z, residue_length=residue_length)
    
    if translate != None:
        conf_final = translate_gromos87(conf_resized, *translate)
    else:
        conf_final = conf_resized

    return conf_final


if __name__ == '__main__':
    parser = ArgumentParser(
        description="""
            Create a system with a fixed size from an input .gro file.
            """
    )

    parser.add_argument('input', 
        type=str, metavar='PATH',
        help="path to original configuration")
    parser.add_argument('x', 
        type=float,
        help="size of system along x")
    parser.add_argument('y', 
        type=float,
        help="size of system along y")
    parser.add_argument('z', 
        type=float,
        help="size of system along z")

    parser.add_argument('--translate',
        default=None, nargs=3, metavar=('DX', 'DY', 'DZ'), type=float,
        help="translate final configuration")

    parser_output = parser.add_argument_group('output options')
    parser_output.add_argument('-o', '--output', 
        type=str, default=None, metavar='PATH',
        help="output path for final system configuration (default: print to stdout)")
    parser_output.add_argument('-t', '--topology', 
        type=str, metavar='PATH', default=None, 
        help="optionally write topology `[ molecules ]` directive to this path")
    parser_output.add_argument('--title', 
        type=str, default="Lennard-Jones slab system", 
        help="set title for output configuration")

    args = parser.parse_args()

    conf = read_gromos87(args.input)
    conf_final = create_conf_with_size(
        conf, args.x, args.y, args.z, translate=args.translate)
    
    if args.output == None:
        write_gromos87_stdout(conf_final)
    else:
        write_gromos87(args.output, conf_final)
    
