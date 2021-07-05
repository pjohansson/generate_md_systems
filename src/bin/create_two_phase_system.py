#!/usr/bin/env python3

import two_phase
import numpy as np
import os
import subprocess

from argparse import ArgumentParser
from gmxutils import *
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory

# Surfactant molecule area density in interface (used in other simulations)
# sur_area_density = (344 / 2) / (40. * 20.)
sur_residue_length = 2


def verify_conf_is_monotype(conf):
    """Asserts that all residues in the configuration are identical.
    
    Raises an `AssertionError` otherwise.
    
    """

    head_residue = conf.atoms[0].residue

    for atom in conf.atoms[1:]:
        assert atom.residue == head_residue


def read_conf(path, default_path, default_dir='include'):
    """Read a configuration either from a given path or the database.
    
    If the input `path` is `None` the `default_path` is used to access
    the files included with the module `two_phase`.

    """

    if not path:
        for module_dir in two_phase.__path__:
            path = os.path.join(module_dir, default_dir, default_path)

            if os.path.exists(path):
                break

        try:
            assert os.path.exists(path)
        except AssertionError as exc:
            stderr.write(
                "error: could not read included data file '{}', "
                "was it installed correctly with the module? "
                "(full path: '{}')\n".format(default_path, path))
            exit(1)
    
    conf = read_gromos87(path)

    try:
        verify_conf_is_monotype(conf)
    except AssertionError:
        stderr.write(
            "error: configuration '{}' includes more than "
            "one residue type\n".format(path))
        exit(1)

    return conf


def calc_volume(conf):
    """Calculate the box volume."""

    dx, dy, dz = conf.box_size
    volume = dx * dy * dz

    return volume


def calc_density(conf, residue_length):
    """Calculate the molecule number density inside the configuration."""

    num_atoms = len(conf.atoms)

    assert num_atoms % residue_length == 0

    num_mols = num_atoms // residue_length

    return num_mols / calc_volume(conf)


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


def cut_conf_to_size(conf, xmax, ymax, zmax, residue_length):
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


def correct_conf_density(conf, residue_length):
    """Correct the configuration to the set molecule density."""

    density = calc_density(conf, residue_length)
    volume = calc_volume(conf)

    target_num_mols = int(density * volume)
    num_mols = len(conf.atoms) // residue_length

    if target_num_mols != num_mols:
        stderr.write("warning: created number of molecules ({}) differs from the target number of molecules ({})\n".format(num_mols, target_num_mols))

    return conf


def create_conf_with_size(conf, size_x, size_y, size_z, residue_length):
    """Expand or cut the given configuration to the set size."""

    conf_stacked = stack_conf_to_minimum_size(conf, size_x, size_y, size_z)
    conf_resized = cut_conf_to_size(conf_stacked, size_x, size_y, size_z, residue_length)
    conf_final = correct_conf_density(conf_resized, residue_length)

    return conf_final


def create_stack(conf_one, conf_two, separation):
    """Translate and stack the two phases with the given separation."""

    def translate_atoms(atoms, dz):
        def move_atom(atom, dz):
            position = Vec3(x=atom.position.x, y=atom.position.y, z=atom.position.z + dz)

            return Atom(
                    name=atom.name,
                    residue=atom.residue,
                    position=position,
                    velocity=atom.velocity
                    )

        return [move_atom(atom, dz) for atom in atoms]

    size_x, size_y, size_z_one = conf_one.box_size
    _, _, size_z_two = conf_two.box_size

    dz_one = size_z_two / 2. + separation
    dz_two = -size_z_two / 2.

    size_z_total = size_z_one + size_z_two + 2. * separation

    conf_one_atoms = translate_atoms(conf_one.atoms, dz_one)
    conf_two_atoms = translate_atoms(conf_two.atoms, dz_two)

    atoms = conf_one_atoms + conf_two_atoms

    return Gromos87(
            title=conf_one.title,
            box_size=(size_x, size_y, size_z_total),
            atoms=atoms,
            )


def get_surfactant_conf(conf_one, conf_two, 
                        separation, conf_stacked, sur_area_density):
    """Generate surfactants at both interfaces."""

    _, _, size_z_one = conf_one.box_size
    _, _, size_z_two = conf_two.box_size
    size_x, size_y, size_z_final = conf_stacked.box_size

    z0 = (size_z_two + separation) / 2.
    z1 = z0 + size_z_one + separation

    area = size_x * size_y
    num_mols = int(sur_area_density * area)

    binary = 'gen_surfactant'
    args = [
            binary,
            '{}'.format(num_mols),
            '-b',
            '{}'.format(size_x),
            '{}'.format(size_y),
            '{}'.format(size_z_final),
            '-z',
            '{}'.format(z0),
            '{}'.format(z1),
            '--full-gromos87',
            ]

    with TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'sur.gro')

        with open(path, 'w') as fp:
            try:
                subprocess.run(args, stdout=fp)
            except FileNotFoundError as exc:
                print("error: could not find required executable '{}'".format(binary))
                exit(1)

        return read_gromos87(path)


def get_final_conf(conf_stacked_phases, conf_surfactants):
    """Add the surfactant atoms to the fluid stack."""

    try:
        atoms = conf_stacked_phases.atoms + conf_surfactants.atoms
    except:
        atoms = conf_stacked_phases.atoms

    return Gromos87(
            title=conf_stacked_phases.title,
            box_size=conf_stacked_phases.box_size,
            atoms=atoms,
            )


def get_posres(conf, mol_len=1., invert_first=True):
    """Generate position restraint coordinates. 
    
    To be honest I don't remember exactly what this function does.
    It's very weird, like this program overall. 
    
    Great coding, past Petter.
    
    """

    def center_z(atoms, invert):
        z = np.mean([atom.position.z for atom in atoms])

        num_atoms = len(atoms)
        dz = mol_len / float(num_atoms)

        if invert:
            z0 = z + mol_len / 2.
            z1 = z - mol_len / 2.
        else:
            z0 = z - mol_len / 2.
            z1 = z + mol_len / 2.

        zs = np.linspace(z0, z1, num_atoms)

        new_atoms = []
        for (atom, z) in zip(atoms, zs):
            position = Vec3(x=atom.position.x, y=atom.position.y, z=z)

            new_atoms.append(Atom(
                name=atom.name,
                residue=atom.residue,
                position=position,
                velocity=atom.velocity,
                ))

        return new_atoms

    if conf == None:
        return None

    num_mols = len(conf.atoms) // sur_residue_length

    atoms = []

    if invert_first:
        invert = True
    else:
        invert = False

    for i in range(num_mols):
        if i == (num_mols // 2):
            invert = not invert

        i0 = sur_residue_length * i
        i1 = (i + 1) * sur_residue_length

        atoms += center_z(conf.atoms[i0:i1], invert)

    return Gromos87(
            title=conf.title,
            box_size=conf.box_size,
            atoms=atoms
            )


def print_topol(fp, 
                conf_one, conf_two, conf_sur, 
                name_one, name_two, name_sur, 
                residue_length):
    """Write topology [ molecules ] directive."""

    def get_num_mols(conf, num_atoms):
        return len(conf.atoms) // num_atoms

    fp.write("{} {}\n".format(name_one, get_num_mols(conf_one, residue_length)))
    fp.write("{} {}\n".format(name_two, get_num_mols(conf_two, residue_length)))

    if conf_sur:
        num_surfactants = get_num_mols(conf_sur, sur_residue_length)

        if num_surfactants > 0:
            fp.write("{} {}\n".format(name_sur, num_surfactants))

    return 


if __name__ == '__main__':
    parser = ArgumentParser('create_system',
            description='Create a two-phase fluid system with surfactants.')

    parser.add_argument('box_size', 
            nargs=3, metavar='SIZE', type=float, 
            help="Size of each fluid phase along x, y and z.")

    parser.add_argument('-d', '--surfactant-density',
            default=0.215, type=float, metavar='VALUE',
            help='area number density of surfactants in each interface (default: %(default)s)')
    parser.add_argument('-s', '--separation', 
            default=5., metavar='DZ', type=float,
            help='separation between phase one and two along z (default: %(default)s)')

    parser_output = parser.add_argument_group('output options',
            description='Options for writing output files.')
    parser_output.add_argument('-o', '--output', 
            type=str, default='conf_two_phase.gro', metavar='PATH', 
            help='write configuration to file at this path (default: %(default)s)')
    parser_output.add_argument('-r', '--posres-output', 
            type=str, default=None, metavar='PATH', 
            help='optionally write position restraints to file at this path')

    parser_input = parser.add_argument_group('input configuration options')
    parser_input.add_argument('--phase-one', 
            type=str, default=None, metavar='PATH', 
            help='path to phase one base configuration file')
    parser_input.add_argument('--phase-two', 
            type=str, default=None, metavar='PATH', 
            help='path to phase two base configuration file')
    parser_input.add_argument('--residue-length', 
            default=2, type=int, metavar='N',
            help='number of atoms per fluid molecule (default: %(default)s)')

    parser_topol = parser.add_argument_group('topology options',
            description="""
                Options which affect the written output in the topology 
                format (Gromacs .top files). This output lists the created 
                molecule groups with their name and number of molecules.
                
                """)
    parser_topol.add_argument('--phase-one-topolname', 
            type=str, default='lj-chain-2',
            metavar='NAME', help='name for phase one molecules (default: %(default)s)')
    parser_topol.add_argument('--phase-two-topolname', type=str, default='lj-chain-2-B',
            metavar='NAME', help='name for phase two molecules (default: %(default)s)')
    parser_topol.add_argument('--surfactant-topolname', type=str, default='lj-surfactant-bond',
            metavar='NAME', help='name for phase two molecules (default: %(default)s)')

    args = parser.parse_args()

    conf_one = read_conf(args.phase_one, 'LJ1.gro')
    conf_two = read_conf(args.phase_two, 'LJ2.gro')

    size_x, size_y, size_z = args.box_size

    conf_one_final = create_conf_with_size(
            conf_one, size_x, size_y, size_z, args.residue_length)
    conf_two_final = create_conf_with_size(
            conf_two, size_x, size_y, size_z, args.residue_length)

    conf_stacked_phases = create_stack(conf_one_final, conf_two_final, args.separation)

    if args.surfactant_density > 0.:
        conf_surfactants = get_surfactant_conf(
                conf_one_final, conf_two_final, args.separation, conf_stacked_phases, args.surfactant_density)
    else:
        conf_surfactants = None

    conf_final = get_final_conf(conf_stacked_phases, conf_surfactants)
    write_gromos87(args.output, conf_final)

    if args.posres_output:
        conf_surfactants_posres = get_posres(conf_surfactants)
        conf_with_posres = get_final_conf(conf_stacked_phases, conf_surfactants_posres)
        write_gromos87(args.posres_output, conf_with_posres)

    print_topol(stdout, 
                conf_one_final, conf_two_final, conf_surfactants,
                args.phase_one_topolname, args.phase_two_topolname, 
                args.surfactant_topolname, args.residue_length)
