#!/usr/bin/env python3

import numpy as np
import os
import subprocess

from generate_md_systems import two_phase
from generate_md_systems.gmx_conf_utils import *

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory

# Surfactant molecule area density in interface (used in other simulations)
# sur_area_density = (344 / 2) / (40. * 20.)
sur_residue_length = 2


def read_conf(path, default_path, default_dir='include'):
    """Read a configuration either from a given path or the database.
    
    If the input `path` is `None` the `default_path` is used to access
    the files included with the module `two_phase`.

    """

    def verify_conf_is_monotype(conf):
        """Asserts that all residues in the configuration are identical.
        
        Raises an `AssertionError` otherwise.
        
        """

        head_residue = conf.atoms[0].residue

        for atom in conf.atoms[1:]:
            assert atom.residue == head_residue

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


def translate_atoms(atoms, shift, axis):
    """Translate a set of `Atom`s by `shift` along `axis`."""

    def move_atom(atom, shift, axis):
        x, y, z = atom.position

        if axis == 'x':
            x += shift
        elif axis == 'y':
            y += shift
        else:
            z += shift

        position = Vec3(x, y, z)

        return Atom(
                name=atom.name,
                residue=atom.residue,
                position=position,
                velocity=atom.velocity
                )

    return [move_atom(atom, shift, axis) for atom in atoms]


def create_stack(conf_one, conf_two, separation, axis):
    """Translate and stack the two phases with the given separation."""

    if axis == 'x':
        d = 0
    elif axis == 'y':
        d = 1
    else:
        d = 2
    
    size1 = conf_one.box_size[d]
    size2 = conf_two.box_size[d]

    shift1 = size2 / 2. + separation
    shift2 = -size2 / 2.

    conf_one_atoms = translate_atoms(conf_one.atoms, shift1, axis)
    conf_two_atoms = translate_atoms(conf_two.atoms, shift2, axis)

    atoms = conf_one_atoms + conf_two_atoms

    box_x, box_y, box_z = conf_one.box_size
    new_size = size1 + size2 + 2. * separation

    if axis == 'x':
        box_x = new_size
    elif axis == 'y':
        box_y = new_size
    else:
        box_z = new_size

    return Gromos87(
            title=conf_one.title,
            box_size=(box_x, box_y, box_z),
            atoms=atoms,
            )


def get_surfactant_conf(conf_one, conf_two, 
                        separation, conf_stacked, 
                        sur_area_density, axis):
    """Generate surfactants at both interfaces."""

    def rotate_conf_to_axis(conf, axis):
        def rotate_atoms(atoms, matrix):
            def rotate_single_atom(atom):
                def multiply_row(row):
                    return row[0] * x0 + row[1] * y0 + row[2] * z0

                x0, y0, z0 = atom.position

                x1 = multiply_row(matrix[0])
                y1 = multiply_row(matrix[1])
                z1 = multiply_row(matrix[2])

                position = Vec3(x1, y1, z1)

                return Atom(
                    position=position,
                    velocity=atom.velocity,
                    name=atom.name,
                    residue=atom.residue,
                )
            
            return [rotate_single_atom(atom) for atom in atoms]

        box_x, box_y, box_z = conf.box_size

        if axis == 'x':
            # Rotate -90 deg. around the Y axis
            matrix = [
                [0, 0, 1],
                [0, 1, 0],
                [-1, 0, 0],
            ]

            atoms = rotate_atoms(conf.atoms, matrix)
            atoms = translate_atoms(atoms, box_x, 'z')
            box_size = box_z, box_y, box_x

        elif axis == 'y':
            # Rotate -90 deg. around the X axis
            matrix = [
                [1, 0, 0],
                [0, 0, 1],
                [0, -1, 0],
            ]

            atoms = rotate_atoms(conf.atoms, matrix)
            atoms = translate_atoms(atoms, box_y, 'z')
            box_size = box_z, box_x, box_y

        else:
            atoms = conf.atoms
            box_size = conf.box_size

        return Gromos87(
            title=conf.title,
            atoms=atoms,
            box_size=box_size,
        )

    # `gen_surfactant` only generates interfaces normal to z. 
    #
    # To get around this we make a transform into our box:
    #   `sz` here contains our box size along the target `axis`,
    #   `sx`, `sy` the box sizes transverse to it.
    #
    # These will be used for the surfactant generation, after 
    # which we rotate the box to our actual axis
    #
    # Note: We here first to the rotation and then translation 
    # separately. This could be done in a single which would 
    # likely save some cycles. But since there generally are 
    # not that many surfactant molecules we leave it like this 
    # for now.

    box_full_x, box_full_y, box_full_z = conf_stacked.box_size

    if axis == 'x': 
        d = 0
        sx = box_full_z
        sy = box_full_y 
        sz = box_full_x
    elif axis == 'y':
        d = 1
        sx = box_full_x 
        sy = box_full_z 
        sz = box_full_y
    else:
        d = 2
        sx = box_full_x 
        sy = box_full_y
        sz = box_full_z

    sz_one = conf_one.box_size[d]
    sz_two = conf_two.box_size[d]

    z0 = (sz_two + separation) / 2.
    z1 = z0 + sz_one + separation

    area = sx * sy
    num_mols = int(sur_area_density * area)

    binary = 'gen_surfactant'
    args = [
            binary,
            '{}'.format(num_mols),
            '-b',
            '{}'.format(sx),
            '{}'.format(sy),
            '{}'.format(sz),
            '-z',
            '{}'.format(z0),
            '{}'.format(z1),
            '--full-gromos87',
            ]

    # `gen_surfactant` can write to standard output, so we 
    # open a temporary file to pipe into. We then read the 
    # configuration of surfactants from that file.
    #
    # Very convoluted.
    with TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'sur.gro')

        with open(path, 'w') as fp:
            try:
                subprocess.run(args, stdout=fp)
            except FileNotFoundError as exc:
                print("error: could not find required executable '{}'".format(binary))
                exit(1)

        conf_surfactants = read_gromos87(path)
    
    return rotate_conf_to_axis(conf_surfactants, axis)


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
    """Write topology `[ molecules ]` directive."""

    def get_num_mols(conf, num_atoms):
        return len(conf.atoms) // num_atoms

    fp.write("{} {}\n".format(name_one, get_num_mols(conf_one, residue_length)))
    fp.write("{} {}\n".format(name_two, get_num_mols(conf_two, residue_length)))

    if conf_sur:
        num_surfactants = get_num_mols(conf_sur, sur_residue_length)

        if num_surfactants > 0:
            fp.write("{} {}\n".format(name_sur, num_surfactants))

    return 


def check_phase_size_axis_alignment(box_size, box_size_2, axis):
    def axis_to_string(ax):
        return ['x', 'y', 'z'][ax]

    def check_along_axes(*axes):
        all_good = True

        for ax in axes:
            box = box_size[ax]
            box2 = box_size_2[ax]

            if not np.isclose(box, box2):
                stderr.write("WARNING: phase sizes are different along axis "
                      "'{}' (1: {}, 2: {})\n".format(
                    axis_to_string(ax),
                    box, 
                    box2,
                ))

                all_good = False
        
        return all_good

    if axis == 'x':
        return check_along_axes(1, 2)
    elif axis == 'y':
        return check_along_axes(0, 2)
    else:
        return check_along_axes(0, 1)


if __name__ == '__main__':
    parser = ArgumentParser(
            description='Create a two-phase fluid system with surfactants.')

    parser.add_argument('box_size', 
            nargs=3, metavar='SIZE', type=float, 
            help="Size of each fluid phase along x, y and z.")

    parser.add_argument('-a', '--axis', 
            choices=['x', 'y', 'z'], default='z', 
            type=str.lower, # makes the choices case insensitive
            help="axis along which to create separate phases")
    parser.add_argument('-d', '--surfactant-density',
            default=0.215, type=float, metavar='VALUE',
            help='area number density of surfactants in each interface (default: %(default)s)')
    parser.add_argument('--box_size_2', 
            type=float, nargs=3, default=None, metavar='SIZE',
            help="set size of second phase instead of using `box_size`")
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
    parser_topol.add_argument('-t', '--topology', 
            type=str, metavar='PATH', default=None, 
            help="optionally write topology `[ molecules ]` directive to this path")
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

    if args.box_size_2 != None:
        size_x2, size_y2, size_z2 = args.box_size_2
        check_phase_size_axis_alignment(args.box_size, args.box_size_2, args.axis)
    else:
        size_x2, size_y2, size_z2 = args.box_size


    conf_one_final = create_gromos87_conf_with_size(
            conf_one, size_x, size_y, size_z, residue_length=sur_residue_length)
    conf_two_final = create_gromos87_conf_with_size(
            conf_two, size_x2, size_y2, size_z2, residue_length=sur_residue_length)

    conf_stacked_phases = create_stack(conf_one_final, conf_two_final, 
            args.separation, args.axis)

    if args.surfactant_density > 0.:
        conf_surfactants = get_surfactant_conf(
                conf_one_final, conf_two_final, args.separation, conf_stacked_phases, args.surfactant_density, args.axis)
    else:
        conf_surfactants = None

    conf_final = get_final_conf(conf_stacked_phases, conf_surfactants)
    write_gromos87(args.output, conf_final)

    if args.posres_output:
        conf_surfactants_posres = get_posres(conf_surfactants)
        conf_with_posres = get_final_conf(conf_stacked_phases, conf_surfactants_posres)
        write_gromos87(args.posres_output, conf_with_posres)

    if args.topology:
        with open(args.topology, 'w') as fp:
            print_topol(fp, 
                        conf_one_final, conf_two_final, conf_surfactants,
                        args.phase_one_topolname, args.phase_two_topolname, 
                        args.surfactant_topolname, args.residue_length)
