#!/usr/bin/env python3

import datetime
import numpy as np
import os
import subprocess

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory

from .gen_two_phase_system import read_conf
from .gen_two_phase_shear_system import *
from generate_md_systems.gmx_conf_utils import write_gromos87, create_gromos87_conf_with_size


def write_topol(fp, conf, name, residue_length):
    """Write topology `[ molecules ]` directive."""


    return


if __name__ == '__main__':
    parser = ArgumentParser(
        description="""
            Generate single-phase couette flow fluid systems with 
            substrates on top and bottom.
            """
    )

    parser.add_argument('x', 
        type=float,
        help="size of liquid phase along x")
    parser.add_argument('y', 
        type=float,
        help="size of liquid phase along y")
    parser.add_argument('z', 
        type=float,
        help="size of liquid phase along z")

    parser_phase = parser.add_argument_group('liquid phase options')
    parser_phase.add_argument('--phase-path', 
            default=None, type=str, metavar='PATH',
            help="optionally set explicit path to phase configuration")
    parser_phase.add_argument('--phase-default', 
            default='LJ1.gro', type=str, metavar='FILE',
            help="library configuration to use (default: %(default)s)")
    parser_phase.add_argument('--phase-residue-length', 
            type=int, default=2,
            help="number of atoms per liquid molecule (default: %(default)s)")
    parser_phase.add_argument('--phase-topolname', 
            type=str, default='lj-chain-2',
            help="number of atoms per liquid molecule (default: %(default)s)")

    parser_fcc = parser.add_argument_group('substrate options')
    parser_fcc.add_argument('--fcc-spacing', 
        type=float, default=1., metavar='DX', 
        help="spacing factor of substrate FCC lattice (default: %(default)s)")
    parser_fcc.add_argument('--fcc-nz', 
        type=int, default=5, metavar='N', 
        help="number of layers for the FCC substrate (default: %(default)s)")
    parser_fcc.add_argument('--fcc-margin', 
        type=float, default=4., metavar='DZ', 
        help="margin between substrate and liquid phases")

    parser_output = parser.add_argument_group('output options')
    parser_output.add_argument('-o', '--output', 
        type=str, default='conf_final.gro', metavar='PATH',
        help="output path for final system configuration (default: %(default)s)")
    parser_output.add_argument('-t', '--topology', 
        type=str, metavar='PATH', default=None, 
        help="optionally write topology `[ molecules ]` directive to this path")
    parser_output.add_argument('--title', 
        type=str, default=None, 
        help="set title for output configuration")

    args = parser.parse_args()

    nx, ny = calc_fcc_num_sites(args.x, args.y, args.fcc_spacing)
    spacing_y = (np.sqrt(3.) / 2.) * args.fcc_spacing 

    fcc_box_x = nx * args.fcc_spacing
    fcc_box_y = ny * spacing_y
    fcc_box_z = calc_fcc_box_z(args.fcc_nz, args.fcc_spacing)

    zshifts = calc_z_shifts(args.z, fcc_box_z, args.fcc_margin)

    # Process all files in a temporary directory which will be 
    # automatically cleaned up afterwards. This is nice.
    #
    # Note: Consider letting the user specify a directory, which 
    # should keep the in-between files.
    with TemporaryDirectory() as tmpdir:
        path_substrate = os.path.join(tmpdir, 'sub_notrans.gro')
        path_phases = os.path.join(tmpdir, 'phases_notrans.gro')

        path_substrate_bottom = os.path.join(tmpdir, 'sub_bottom.gro')
        path_substrate_top = os.path.join(tmpdir, 'sub_top.gro')
        path_phases_trans = os.path.join(tmpdir, 'phases_trans.gro')

        path_topol_substrate = os.path.join(tmpdir, 'topol-sub.top')
        path_topol_phases = os.path.join(tmpdir, 'topol-phases.top')

        # Substrate configuration and all translations are done using
        # external commands. Prepare the calls here, then call below
        # in the designated section.
        #
        # Translations are done using `gmx editconf`, which has to be
        # available. 
        mksubstrate_args = [
            'make_LJ_substrate.py',
            str(nx), str(ny), str(args.fcc_nz),
            '--output', path_substrate,
            '--spacing', str(args.fcc_spacing),
            '--topology', path_topol_substrate,
        ]

        editconf_sub_bottom_args = get_editconf_translate_args(
            path_substrate, 
            path_substrate_bottom, 
            zshifts['sub_bottom'])

        editconf_sub_top_args = get_editconf_translate_args(
            path_substrate, 
            path_substrate_top, 
            zshifts['sub_top'])

        editconf_phase_args = get_editconf_translate_args(
            path_phases, 
            path_phases_trans, 
            zshifts['phase'])

        combine_args = [ 
            'cat_gro_files.py',
            path_substrate_bottom, 
            path_substrate_top, 
            path_phases_trans,
            '--output', args.output,
            '--box_size', '-1', '-1', str(zshifts['box_z']),
            '--title', get_output_title(args.title),
            '--quiet',
        ]

        suppress_stdout = { 'stdout': subprocess.DEVNULL }

        try:
            # Creating and writing the single liquid phase is done internally,
            # without calling an external command as prepared above.
            conf_orig = read_conf(args.phase_path, args.phase_default)
            conf_final = create_gromos87_conf_with_size(
                conf_orig, fcc_box_x, fcc_box_y, args.z, args.phase_residue_length)
            write_gromos87(path_phases, conf_final)

            subprocess.run(mksubstrate_args)

            subprocess.run(editconf_sub_bottom_args, **suppress_stdout)
            subprocess.run(editconf_sub_top_args, **suppress_stdout)
            subprocess.run(editconf_phase_args, **suppress_stdout)

            subprocess.run(combine_args, **suppress_stdout)

            if args.topology:
                combine_topol_args = [
                    'cat',
                    path_topol_substrate, 
                    path_topol_substrate, 
                ]

                with open(args.topology, 'w') as fp:
                    subprocess.run(combine_topol_args, stdout=fp)

                    num_mols = len(conf_final.atoms) // args.phase_residue_length
                    fp.write(f"{args.phase_topolname} {num_mols}\n")

        except FileNotFoundError as exc:
            print("error: could not find required executable ({})".format(exc))
            exit(1)
