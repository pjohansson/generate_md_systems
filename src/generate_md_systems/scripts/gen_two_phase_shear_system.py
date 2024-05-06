#!/usr/bin/env python3

import datetime
import numpy as np
import os
import subprocess

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory


def calc_fcc_num_sites(x, y, spacing):
    """Calculate the number of sites along x and y."""

    spacing_y = (np.sqrt(3.) / 2.) * spacing

    nx = int(np.ceil(x / spacing))
    ny = int(np.ceil(y / spacing_y))

    return nx, ny


def calc_fcc_box_z(nz, spacing):
    """Calculate the height of the FCC lattice."""

    spacing_z = (np.sqrt(6.) / 3.) * spacing

    return nz * spacing_z


def calc_z_shifts(zphase, zsubstrate, margin):
    """Calculate shifts from z = 0 for each phase."""

    shifts = {}

    shifts['sub_bottom'] = margin

    sub_bottom_upper_margin = shifts['sub_bottom'] + zsubstrate
    shifts['phase'] = sub_bottom_upper_margin + margin

    phase_upper_margin = shifts['phase'] + zphase
    shifts['sub_top'] = phase_upper_margin + margin

    sub_top_upper_margin = shifts['sub_top'] + zsubstrate
    shifts['box_z'] = sub_top_upper_margin + margin

    return shifts


def get_editconf_translate_args(from_path, to_path, shift):
    """Create cmd command to translate a configuration."""

    return [
        'gmx', 'editconf',
        '-f', from_path,
        '-o', to_path,
        '-translate', '0',  '0', str(shift),
        '-quiet',
    ]


def get_output_title(title,
                     default_title="Two-phase system with substrates"):
    """Get set title for output configuration or generate a default."""

    if title:
        return title
    else:
        now = datetime.datetime.now()
        return default_title + ', generated at {}'.format(now)


def main():
    parser = ArgumentParser(
        description="""
            Generate two-phase fluid systems with surfactants at
            interfaces and substrates on top and bottom."""
    )

    parser.add_argument('x',
        type=float,
        help="size of each liquid phase along x")
    parser.add_argument('y',
        type=float,
        help="size of each liquid phase along y")
    parser.add_argument('z',
        type=float,
        help="size of each liquid phase along z")

    # Disable this option: at least for now we only create
    # shear systems directed along the x axis. Otherwise
    # we need to decide along which axis the substrates are
    # defined and rotate them properly to match.
    # parser.add_argument('-a', '--axis',
    #         choices=['x', 'y', 'z'], default='x', type=str.lower,
    #         help="axis along which liquid phases are set (default: %(default)s)")

    parser_phase = parser.add_argument_group('liquid phase options')
    parser.add_argument('--phase2-size',
            type=float, default=None, metavar=('X'),
            help="optionally set size of second phase separately to the first")
    parser_phase.add_argument('-s', '--separation',
        type=float, default=5., metavar='DIST',
        help="separation distance between liquid phases (default: %(default)s)")
    parser_phase.add_argument('-d', '--surfactant-density',
        type=float, default=0.215, metavar='RHO',
        help="areal atom number density of surfactants at each interface (default: %(default)s)")

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

    if args.phase2_size:
        total_phase_size = args.x + args.phase2_size
    else:
        total_phase_size = 2. * args.x

    box_x = total_phase_size + 2. * args.separation

    nx, ny = calc_fcc_num_sites(box_x, args.y, args.fcc_spacing)
    fcc_box_z = calc_fcc_box_z(args.fcc_nz, args.fcc_spacing)

    spacing_y = (np.sqrt(3.) / 2.) * args.fcc_spacing
    fcc_box_y = ny * spacing_y

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

        mkphase_args = [
            'gen_two_phase_system.py',
            str(args.x), str(fcc_box_y), str(args.z),
            '--output', path_phases,
            '--surfactant-density', str(args.surfactant_density),
            '--separation', str(args.separation),
            '--axis', 'x',
            '--topology', path_topol_phases,
        ]

        if args.phase2_size:
            mkphase_args += [
                '--box_size_2',
                str(args.phase2_size), str(fcc_box_y), str(args.z),
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

        supress_stdout = { 'stdout': subprocess.DEVNULL }

        try:
            subprocess.run(mksubstrate_args)
            subprocess.run(mkphase_args)

            subprocess.run(editconf_sub_bottom_args, **supress_stdout)
            subprocess.run(editconf_sub_top_args, **supress_stdout)
            subprocess.run(editconf_phase_args, **supress_stdout)

            subprocess.run(combine_args, **supress_stdout)

            if args.topology:
                combine_topol_args = [
                    'cat',
                    path_topol_substrate,
                    path_topol_substrate,
                    path_topol_phases,
                ]

                with open(args.topology, 'w') as fp:
                    subprocess.run(combine_topol_args, stdout=fp)

        except FileNotFoundError as exc:
            print("error: could not find required executable ({})".format(exc))
            exit(1)
