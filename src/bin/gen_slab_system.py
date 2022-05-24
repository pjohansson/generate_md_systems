#!/usr/bin/env python3

import datetime
import numpy as np
import os
import subprocess

from argparse import ArgumentParser
from dataclasses import dataclass
from sys import exit, stderr
from tempfile import TemporaryDirectory
from collections.abc import Iterable, Sequence

from generate_md_systems import single_phase
from generate_md_systems.gmx_conf_utils import read_gromos87, write_gromos87
from .gen_two_phase_shear_system import *
from generate_md_systems.gmx_conf_utils import Gromos87, create_gromos87_conf_with_size


@dataclass
class Slab:
    thickness: float
    z0: float


def get_slab_definitions(
    thicknesses: Sequence[float],
    separations: Sequence[float],
    z0: float = 0.,
) -> list[Slab]:
    if len(thicknesses) > 1 and len(separations) == 1:
        separations = separations * len(thicknesses)
    elif len(separations) < len(thicknesses) - 1:
        raise ValueError(
            f"number of `--separation` values must be 1 "
            f"or equal to number of z values - 1 (was {len(separations)})")

    slabs = []

    for i, thickness in enumerate(thicknesses):
        slabs.append(Slab(thickness=thickness, z0=z0))

        try:
            z0 += thickness + separations[i]
        except IndexError:
            break

    return slabs


def get_slabs_total_height(slabs: Sequence[Slab]) -> float:
    if slabs == []:
        return 0.
    else:
        return slabs[-1].z0 + slabs[-1].thickness - slabs[0].z0


def adjust_slabs_z0(slabs: Iterable[Slab], zadd: float) -> list[Slab]:
    for slab in slabs:
        slab.z0 += zadd

    return slabs


def calc_z_shifts(slabs: Iterable[Slab],
                  fcc_thickness: float,
                  fcc_margin: float,
                  box_margin: float,
                  ) -> dict[str, float]:
    bottom_z0 = box_margin
    slabs_z0 = bottom_z0 + fcc_thickness + fcc_margin
    top_z0 = slabs_z0 + get_slabs_total_height(slabs) + fcc_margin
    box_height = top_z0 + fcc_thickness + box_margin

    return {
        'sub_bottom': bottom_z0,
        'sub_top': top_z0,
        'slabs': slabs_z0,
        'box_height': box_height,
    }


def read_conf(path: str | None,
              default_path: str,
              default_dir: str = 'include',
              ) -> Gromos87:
    """Read a configuration either from a given path or the database.

    If the input `path` is `None` the `default_path` is used to access
    the files included with the module `single_phase`.

    """

    def verify_conf_is_monotype(conf):
        """Asserts that all residues in the configuration are identical.

        Raises an `AssertionError` otherwise.

        """

        head_residue = conf.atoms[0].residue

        for atom in conf.atoms[1:]:
            assert atom.residue == head_residue

    if not path:
        for module_dir in single_phase.__path__:
            path = os.path.join(module_dir, default_dir, default_path)

            if os.path.exists(path):
                break

        try:
            assert os.path.exists(path)
        except AssertionError as exc:
            stderr.write(
                f"error: could not read included data file '{default_path}', "
                "was it installed correctly with the module? "
                f"(full path: '{path}')\n")
            exit(1)

    conf = read_gromos87(path)

    try:
        verify_conf_is_monotype(conf)
    except AssertionError:
        stderr.write(
            f"error: configuration '{path}' includes more than "
            "one residue type\n")
        exit(1)

    return conf


if __name__ == '__main__':
    parser = ArgumentParser(
        description="""
            Generate a slab system for Poiseuille flow simulations.
            """
    )

    parser.add_argument('x',
                        type=float,
                        help="size of system along x")
    parser.add_argument('y',
                        type=float,
                        help="size of system along y")
    parser.add_argument('z',
                        type=float, nargs='*',
                        help="thickness of each liquid slab (along z)")

    parser_phase = parser.add_argument_group('liquid phase options')
    parser_phase.add_argument('-s', '--separation',
                              nargs='+', type=float, default=[1.], metavar='DZ',
                              help="separation distance between liquid slabs (default: %(default)s)")
    parser_phase.add_argument('--liquid-conf',
                              type=str, default=None, metavar='PATH',
                              help="optional path to liquid slab base configuration")
    parser_phase.add_argument('--liquid-topolname',
                              default='lj-fluid', type=str, metavar='NAME',
                              help="name of liquid molecule in topology file (default: %(default)s)")
    parser_phase.add_argument('--liquid-residue-length',
                              default=1, type=str, dest='residue_length', metavar='N',
                              help="length of residue/molecule in liquid phase (default: %(default)s)")

    parser_fcc = parser.add_argument_group('substrate options')
    parser_fcc.add_argument('--fcc-spacing',
                            type=float, default=1., metavar='DX',
                            help="spacing factor of substrate FCC lattice (default: %(default)s)")
    parser_fcc.add_argument('--fcc-nz',
                            type=int, default=5, metavar='N',
                            help="number of layers for the FCC substrate (default: %(default)s)")
    parser_fcc.add_argument('--fcc-margin',
                            type=float, default=1., metavar='DZ',
                            help="margin between substrate and liquid phases")
    parser_fcc.add_argument('--box-margin',
                            type=float, default=5., metavar='DZ',
                            help="margin between box edges and substrates (along z)")

    parser_output = parser.add_argument_group('output options')
    parser_output.add_argument('-o', '--output',
                               type=str, default='conf_final.gro', metavar='PATH',
                               help="output path for final system configuration (default: %(default)s)")
    parser_output.add_argument('-t', '--topology',
                               type=str, metavar='PATH', default=None,
                               help="optionally write topology `[ molecules ]` directive to this path")
    parser_output.add_argument('--title',
                               type=str, default="Lennard-Jones slab system",
                               help="set title for output configuration")

    args = parser.parse_args()

    slabs = get_slab_definitions(args.z, args.separation, 0.)

    nx, ny = calc_fcc_num_sites(args.x, args.y, args.fcc_spacing)
    spacing_y = (np.sqrt(3.) / 2.) * args.fcc_spacing

    fcc_box_x = nx * args.fcc_spacing
    fcc_box_y = ny * spacing_y
    fcc_box_z = calc_fcc_box_z(args.fcc_nz, args.fcc_spacing)

    zshifts = calc_z_shifts(slabs, fcc_box_z, args.fcc_margin, args.box_margin)
    slabs = adjust_slabs_z0(slabs, zshifts['slabs'])

    conf_orig = read_conf(args.liquid_conf, 'LJ1.gro')

    # Process all files in a temporary directory which will be
    # automatically cleaned up afterwards. This is nice.
    #
    # Note: Consider letting the user specify a directory, which
    # should keep the in-between files.
    with TemporaryDirectory() as tmpdir:
        path_substrate = os.path.join(tmpdir, 'sub_notrans.gro')

        path_substrate_bottom = os.path.join(tmpdir, 'sub_bottom.gro')
        path_substrate_top = os.path.join(tmpdir, 'sub_top.gro')
        paths_slabs = [
            os.path.join(tmpdir, f'phases_{i}.gro')
            for i in range(len(slabs))
        ]

        path_topol_substrate = os.path.join(tmpdir, 'topol-sub.top')
        paths_slab_topols = [
            os.path.join(tmpdir, f'phases_top_{i}.top')
            for i in range(len(slabs))
        ]

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

        topdata_slabs = []

        for slab, fnslab in zip(slabs, paths_slabs):
            conf = create_gromos87_conf_with_size(
                conf_orig, fcc_box_x, fcc_box_y, slab.thickness,
                residue_length=args.residue_length,
                translate=[0., 0., slab.z0])
            write_gromos87(fnslab, conf)

            topdata_slabs.append([
                args.liquid_topolname,
                len(conf.atoms) // args.residue_length
            ])

        combine_args = [
            'cat_gro_files.py',
            path_substrate_bottom,
            path_substrate_top,
            *paths_slabs,
            '--output', args.output,
            '--box_size', '-1', '-1', str(zshifts['box_height']),
            '--title', get_output_title(args.title),
            '--quiet',
        ]

        supress_stdout = {'stdout': subprocess.DEVNULL}

        try:
            subprocess.run(mksubstrate_args)
            subprocess.run(editconf_sub_bottom_args, **supress_stdout)
            subprocess.run(editconf_sub_top_args, **supress_stdout)

            subprocess.run(combine_args, **supress_stdout)

            if args.topology:
                combine_topol_args = [
                    'cat',
                    path_topol_substrate,
                    path_topol_substrate,
                ]

                with open(args.topology, 'w') as fp:
                    subprocess.run(combine_topol_args, stdout=fp)

                    for name, num in topdata_slabs:
                        fp.write(f"{name} {num}\n")

        except FileNotFoundError as exc:
            print("error: could not find required executable ({})".format(exc))
            exit(1)
