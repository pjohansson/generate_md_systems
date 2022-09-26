#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

from generate_md_systems.gmx_conf_utils import (
    Atom,
    Gromos87,
    create_gromos87_conf_with_size,
    read_gromos87,
    write_gromos87,
    write_gromos87_stdout,
)


def cut_conf(
    conf: Gromos87,
    origin: tuple[float, float, float],
    end: tuple[float, float, float],
) -> Gromos87:
    def is_inside(position: tuple[float, float, float]) -> bool:
        x, y, z = position

        return (
            (x >= x0) and (x <= x1)
            and (y >= y0) and (y <= y1)
            and (z >= z0) and (z <= z1)
        )

    x0, y0, z0 = origin
    x1, y1, z1 = end

    atoms = []

    for atom in conf.atoms:
        if not is_inside(atom.position):
            atoms.append(atom)

    return Gromos87(
        title=conf.title,
        atoms=atoms,
        box_size=conf.box_size,
    )


def calc_extent(
    box_size: tuple[float, float, float],
    extent: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        box if ext <= 0. else ext for box, ext in zip(box_size, extent)
    )


def calc_cut_with_padding(
    origin: tuple[float, float, float],
    end: tuple[float, float, float],
    padding: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    origin_cut = tuple(v - padding for v in origin)
    end_cut = tuple(v + padding for v in end)

    return origin_cut, end_cut


def add_confs(*configurations: Gromos87) -> Gromos87:
    title = configurations[0].title
    box_size = configurations[0].box_size
    atoms = []

    for conf in configurations:
        atoms += conf.atoms

    return Gromos87(
        title=title,
        atoms=atoms,
        box_size=box_size,
    )

    

if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument(
        'conf_target',
        metavar='TARGET',
        help="target configuration to insert into",
    )
    parser.add_argument(
        'conf_insert',
        metavar='INSERT',
        help="inserted configuration",
    )
    parser.add_argument(
        '-o', '--output',
        default=None, metavar='PATH', type=str,
        help="write final configuration to disk instead of stdout",
    )
    parser.add_argument(
        '-0', '--origin',
        nargs=3, default=(0., 0., 0.), metavar=('X0', 'Y0', 'Z0'), type=float,
        help="corner to insert from (default: %(default)s)",
    )
    parser.add_argument(
        '-1', '--extent',
        nargs=3, default=(-1., -1., -1.), metavar=('DX', 'DY', 'DZ'), type=float,
        help="size of inserted conf, =0 for original size (default: %(default)s)",
    )
    parser.add_argument(
        '-p', '--padding',
        type=float, default=0., metavar='PAD',
        help="padding around insertet phase (default: %(default)s)",
    )

    args = parser.parse_args()

    conf_target = read_gromos87(args.conf_target)
    conf_insert = read_gromos87(args.conf_insert)

    origin = args.origin
    extent = calc_extent(conf_insert.box_size, args.extent)
    end = tuple(o + e for o, e in zip(origin, extent))

    origin_cut, end_cut = calc_cut_with_padding(origin, end, args.padding)
    conf_target_cut = cut_conf(conf_target, origin_cut, end_cut)

    conf_insert_fill = create_gromos87_conf_with_size(
        conf_insert, 
        *extent, 
        residue_length=1, 
        translate=origin,
    )

    conf_final = add_confs(conf_target_cut, conf_insert_fill)

    if args.output == None:
        write_gromos87_stdout(conf_final)
    else:
        write_gromos87(args.output, conf_final)
