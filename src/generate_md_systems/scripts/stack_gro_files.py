#!/usr/bin/env python3

from argparse import ArgumentParser
from generate_md_systems.gmx_conf_utils import *


def get_box_size(
    box_size_per_conf: list[tuple[float, float, float]],
    axis: str,
    separation: float,
    add_edges: float, 
) -> tuple[float, float, float]:
    max_x = 0.
    max_y = 0.
    max_z = 0.

    num_separations = len(box_size_per_conf) - 1
    height = 2. * add_edges + num_separations * separation

    for box_x, box_y, box_z in box_size_per_conf:
        if box_x > max_x: 
            max_x = box_x 
        if box_y > max_y: 
            max_y = box_y 
        if box_z > max_z: 
            max_z = box_z 

        if axis == 'x':
            height += box_x
        elif axis == 'y':
            height += box_y
        else:
            height += box_z

    if axis == 'x':
        max_x = height
    elif axis == 'y':
        max_y = height
    else:
        max_z = height
    
    return Vec3(x=max_x, y=max_y, z=max_z)


def main():
    parser = ArgumentParser()

    parser.add_argument(
        'paths', 
        nargs='+', metavar='PATH',
        help=".gro files to stack",
    )
    parser.add_argument(
        '-a', '--axis',
        type=str.lower, choices=['x', 'y', 'z'], default='z',
        help="axis along which to stack the files (default: %(default)s)",
    )
    parser.add_argument(
        '-s', '--separation',
        type=float, default=0., metavar='VALUE',
        help="separation between stacked phases",
    )
    parser.add_argument(
        '-o', '--output',
        metavar='PATH',
        help="write output to given path instead of to stdout",
    )
    parser.add_argument(
        '-0', '--add_edges',
        type=float, default=0., metavar='DX',
        help="add extra space to bottom and top",
    )
    parser.add_argument(
        '-t', '--title',
        default=None, 
        help="set title of output configuration",
    )

    args = parser.parse_args()

    atoms = []
    box_size_per_conf = []

    pos = 0.
    translate = [0., 0., 0.]
    title = args.title

    if args.axis == 'x':
        translate[0] = args.add_edges 
    elif args.axis == 'y':
        translate[1] = args.add_edges 
    else:
        translate[2] = args.add_edges 

    for path in args.paths:
        conf = translate_gromos87(read_gromos87(path), *translate)

        atoms += conf.atoms
        box_size_per_conf.append(conf.box_size)

        if title == None:
            title = conf.title

        if args.axis == 'x':
            translate[0] += conf.box_size.x + args.separation
        elif args.axis == 'y':
            translate[1] += conf.box_size.y + args.separation
        else:
            translate[2] += conf.box_size.z + args.separation
    
    box_size = get_box_size(box_size_per_conf, args.axis, args.separation, args.add_edges)
    conf_stacked = Gromos87(title, atoms, box_size)

    if args.output == None:
        write_gromos87_stdout(conf_stacked)
    else:
        write_gromos87(args.output, conf_stacked)
