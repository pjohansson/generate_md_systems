#!/usr/bin/env python3

import datetime

from generate_md_systems.gmx_conf_utils import * 

from argparse import ArgumentParser
from sys import stderr, stdout 


def get_largest_box_size(confs, use_box_size):
    xmax, ymax, zmax = 0., 0., 0.
    use_xmax, use_ymax, use_zmax = use_box_size

    for conf in confs:
        box_x, box_y, box_z = conf.box_size 

        if box_x > xmax: 
            xmax = box_x 
        
        if box_y > ymax: 
            ymax = box_y 

        if box_z > zmax: 
            zmax = box_z
    
    if use_xmax >= 0.: 
        xmax = use_xmax

    if use_ymax >= 0.: 
        ymax = use_ymax

    if use_zmax >= 0.: 
        zmax = use_zmax
        
    return xmax, ymax, zmax


def get_title(confs, use_title):
    """Return the explicit title if present, otherwise the first."""

    if use_title != None:
        return use_title 

    try:
        title = confs[0].title
    # This should never trigger but why not catch it just because
    except IndexError:
        now = datetime.datetime.now()
        title = "Conf generated at {}".format(now)
    
    return title


def main():
    parser = ArgumentParser(
        description="Concatenate Gromos87 configuration (.gro) files."
    )

    parser.add_argument('files', 
        type=str, nargs='+', metavar='PATH',
        help=".gro files to concatenate")
    parser.add_argument('-o', '--output', 
        type=str, default=None, metavar='PATH', 
        help="write final configuration to this file instead of to stdout")
    parser.add_argument('-t', '--title', 
        type=str, default=None, 
        help="set output configuration title")
    parser.add_argument('--box_size', 
        nargs=3, type=float, metavar=('X', 'Y', 'Z'), default=(-1, -1, -1),
        help="set output box size (use <0 to not change)")
    parser.add_argument('-q', '--quiet', 
        action='store_true', 
        help="be less loud and noisy")

    args = parser.parse_args()

    confs = [read_gromos87(path) for path in args.files]

    atoms = []
    for conf in confs:
        atoms += conf.atoms

    title = get_title(confs, args.title)
    box_size = get_largest_box_size(confs, args.box_size)

    conf_final = Gromos87(title, atoms, box_size)

    if not args.quiet:
        stderr.write("Input:\n")
        for fn, conf in zip(args.files, confs):
            stderr.write("  {}: {} atoms\n".format(fn, len(conf.atoms)))

        stderr.write("\n")
        stderr.write("Output:\n")
        stderr.write("  {}\n".format(conf_final.title))
        stderr.write("  {} atoms\n".format(len(conf_final.atoms)))
        stderr.write("  {:g} x {:g} x {:g}\n".format(*box_size))

    if args.output:
        write_gromos87(args.output, conf_final)
    else:
        try:
            write_gromos87_stdout(conf_final)
        except BrokenPipeError:
            pass
