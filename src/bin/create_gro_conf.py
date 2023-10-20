#!/usr/bin/env python3

import datetime
import numpy as np
import os
import subprocess

from generate_md_systems.gmx_conf_utils import *

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from tempfile import TemporaryDirectory


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
    conf_final = create_gromos87_conf_with_size(
        conf, args.x, args.y, args.z, translate=args.translate)
    
    if args.output == None:
        write_gromos87_stdout(conf_final)
    else:
        write_gromos87(args.output, conf_final)
