#!/usr/bin/env python3

import datetime
import numpy as np

from argparse import ArgumentParser
from sys import stderr, stdout


def gen_fcc(nx, ny, nz, spacing):
    """Generate points in an FCC lattice with the given spacing."""

    points = []

    dx = spacing

    dy = (np.sqrt(3.) / 2.) * spacing
    dx_y = 0.5 * spacing

    dz = (np.sqrt(6) / 3.) * spacing
    dx_z = 0.5 * spacing
    dy_z = (1. / 3.) * spacing

    box_size = [dx * nx, dy * ny, dz * nz]

    x0 = dx / 2.
    y0 = dy / 2.
    z0 = dz / 2.

    z = z0

    for iz in range(nz):
        y = y0 + iz * dy_z

        for iy in range(ny):
            x = x0 + iy * dx_y + iz * dx_z

            for ix in range(nx):
                points.append([x, y, z])

                x += dx
            y += dy
        z += dz

    return points, box_size


def apply_pbc(points, box_size):
    """Apply periodic boundary conditions to keep points in the box."""

    def pbc_1d(v, vmax):
        while v < 0.:
            v += vmax

        return v % vmax

    points_pbc = []
    xmax, ymax, zmax = box_size

    for x, y, z in points:
        point = [
            pbc_1d(x, xmax),
            pbc_1d(y, ymax),
            pbc_1d(z, zmax),
        ]

        points_pbc.append(point)

    return points_pbc


def print_conf(fp, points, atom_name, residue_name, box_size, title=None):
    """Write points as a .gro configuration file.

    In case no title is specified one will be generated with a timestamp.

    """

    if title:
        fp.write('{}\n'.format(title))
    else:
        now = datetime.datetime.now()
        fp.write('FCC lattice generated at {}\n'.format(now))

    fp.write('{}\n'.format(len(points)))

    for i, (x, y, z) in enumerate(points):
        fp.write(
            '{num:>5}{atom:<5}{residue:>5}{num:<5}'
            '{x:>8.3f}{y:>8.3f}{z:>8.3f}\n'.format(
                x=x, y=y, z=z,
                num=(i+1) % 100000,
                atom=atom_name[:5],
                residue=residue_name[:5]
                ))

    fp.write('{:9.5f} {:9.5f} {:9.5f}\n'.format(*box_size))

    return


def print_topol(fp, points, topolname):
    """Write topology `[ molecules ]` directive."""

    fp.write("{} {}\n".format(topolname, len(points)))

    return


def main():
    parser = ArgumentParser()

    parser.add_argument(
        'nx', type=int, metavar='NX',
        help="Number of sites along the x axis")
    parser.add_argument(
        'ny', type=int, metavar='NY',
        help="Number of sites along the y axis")
    parser.add_argument(
        'nz', type=int, metavar='NZ', default=5, nargs='?',
        help="Number of sites along the z axis (default: %(default)s)")

    parser.add_argument(
        '-o', '--output', type=str, metavar='PATH', default=None,
        help="Write output to file instead of stdout")
    parser.add_argument(
        '-s', '--spacing', metavar='DX', type=float, default=1.0,
        help="Set spacing between sites (default: %(default)s)")
    parser.add_argument(
        '--apply_pbc', action='store_true',
        help="Apply periodic boundary conditions to coordinates")

    parser.add_argument(
        '-n', '--atom_name', default='CUB', metavar="NAME",
        help="Set atom name for points (default: %(default)s)")
    parser.add_argument(
        '-r', '--residue_name', default='CUB', metavar="NAME",
        help="Set residue name for points (default: %(default)s)")

    parser.add_argument(
        '-t', '--topology', type=str, metavar='PATH', default=None,
        help="Optionally write topology `[ molecules ]` directive to this path")
    parser.add_argument(
        '--topology-name', type=str, metavar='NAME', default='lj-substrate',
        help="Name of substrate in `[ molecules ]` topology directive (default: %(default)s)")

    parser.add_argument(
        '-q', '--quiet', action='store_true',
        help="Be quiet")

    args = parser.parse_args()

    points, box_size = gen_fcc(args.nx, args.ny, args.nz, args.spacing)

    if args.apply_pbc:
        points = apply_pbc(points, box_size)

    if args.output:
        fp = open(args.output, 'w')
    else:
        fp = stdout

    print_conf(fp, points, args.atom_name, args.residue_name, box_size)

    if args.output:
        fp.close()

    if args.topology:
        with open(args.topology, 'w') as fp:
            print_topol(fp, points, args.topology_name)

    if not args.quiet:
        stderr.write(
            "generated {} points in a box of size "
            "{:g} x {:g} x {:g}\n".format(len(points), *box_size))
