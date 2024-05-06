#!/usr/bin/env python3

import numpy as np
import toml

from argparse import ArgumentParser
from collections import namedtuple
from dataclasses import dataclass
from sys import exit, stderr
from typing import Any, Mapping, Sequence

from generate_md_systems.gmx_conf_utils import (
    Atom,
    Gromos87,
    write_gromos87,
    write_gromos87_stdout,
)


Spacing = namedtuple('FccSpacing', ['dx', 'dy', 'dz', 'dx_y', 'dx_z', 'dy_z'])

@dataclass
class AtomSpec:
    name: str
    dx: tuple[float, float, float]

@dataclass
class LatticeSpec:
    title: str
    nx: int
    ny: int
    nz: int
    spacing: float
    residue_name: str
    atoms: list[AtomSpec]
    layers: dict[int, str]
    translate: tuple[float, float, float]

def read_fcc_spec(path: str) -> LatticeSpec:
    def read_value(category: str, key: str, of_type):
        try:
            return of_type(fcc_toml[category][key])
        except KeyError:
            raise Exception(f"missing field `{key} = {of_type.__name__}` in [ {category} ]")

    def check_vec3(value: Sequence[Any]) -> tuple[float, float, float]:
        x, y, z = value
        return (float(x), float(y), float(z))

    def read_site(site_def) -> AtomSpec:
        try:
            name = site_def['name']
        except KeyError:
            raise Exception(f"missing field `name = str` in atoms.sites")

        try:
            dx, dy, dz = check_vec3(site_def['dx'])
        except KeyError:
            raise Exception(f"missing field `dx = [float, float, float]` in atoms.sites (atom name  = {name})")
        except ValueError as exc:
            raise Exception(f"expected `dx = [float, float, float] in atoms.sites (atom name = {name})")

        return AtomSpec(
            name,
            dx=(dx, dy, dz),
        )

    try:
        fcc_toml = toml.load(path)
    except Exception:
        raise Exception("could not read file as .toml, is the path correct?")

    title = read_value('system', 'title', str)
    nx = read_value('system', 'nx', int)
    ny = read_value('system', 'ny', int)
    nz = read_value('system', 'nz', int)
    spacing = read_value('system', 'spacing', float)

    residue_name = read_value('atoms', 'residue', str)

    try:
        atom_name = read_value('atoms', 'name', str)
    except:
        atom_name = residue_name

    atoms = [read_site(v) for v in read_value('atoms', 'sites', list)]

    if atoms == []:
        atoms = [AtomSpec(name=atom_name, dx=(0., 0., 0.))]

    layers_spec = fcc_toml.get('layers', {})

    try:
        layers = {int(layer): spec for layer, spec in layers_spec.items()}
    except ValueError as exc:
        raise Exception(f"invalid layer index in [ layers ] ({exc})")

    try:
        translate = check_vec3(read_value('system', 'translate', list))
    except KeyError:
        raise Exception(f"missing field `translate = [float, float, float]`")
    except ValueError as exc:
        raise Exception(f"expected `translate = [float, float, float])")
    except Exception:
        translate = (0., 0., 0.)

    return LatticeSpec(
        title,
        nx,
        ny,
        nz,
        spacing,
        residue_name,
        atoms,
        layers,
        translate,
    )


def get_layer_inds_from_spec(
    index_layer: int,
    nx: int,
    ny: int,
    fcc_layers: Mapping[int, str],
    chars_full: Sequence[str] = ['+'],
    chars_empty: Sequence[str] = ['-'],
) -> np.ndarray:
    def layer_is_match(layer: str, haystack: Sequence[str]) -> bool:
        return np.all([c in haystack for c in layer])  # type: ignore

    layer_str = fcc_layers.get(index_layer, '-')

    if layer_str == "" or layer_is_match(layer_str, chars_empty):
        return np.zeros((nx, ny), dtype=bool)
    elif layer_is_match(layer_str, chars_full):
        return np.ones((nx, ny), dtype=bool)

    layer_1d = np.array(
        [int(i) for i in layer_str.split()],
        dtype=bool,
    )

    if layer_1d.size != nx:
        raise ValueError(
            f"in FCC specification: layer {index_layer} "
            f"did not have {nx = } sites (had {layer_1d.size})"
        )

    layer_2d: np.ndarray = np.tile(layer_1d, (ny, 1)).T

    return layer_2d


def calc_box_size(nx: int, ny: int, nz: int, spacing: Spacing) -> tuple[float, float, float]:
    return (
        nx * spacing.dx,
        ny * spacing.dy,
        nz * spacing.dz,
    )


def calc_fcc_spacings(spacing: float) -> Spacing:
    dx = spacing
    dy = (np.sqrt(3.) / 2.) * spacing
    dz = (np.sqrt(6.) / 3.) * spacing

    dx_y = 0.5 * spacing
    dx_z = 0.5 * spacing
    dy_z = (1. / 3.) * spacing

    return Spacing(
        dx=dx,
        dy=dy,
        dz=dz,
        dx_y=dx_y,
        dx_z=dx_z,
        dy_z=dy_z,
    )


def gen_fcc_layer(
    x0: float,
    y0: float,
    z: float,
    spacing: Spacing,
    inds: np.ndarray,
    box_size: tuple[float, float, float],
) -> list[tuple[float, float, float]]:
    if inds.ndim != 2:
        raise ValueError(
            f"`inds` array must be 2-dimensional (was {inds.ndim}-dim)"
        )

    lattice_points = []

    nx, ny = inds.shape
    xs = x0 + spacing.dx * np.arange(nx)
    ys = y0 + spacing.dy * np.arange(ny)

    box_x, box_y, _ = box_size

    for ix in range(nx):
        for iy in range(ny):
            add_point = inds[ix, iy]

            if add_point:
                shift_x_y = (iy % 2) * spacing.dx_y

                x = np.fmod(xs[ix] + shift_x_y, box_x)
                y = np.fmod(ys[iy], box_y)

                lattice_points.append((x, y, z))

    return lattice_points


def map_atom_data_to_lattice_points(
    lattice_points: list[tuple[float, float, float]],
    sites: Sequence[AtomSpec],
    residue_name: str,
    translate: tuple[float, float, float]
) -> list[Atom]:
    atoms = []
    x0, y0, z0 = translate

    for x, y, z in lattice_points:
        x1 = x0 + x
        y1 = y0 + y
        z1 = z0 + z

        for site in sites:
            dx, dy, dz = site.dx
            x = (x1 + dx, y1 + dy, z1 + dz)

            atoms.append(Atom(x, None, site.name, residue_name))

    return atoms


def main():
    parser = ArgumentParser(
        description="Generate a FCC lattice from a specification file.",
    )

    parser.add_argument(
        'conf',
        help="path to FCC lattice specification file",
    )
    parser.add_argument(
        '-o', '--output',
        default=None, metavar='PATH',
        help="write lattice to file instead of stdout",
    )

    args = parser.parse_args()

    try:
        fcc_spec = read_fcc_spec(args.conf)
    except Exception as exc:
        print(f"error: could not read FCC specification from `{args.conf}` ({exc})")
        exit(1)

    spacing = calc_fcc_spacings(fcc_spec.spacing)
    x0 = spacing.dx / 4.
    y0 = spacing.dy / 6.
    z0 = spacing.dz / 2.

    nx = fcc_spec.nx
    ny = fcc_spec.ny
    nz = fcc_spec.nz
    box_size = calc_box_size(nx, ny, nz, spacing)

    lattice_points = []

    for iz in range(nz):
        try:
            inds = get_layer_inds_from_spec(iz, nx, ny, fcc_spec.layers)
        except ValueError as exc:
            stderr.write(f"error: {exc}\n")
            exit(1)

        # When the current layer's shift along x (or y) is >= spacing,
        # shift it to its "lowest" position. The added number 1e-7
        # is to handle floating point errors when the two are
        # "exactly" matching (dx = 0!).
        dx = np.fmod(iz * spacing.dx_z + 1e-7, spacing.dx)
        dy = np.fmod(iz * spacing.dy_z + 1e-7, spacing.dy)

        x = x0 + dx
        y = y0 + dy
        z = z0 + iz * spacing.dz

        lattice_points += gen_fcc_layer(x, y, z, spacing, inds, box_size)

    atoms = map_atom_data_to_lattice_points(
        lattice_points,
        fcc_spec.atoms,
        fcc_spec.residue_name,
        fcc_spec.translate
    )

    conf = Gromos87(
        title=fcc_spec.title,
        atoms=atoms,
        box_size=box_size,
    )

    if args.output:
        write_gromos87(args.output, conf)
    else:
        write_gromos87_stdout(conf)
