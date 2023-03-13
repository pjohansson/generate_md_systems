import numpy as np
from .io import Atom, Gromos87, Vec3


def translate_gromos87(
    conf: Gromos87,
    dx: float,
    dy: float,
    dz: float,
) -> Gromos87:
    def translate_atom(atom):
        x, y, z = atom.position

        atom = Atom(
            position=Vec3(x + dx, y + dy, z + dz),
            velocity=atom.velocity,
            name=atom.name,
            residue=atom.residue,
        )

        return atom

    if dx == None:
        dx = 0.
    if dy == None:
        dy = 0.
    if dz == None:
        dz = 0.

    atoms = [translate_atom(a) for a in conf.atoms]

    conf = Gromos87(
        title=conf.title,
        atoms=atoms,
        box_size=conf.box_size,
    )

    return conf


def stack_gromos87_conf_to_minimum_size(
    conf: Gromos87, 
    to_size_x: float, 
    to_size_y: float, 
    to_size_z: float,
) -> Gromos87:
    """Stack multiples of the configuration to a minimum final size.

    Note that this only stacks the boxes, it does not cut the 
    configuration to the correct size. It only guarantuees that
    the output size is at least this large.

    """

    def move_atom(atom, ix, iy, iz, dx, dy, dz):
        x = atom.position.x + ix * dx
        y = atom.position.y + iy * dy
        z = atom.position.z + iz * dz

        return Atom(
            name=atom.name,
            residue=atom.residue,
            position=Vec3(x, y, z),
            velocity=atom.velocity,
        )

    size_x, size_y, size_z = conf.box_size

    nx = int(np.ceil(to_size_x / size_x))
    ny = int(np.ceil(to_size_y / size_y))
    nz = int(np.ceil(to_size_z / size_z))

    atoms = conf.atoms.copy()

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                if ix > 0 or iy > 0 or iz > 0:
                    atoms += [
                        move_atom(atom, ix, iy, iz, size_x, size_y, size_z)
                        for atom in conf.atoms]

    box_size = size_x * nx, size_y * ny, size_z * nz

    return Gromos87(
        title=conf.title,
        box_size=box_size,
        atoms=atoms,
    )


def cut_gromos87_conf_to_size(
    conf: Gromos87, 
    xmax: float, 
    ymax: float, 
    zmax: float, 
    residue_length: int = 1,
) -> Gromos87:
    """Cut a configuration to the set size.

    Molecules are kept if any of their atoms lie inside the box.

    """

    def check_atom(atom, xmax, ymax, zmax):
        return (atom.position.x <= xmax
                and atom.position.y <= ymax
                and atom.position.z <= zmax)

    # Keep molecule if any atom is inside, else discard
    def check_molecule(atoms, xmax, ymax, zmax):
        return np.any(
            [check_atom(atom, xmax, ymax, zmax) for atom in atoms])

    num_mols = len(conf.atoms) // residue_length
    kept_atoms = []

    for i in range(num_mols):
        i0 = i * residue_length
        i1 = (i + 1) * residue_length

        atoms = conf.atoms[i0:i1]

        if check_molecule(atoms, xmax, ymax, zmax):
            kept_atoms += atoms

    box_size = xmax, ymax, zmax

    return Gromos87(
        title=conf.title,
        box_size=box_size,
        atoms=kept_atoms
    )


def create_gromos87_conf_with_size(
    conf: Gromos87, 
    size_x: float, 
    size_y: float, 
    size_z: float, 
    residue_length: int = 1, 
    translate: tuple[float, float, float] | None = None,
) -> Gromos87:
    """Expand or cut the given configuration to the set size."""

    conf_stacked = stack_gromos87_conf_to_minimum_size(
        conf, size_x, size_y, size_z)
    conf_resized = cut_gromos87_conf_to_size(
        conf_stacked, size_x, size_y, size_z, residue_length=residue_length)

    if translate != None:
        conf_final = translate_gromos87(conf_resized, *translate)
    else:
        conf_final = conf_resized

    return conf_final
