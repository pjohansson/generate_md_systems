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


def put_in_box(conf: Gromos87, residue_length: int) -> Gromos87:
    """Put all residues in the simulation box, effectively removing PBC."""

    def is_outside(position: Vec3) -> bool:
        x, y, z = position

        return (
            (x < 0.)
            or (y < 0.)
            or (z < 0.)
            or (x > box_x)
            or (y > box_y)
            or (z > box_z)
        )

    def rmpbc_shift(position: Vec3) -> Vec3:
        def rmpbc_1d(x: float, box: float) -> float:
            n = 0

            while x > box:
                n -= 1
                x -= box
            while x < 0.:
                n += 1
                x += box

            return float(n) * box

        x, y, z = position

        return Vec3(
            x=rmpbc_1d(x, box_x),
            y=rmpbc_1d(y, box_y),
            z=rmpbc_1d(z, box_z)
        )

    box_x, box_y, box_z = conf.box_size

    num_residues = len(conf.atoms) // residue_length

    atoms = []

    for i in range(num_residues):
        n = i * residue_length

        if is_outside(conf.atoms[n].position):
            # Use a common translation for all atoms of the same residue,
            # calculated from the "head" atom
            shift = rmpbc_shift(conf.atoms[n].position)

            for j in range(n, n + residue_length):
                atom = conf.atoms[j]
                position = Vec3(
                    x=atom.position.x + shift.x,
                    y=atom.position.y + shift.y,
                    z=atom.position.z + shift.z
                )

                atoms.append(Atom(
                    name=atom.name,
                    residue=atom.residue,
                    position=position,
                    velocity=atom.velocity
                ))

        else:
            for j in range(n, n + residue_length):
                atoms.append(conf.atoms[j])

    return Gromos87(
        title=conf.title,
        atoms=atoms,
        box_size=conf.box_size
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

    conf_rmpbc = put_in_box(conf, residue_length=residue_length)
    conf_stacked = stack_gromos87_conf_to_minimum_size(
        conf_rmpbc, size_x, size_y, size_z)
    conf_resized = cut_gromos87_conf_to_size(
        conf_stacked, size_x, size_y, size_z, residue_length=residue_length)

    if translate != None:
        conf_final = translate_gromos87(conf_resized, *translate)
    else:
        conf_final = conf_resized

    return conf_final
