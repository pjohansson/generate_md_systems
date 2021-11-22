from .io import Atom, Gromos87, Vec3

def translate_gromos87(conf, dx, dy, dz):
    def translate_atom(atom):
        x, y, z = atom.position

        atom = Atom(
            position = Vec3(x + dx, y + dy, z + dz),
            velocity=atom.velocity,
            name=atom.name,
            residue=atom.residue,
        )

        return atom

    if dx == None: dx = 0.
    if dy == None: dy = 0.
    if dz == None: dz = 0.
    
    atoms = [translate_atom(a) for a in conf.atoms]

    conf = Gromos87(
        title=conf.title,
        atoms=atoms,
        box_size=conf.box_size,
    )

    return conf