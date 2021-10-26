import os
from collections import namedtuple
from sys import stderr, stdout

Vec3 = namedtuple('Vec3', ['x', 'y', 'z'])
Atom = namedtuple('Atom', ['position', 'velocity', 'name', 'residue'])
Gromos87 = namedtuple('Gromos87Conf', ['title', 'atoms', 'box_size'])


def read_gromos87(fn):
    """Read a Gromos87 formatted atomic configuration file."""

    def read_atom_line(line):
        residue = line[5:10].strip()
        atom_name = line[10:15].strip()

        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])

        position = Vec3(x, y, z)

        if line[44:].strip() != "":
                vx = float(line[44:52])
                vy = float(line[52:60])
                vz = float(line[60:68])
                velocity = Vec3(vx, vy, vz)
        else:
            velocity = None

        return Atom(position, velocity, name=atom_name, residue=residue)

    with open(fn) as fp:
        title = fp.readline().strip()
        num_atoms = int(fp.readline().strip())
        atoms = [read_atom_line(fp.readline()) for _ in range(num_atoms)]
        box_size = Vec3(*[float(v) for v in fp.readline().strip().split()])

    return Gromos87(title, atoms, box_size)


def write_gromos87(fn, conf):
    """Write a configuration to a Gromos87 formatted file."""

    def backup_file(fn):
        original_fn = fn

        dirname, basename = os.path.split(fn)

        i = 1

        while os.path.isfile(fn):
            current_basename = "#{}.{}#".format(basename, i)
            fn = os.path.join(dirname, current_basename)
            i += 1

        if fn != original_fn:
            os.rename(original_fn, fn)
            stderr.write("backed up '{}' to '{}'\n".format(original_fn, fn))

        return

    backup_file(fn)

    with open(fn, 'w') as fp:
        _write_gromos87_to_filepointer(fp, conf)


def write_gromos87_stdout(conf):
    """Write a configuration to stdout in Gromos87 format."""

    _write_gromos87_to_filepointer(stdout, conf)


def write_index_group(fp, group_name, indices, begin=None, end="\n"):
    """Write an index group to a buffer as a Gromacs index format."""

    if type(begin) == str:
        fp.write(begin)

    fp.write("[ {} ]\n".format(group_name))

    for i, index in enumerate(indices):
        fp.write("{} ".format(index))

        if i > 0 and i % 15 == 0:
            fp.write("\n")

    if i % 15 != 0:
        fp.write("\n")

    if type(end) == str:
        fp.write(end)


def _write_gromos87_to_filepointer(fp, conf):
    fp.write("{}\n".format(conf.title))
    fp.write("{:5}\n".format(len(conf.atoms)))

    for i, atom in enumerate(conf.atoms):
        residue_num = (i + 1) % 100000
        atom_num = residue_num

        fp.write("{:>5}{:<5}{:>5}{:>5}".format(
            residue_num, atom.residue, atom.name, atom_num))

        x, y, z = atom.position
        fp.write("{:8.3f}{:8.3f}{:8.3f}".format(x, y, z))

        if atom.velocity != None:
            vx, vy, vz = atom.velocity
            fp.write("{:8.4f}{:8.4f}{:8.4f}".format(vx, vy, vz))

        fp.write("\n")

    dx, dy, dz = conf.box_size
    fp.write("{:10.5f}{:10.5f}{:10.5f}\n".format(dx, dy, dz))

