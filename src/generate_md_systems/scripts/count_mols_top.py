import tomllib

from argparse import ArgumentParser
from sys import exit, stderr
from typing import Any, Mapping

from generate_md_systems.gmx_conf_utils import read_gromos87


DEFAULT_TOPOL_SPEC = {
    'SOL': {'length': 3},
}


def create_record(
    residue: str,
    count: int,
    residue_length: int,
    spec: Mapping[str, Any],
    use_default_spec: bool = True,
) -> tuple[str, int]:
    def read_value(field: str, default: Any) -> Any:
        try:
            value = spec[residue][field]
        except Exception as exc:
            if use_default_spec:
                try:
                    value = DEFAULT_TOPOL_SPEC[residue][field]
                except Exception:
                    value = default
            else:
                value = default

        return value

    name = read_value('name', residue)
    count = count // read_value('length', residue_length)

    return name, count


def main():
    parser = ArgumentParser()

    parser.add_argument(
        'conf',
        help='.gro file to count .top [molecules] in',
    )

    parser.add_argument(
        '-l', '--residue-length',
        default=1, type=int, metavar='N',
        help='number of atoms per residue (default: %(default)s)',
    )
    parser.add_argument(
        '-m', '--prepend-molecule-directive',
        action='store_true',
        help='prepend the counter with the `[ molecules ]` directive',
    )
    parser.add_argument(
        '-s', '--specification',
        default=None, metavar='PATH',
        help='path to file with residue specification in .toml format',
    )
    parser.add_argument(
        '-i', '--ignore-default-spec',
        action='store_false', dest='use_default_spec',
        help='do not read specifications from the default dictionary',
    )

    args = parser.parse_args()

    conf = read_gromos87(args.conf)

    if len(conf.atoms) == 0:
        stderr.write(f"note: no atoms in `{args.conf}`")
        exit(0)

    if args.specification is not None:
        try:
            with open(args.specification, 'rb') as fp:
                spec = tomllib.load(fp)
        except Exception as exc:
            print(
                f"error: could not read specification in "
                f"`{args.specification}` ({exc})"
            )
            exit(1)
    else:
        spec = {}

    current_residue = conf.atoms[0].residue
    current_num = 0
    counter = []

    def add_to_counter():
        record = create_record(
            current_residue,
            current_num,
            args.residue_length,
            spec,
            use_default_spec=args.use_default_spec,
        )

        counter.append(record)

    for atom in conf.atoms:
        if atom.residue == current_residue:
            current_num += 1
        else:
            add_to_counter()

            current_residue = atom.residue
            current_num = 1

    add_to_counter()

    max_len_residue = 0
    max_len_count = 0

    for residue, count in counter:
        if len(residue) > max_len_residue:
            max_len_residue = len(residue)

        if len(str(count)) > max_len_count:
            max_len_count = len(str(count))

    if args.prepend_molecule_directive:
        print(f"[ molecules ]")

    for residue, count in counter:
        print(f"{residue:{max_len_residue}} {count:{max_len_count}}")
