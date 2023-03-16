#!/usr/bin/env python3

from argparse import ArgumentParser 
from gmxutils import *

if __name__ == '__main__':
    parser = ArgumentParser() 

    parser.add_argument('file', 
        type=str, metavar='PATH', 
        help=".gro file to create restraints for")

    parser.add_argument('-r', '--residue', 
        type=str, default='CUB', metavar='RES',
        help="residue name of restrained molecules")

    args = parser.parse_args()

    
