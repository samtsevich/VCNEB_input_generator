'''
Code for creating linear interpolation between 2 structures.


IMPORTANT:
initial and final structures should be POSCAR-files of VASP5
'''

from __future__ import division

__author__ = 'asamtsevich'

from ase.atoms import Atoms
from ase.io.vasp import read_vasp, write_vasp

import argparse
import os
import numpy as np
import sys



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCNEB_tool')

    # parser.add_argument("-v", "--version", dest="version", action="store_true",
    #                     help="show program's version number and exit")

    parser.add_argument('-i', dest='initStructure', type=str,
                        help='path to initial POSCAR file')

    parser.add_argument('-f', dest='finalStructure', type=str,
                        help='path to final POSCAR file')

    parser.add_argument('-o', dest='output', type=str,
                        help='path to output file')

    parser.add_argument('-N', dest='N', type=int,
                        help='num of images that will be in output file')

    # parser.add_argument("-p", "--parameter", action='callback', dest="parm",
    #                     help="specify parameter to get help. If no value or 'all' value is specified, all INPUT.txt parameters will be shown",
    #                     metavar="PARM")


    args = parser.parse_args()

    if not args.initStructure:
        print('Please provide initStructure path')
        parser.print_help()
    if not args.finalStructure:
        print('Please provide finalStructure path')
        parser.print_help()
    if not args.output:
        print('Please provide output path')
        parser.print_help()
    if args.N <= 2:
        print('Please provide number of images > 2')
        parser.print_help()

    initStructure = read_vasp(args.initStructure)
    finalStructure = read_vasp(args.finalStructure)

    symbols = initStructure.get_chemical_symbols()

    initCoord = initStructure.get_scaled_positions()
    initCell = initStructure.get_cell()

    assert len(initStructure) == len(finalStructure)
    N_system = len(initStructure)

    N = args.N - 1

    diffCoord = finalStructure.get_scaled_positions() - initCoord
    for i in range(N_system):
        for j in range(0,3):
            if diffCoord[i, j] > 0.5:
                initCoord[i,j] += 1.0
            elif diffCoord[i, j] < -0.5:
                initCoord[i, j] -= 1.0

    diffCoord = finalStructure.get_scaled_positions() - initCoord

    assert (np.abs(diffCoord) < 0.5).all()

    # for atom in diffCoord:
    #     for x in atom:
    #         if x > 0.5:
    #             x -= 1.0
    #         elif x < -0.5:
    #             x += 1.0
    # stepCoord = diffCoord / N

    stepCoord = diffCoord / N
    stepCell = (finalStructure.get_cell() - initCell) / N

    TMP_STRUCTURE_PATH = 'tmp.POSCAR'

    with open(args.output, 'w') as fp:
        write_vasp(TMP_STRUCTURE_PATH, initStructure, label='Image_ini ', direct=True, vasp5=True)
        with open(TMP_STRUCTURE_PATH, 'r') as f:
            fp.write(f.read())

        for i in range(N-1):
            current = Atoms(symbols=symbols, scaled_positions=initCoord + (i+1)*stepCoord, cell=initCell + (i+1)*stepCell)
            write_vasp(TMP_STRUCTURE_PATH, current, label='Image ' + str(i+2), direct=True, vasp5=True)
            with open(TMP_STRUCTURE_PATH, 'r') as f:
                fp.write(f.read())

        write_vasp(TMP_STRUCTURE_PATH, finalStructure, label='Image_end ', direct=True, vasp5=True)
        with open(TMP_STRUCTURE_PATH, 'r') as f:
            fp.write(f.read())

    # cleaning unnecessary files
    os.remove(TMP_STRUCTURE_PATH)

    print('Building complete.')
