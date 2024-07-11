#!/usr/bin/env python
# coding: utf-8

"""workflow2_make_perturbations.py: Creates hybrid structures and topologies."""

import os
from pmx.workflow import pmxworkflow
import argparse
import subprocess


from openff.toolkit.topology import Molecule

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def atomsToMorph(pwf):
    pmxworkflow.printInfo(runtype='atomsToMorph',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        lig1 = item['ligand_a']
        lig2 = item['ligand_b']

        # make temporary directory tmp
        os.makedirs(f'{pwf.hybPath}/{edge}/water/tmp/', exist_ok=True)

        # input
        pdb1 = f'{pwf.hybPath}/{edge}/water/tmp/{lig1}.pdb'  # pdb structure lig1
        pdb2 = f'{pwf.hybPath}/{edge}/water/tmp/{lig2}.pdb'  # pdb structure lig2

        # conversion from sdf to pdb
        ligPath1 = f'{pwf.ligPath}/{lig1}/crd/{lig1}.sdf'
        ligand = Molecule.from_file(f'{ligPath1}', allow_undefined_stereo=True)
        ligand.to_file(pdb1, 'pdb')

        ligPath2 = f'{pwf.ligPath}/{lig2}/crd/{lig2}.sdf'
        ligand = Molecule.from_file(f'{ligPath2}', allow_undefined_stereo=True)
        ligand.to_file(pdb2, 'pdb')

        # output
        o1 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb1.pdb'  # processed structure lig1; temporary file
        o2 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb2.pdb'  # processed structure lig2; temporary file
        pairs = f'{pwf.hybPath}/{edge}/water/crd/pairs.dat'  # atom mapping
        score = f'{pwf.hybPath}/{edge}/water/crd/score.dat'  # score of the mapping

        # make coordinate directory crd
        os.makedirs(f'{pwf.hybPath}/{edge}/water/crd/', exist_ok=True)

        command = ['python3', pwf.scriptpath + '/atoms_to_morph.py',
                                    '-i1', pdb1,
                                    '-i2', pdb2,
                                    '-opdb1', o1,
                                    '-opdb2', o2,
                                    '-score', score,
                                    '-o', pairs,
                                    '-H2H',
                                    '-alignment',
                                    '-timeout', str(30)]
        print(f'COMMAND: {" ".join(command)}')
        # run atoms_to_morph.py script
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


# create hybrid structures/topologies
def makeHybrid(pwf):
    pmxworkflow.printInfo(runtype='makeHybrid',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        lig1 = item['ligand_a']
        lig2 = item['ligand_b']
        # make temporary directory tmp
        os.makedirs(f'{pwf.hybPath}/{edge}/water/tmp/', exist_ok=True)
        # input
        l1 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb1.pdb'  # structure lig1; temporary file
        l2 = f'{pwf.hybPath}/{edge}/water/tmp/out_pdb2.pdb'  # structure lig2; temporary file
        pairs = f'{pwf.hybPath}/{edge}/water/crd/pairs.dat'  # atom mapping
        itp1 = f'{pwf.ligPath}/{lig1}/top/{pwf.forcefield}/{lig1}.itp'  # topology lig1
        itp2 = f'{pwf.ligPath}/{lig2}/top/{pwf.forcefield}/{lig2}.itp'  # topology lig2
        # output
        opdbA = f'{pwf.hybPath}/{edge}/water/crd/mergedA.pdb'  # hybrid structure based on lig1
        opdbB = f'{pwf.hybPath}/{edge}/water/crd/mergedB.pdb'  # hybrid structure based on lig2
        oitp = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/merged.itp' # hybrid topology
        offitp = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffmerged.itp'  # force field parameters for dummies
        olog = f'{pwf.hybPath}/{edge}/water/tmp/hybrid.log'  # output log


        if not os.path.exists(l1):
            ligPath1 = f'{pwf.ligPath}/{lig1}/crd/{lig1}.sdf'
            ligand = Molecule.from_file(f'{ligPath1}', allow_undefined_stereo=True)
            ligand.to_file(l1, 'pdb')
        if not os.path.exists(l2):
            ligPath2 = f'{pwf.ligPath}/{lig2}/crd/{lig2}.sdf'
            ligand = Molecule.from_file(f'{ligPath2}', allow_undefined_stereo=True)
            ligand.to_file(l2, 'pdb')

        # make top directory
        os.makedirs(f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/', exist_ok=True)

        command = ['python', pwf.scriptpath + '/make_hybrid.py',
                                    '-pairs', pairs,
                                    '-l1', l1,
                                    '-l2', l2,
                                    '-itp1', itp1,
                                    '-itp2', itp2,
                                    '-oitp', oitp,
                                    '-oa', opdbA,
                                    '-ob', opdbB,
                                    '-ffitp', offitp,
                                    '-log', olog,
                                    '-scDUMm', str(0.001)]

        print(f'COMMAND: {" ".join(command)}')
        # creates hybrid structure/topology
        process = subprocess.Popen(command, 
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


# deal with the atomtype topology files
def oneffFile(pwf):
    pmxworkflow.printInfo(runtype='oneffFile',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        lig1 = item['ligand_a']
        lig2 = item['ligand_b']

        # input
        ff1 = f'{pwf.ligPath}/{lig1}/top/{pwf.forcefield}/ff{lig1}.itp'  # atomtypes lig1
        ff2 = f'{pwf.ligPath}/{lig2}/top/{pwf.forcefield}/ff{lig2}.itp'  # atomtypes lig2
        ffmerged = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffmerged.itp'  # force field parameters for dummies
        # output
        ffout = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffMOL.itp'  # atomtypes out

        command = ['python', pwf.scriptpath + '/one_ff_file.py',
                                    '-ffitp', ff1, ff2, ffmerged,
                                    '-ffitp_out', ffout]

        print(f'COMMAND: {" ".join(command)}')
        process = subprocess.Popen(command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


# clean unnecessary files
def cleanUp(pwf):
    pmxworkflow.printInfo(runtype='Cleaning up',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')
        # clean up
        toclean = os.listdir(f'{pwf.hybPath}/{edge}/water/tmp/')
        for file in toclean:
            os.remove(f'{pwf.hybPath}/{edge}/water/tmp/{file}')
        os.rmdir(f'{pwf.hybPath}/{edge}/water/tmp/')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',
                        '--target',
                        metavar='TARGET',
                        type=str,
                        default='01_jnk1',
                        help='The target protein.')
    parser.add_argument('-f',
                        '--forcefield',
                        metavar='FORCEFIELD',
                        type=str,
                        default='openff-2.0.0-rc.2.offxml',
#                        choices=['smirnoff99Frosst-1.1.0.offxml', 'openff-1.0.0.offxml', 'openff-1.2.0.offxml', 'gaff2', 'openff-2.0.0-rc.2.offxml'],
                        help='The force field used.')
    parser.add_argument('-p',
                        '--path',
                        metavar='PATH',
                        type=str,
                        default='./',
                        help='The path to the data.')
    parser.add_argument('-e',
                        '--edges',
                        metavar='EDGES',
                        nargs='+',
                        type=str,
                        default=['all'],
                        help='Either "all" or a list of edges to be calculated.')
    parser.add_argument('-v',
                        '--verbose',
                        type=bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Turn on verbose output.')
    args = parser.parse_args()

    pwf = pmxworkflow.pmxvariables(target=args.target,
                                   forcefield=args.forcefield,
                                   path=args.path,
                                   replicates=[1],
                                   verbose=args.verbose)

    if args.edges[0] != 'all':
        edges = dict()
        for edge in args.edges:
            edges[edge] = pwf.edges[edge]
        pwf.edges = edges

    atomsToMorph(pwf)
    makeHybrid(pwf)
    oneffFile(pwf)
    cleanUp(pwf)
    print("Process completed. Make sure to parameterize your protein and cofactors before executing workflow 3. "\
    "The topology files should be stored in the 01/protein directory as 'protein.top' and 'cofactors_crystalwater.top'. "\
    "These files will be used to parse the number of molecules present under [ molecules ]. Also make separate itp files"\
    " for protein (and all cofactors if necessary) in the same directory. Finally, replace the"\
    " 'protein.pdb' file in the 01/protein directory by the one generated by pdb2gmx (with hydrogens).")

