#!/usr/bin/env python
# coding: utf-8

"""workflow3_solvate.py: Solvates simulation scripts and add ions."""

import os
import glob
import shutil
from pmx.workflow import pmxworkflow
import argparse
import subprocess

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

################################################
# prepare water boxes for simulations in water #
################################################
def solvateLigand(pwf):
    pmxworkflow.printInfo(runtype='solvateLigand',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        lig1 = item['ligand_a']
        lig2 = item['ligand_b']
        # make temporary directory tmp
        os.makedirs(f'{pwf.hybPath}/{edge}/water/tmp/', exist_ok=True)

        # input
        topology = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/topol.top'
        pmxworkflow.create_top(fname=topology, itp=['ffMOL.itp', 'merged.itp'])

        #copy force field files
        forcefield_path = f'{pwf.ffPath}/mutff45/amber14sb_OL15.ff'
        top_path = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/amber14sb_OL15.ff'
        try:
            shutil.copytree(forcefield_path, top_path)
        except FileExistsError:
            pass
 
        # create box
        # input
        inp = f'{pwf.hybPath}/{edge}/water/crd/mergedA.pdb'  # input ligand structure
        # output
        box = f'{pwf.hybPath}/{edge}/water/tmp/box.pdb'  # ligand placed in box; temporary file
        process = subprocess.Popen(['gmx', 'editconf',
                                    '-f', inp,
                                    '-o', box,
                                    '-bt', 'dodecahedron',
                                    '-d', str(1.5)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))

        # solvate
        # output
        water = f'{pwf.hybPath}/{edge}/water/tmp/water.pdb'  # solvated ligand; temporary file
        process = subprocess.Popen(['gmx', 'solvate',
                                    '-cp', box,
                                    '-o', water,
                                    '-cs', 'spc216.gro',
                                    '-scale', str(1.0),
                                    '-p', topology],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))

        # add ions, 2 sub-steps: a) grompp, b) add ions for 3 replicas independently
        # grompp
        mdpemA = f'{pwf.mdpPath}/em_stateA.mdp'
        mdout = f'{pwf.hybPath}/{edge}/water/tmp/mdout.mdp'  # temporary mdout file (not important)
        tprfile = f'{pwf.hybPath}/{edge}/water/tmp/tpr.tpr'  # temporary tpr file
        process = subprocess.Popen(['gmx', 'grompp',
                                    '-p', topology,
                                    '-c', water,
                                    '-r', water,
                                    '-o', tprfile,
                                    '-f', mdpemA,
                                    '-po', mdout,
                                    '-maxwarn', str(2)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDERR{} '.format(out[1].decode("utf-8")))
            print('STDOUT{} '.format(out[0].decode("utf-8")))


################################################
# add ions for water boxes #
################################################
def addIonsLigand(pwf):
    pmxworkflow.printInfo(runtype='addIonsLigand',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge in pwf.edges.keys():
        topology = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/topol.top'
        tprfile = f'{pwf.hybPath}/{edge}/water/tmp/tpr.tpr'  # temporary tpr file

        # add ions
        for i in pwf.replicates:
            ions = f'{pwf.hybPath}/{edge}/water/crd/ions{i}.pdb'  # ion file
            topfile = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/topol{i}.top'  # top file
            try:
                shutil.copy(topology, topfile)
            except FileExistsError:
                pass
            process = subprocess.Popen(['gmx', 'genion',
                                        '-s', tprfile,
                                        '-p', topfile,
                                        '-conc', str(0.15),
                                        '-neutral',
                                        '-nname', 'CL',
                                        '-pname', 'NA',
                                        '-o', ions],
                                        stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            process.stdin.write('SOL'.encode())
            out = process.communicate()

            if pwf.verbose:
                print('STDERR{} '.format(out[1].decode("utf-8")))
                print('STDOUT{} '.format(out[0].decode("utf-8")))

def assembleComplex(file, molfiles=[]):
    fpout = open(file, 'w')
    for molfile in molfiles:
        fp = open(molfile, 'r')
        lines = fp.readlines()
        fp.close()

        for l in lines:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                fpout.write(l)
    fpout.close()


################################################
# assemble protein and ligand                  #
################################################
def combineProteinLigand(pwf):
    pmxworkflow.printInfo(runtype='Combine Protein and Ligand',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        # make directory for complex
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/crd/', exist_ok=True)
        # make temporary directory
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/tmp/', exist_ok=True)
        # make directory for complex topology
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}', exist_ok=True)

        # assemble coordinates
        # input
        proteinFile = f'{pwf.protPath}/crd/protein.pdb'  # input protein structure with ions, cofactors and crystal waters
        ligandFile = f'{pwf.hybPath}/{edge}/water/crd/mergedA.pdb'  # input ligand structure
        waterFile = f'{pwf.protPath}/crd/cofactors_crystalwater.pdb'  # input protein structure with ions, cofactors and crystal waters
        #output
        complexFile = f'{pwf.hybPath}/{edge}/complex/crd/complex.pdb' # complex of the former two structures

        if os.path.isfile(ligandFile):
            assembleComplex(complexFile, molfiles=[proteinFile, ligandFile, waterFile])
        else:
            print(f"No hybrid structure available. File {ligandFile} missing.")
            continue
 
        # assemble topologies
        forcefield_path = f'{pwf.ffPath}/mutff45/amber14sb_OL15.ff'
        top_path = f'{pwf.protPath}/top/amber14sb_OL15.ff'
        try:
            shutil.copytree(forcefield_path, top_path)
        except FileExistsError:
            pass
        PToporig = f'{pwf.protPath}/protein.top'
        CCToporig = f'{pwf.protPath}/cofactors_crystalwater.top'
        proteinTop = f'{pwf.protPath}/top/amber14sb_OL15.ff/protein.top' # protein topol        
        cofactors_crystalwaterTop = f'{pwf.protPath}/top/amber14sb_OL15.ff/cofactors_crystalwater.top' # protein topology
        try:
            shutil.copy(PToporig, proteinTop)
            shutil.copy(CCToporig, cofactors_crystalwaterTop)
        except FileExistsError:
            pass
        # copy itp files
        for filename in os.listdir(pwf.protPath):
            if filename.endswith('.itp'):
                filepath = os.path.join(pwf.protPath, filename)
                shutil.copy(filepath, top_path)

        # output
        complexTop = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top'  # complex topology

        if os.path.exists(f'{pwf.protPath}/top/amber14sb_OL15.ff/ffcofactors_crystalwater.itp'):
            # generate oneff_fuke
            # input
            ff1 = f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/ffMOL.itp'
            ff2 = f'{pwf.protPath}/top/amber14sb_OL15.ff/ffcofactors_crystalwater.itp'
            # output
            ffout = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/ffcomplex.itp'


            command = ['python', pwf.scriptpath + '/one_ff_file.py',
                       '-ffitp', ff1, ff2,
                       '-ffitp_out', ffout]

            print(f'COMMAND: {" ".join(command)}')
            process = subprocess.Popen(command,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.wait()

            if pwf.verbose:
                out = process.communicate()
                print('STDERR{} '.format(out[1].decode("utf-8")))
                print('STDOUT{} '.format(out[0].decode("utf-8")))

            # included itp files
            # merged molecules itp
            includes = ['ffcomplex.itp', 'merged.itp']
        else:
            # included itp files
            # merged molecules itp
            includes = ['ffMOL.itp', 'merged.itp']

        # protein
        includes += list(filter(lambda x: x.endswith('.itp') and not x.startswith('ff') and not x.startswith('tip') and not x.startswith('ions') and not x.startswith('spc') and not x.startswith('posre') and not x.startswith('gbsa') and not x.startswith('forcefield'), os.listdir(f'{pwf.protPath}/top/amber14sb_OL15.ff/')))

        # molecules to be included
        molecules = []
        # protein system parts
        with open(proteinTop, 'r') as fp:
            lines = fp.readlines()
            for i, l in enumerate(lines):
                if l.startswith('[ molecules ]'):
                    break
            for l in lines[i+1:]:
                if l.startswith(';'):
                    continue
                molecules.append(l.split())
        # add ligand molecule
        molecules.append(['MOL', 1])
        # add cofactor and water system parts
        with open(cofactors_crystalwaterTop, 'r') as fp:
            lines = fp.readlines()
            for i, l in enumerate(lines):
                if l.startswith('[ molecules ]'):
                    break
            for l in lines[i+1:]:
                if l.startswith(';'):
                    continue
                molecules.append(l.split())

        pmxworkflow.create_top(fname=f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top', ff='amber14sb_OL15.ff', water='tip3p',
                   itp=includes, mols=molecules,
                   destination=f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/', toppaths=[f'{pwf.protPath}/top/amber14sb_OL15.ff/', f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/', f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/'])

################################################
# prepare water boxes for complex simulations  #
################################################
def solvateComplex(pwf):
    pmxworkflow.printInfo(runtype='solvateComplex',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        # make directory for complex
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/crd/', exist_ok=True)
        # make temporary directory
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/tmp/', exist_ok=True)

        forcefield_path = f'{pwf.protPath}/top/amber14sb_OL15.ff'
        top_path = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/amber14sb_OL15.ff'
        try:
            shutil.copytree(forcefield_path, top_path)
        except FileExistsError:
            pass

        # create box
        # input
        inp = f'{pwf.hybPath}/{edge}/complex/crd/complex.pdb'  # input complex structure
        topology = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top'
        # output
        box = f'{pwf.hybPath}/{edge}/complex/tmp/box.pdb'  # complex placed in box, temporary file
        process = subprocess.Popen(['gmx', 'editconf',
                                    '-f', inp,
                                    '-o', box,
                                    '-bt', 'dodecahedron',
                                    '-d', str(1.5)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDOUT{} '.format(out[0].decode("utf-8")))
            print('STDERR{} '.format(out[1].decode("utf-8")))

        # solvate
        water = f'{pwf.hybPath}/{edge}/complex/tmp/water.pdb'  # solvated complex, temporary file
        process = subprocess.Popen(['gmx', 'solvate',
                                    '-cp', box,
                                    '-o', water,
                                    '-cs', 'spc216.gro',
                                    '-scale', str(1.0),
                                    '-p', topology],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDOUT{} '.format(out[0].decode("utf-8")))
            print('STDERR{} '.format(out[1].decode("utf-8")))

        # grompp
        mdpemA = f'{pwf.mdpPath}/em_stateA.mdp'
        mdout = f'{pwf.hybPath}/{edge}/complex/tmp/mdout.mdp'  # temporary mdout file (not important)
        tprfile = f'{pwf.hybPath}/{edge}/complex/tmp/tpr.tpr'  # temporary tpr file
        process = subprocess.Popen(['gmx', 'grompp',
                                    '-p', topology,
                                    '-c', water,
                                    '-r', water,
                                    '-o', tprfile,
                                    '-f', mdpemA,
                                    '-po', mdout,
                                    '-maxwarn', str(3)],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out = process.communicate()

        if pwf.verbose:
            print('STDOUT{} '.format(out[0].decode("utf-8")))
            print('STDERR{} '.format(out[1].decode("utf-8")))


################################################
# add ions for protein boxes #
################################################
def addIonsComplex(pwf):
    pmxworkflow.printInfo(runtype='addIonsComplex',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge in pwf.edges.keys():
        # make directory for complex
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/', exist_ok=True)
        # make temporary directory
        os.makedirs(f'{pwf.hybPath}/{edge}/complex/tmp/', exist_ok=True)

        topology = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol.top'
        tprfile = f'{pwf.hybPath}/{edge}/complex/tmp/tpr.tpr'  # temporary tpr file
        if not os.path.isfile(topology):
            print(f"Topology file {topfile} does not exist.")
            continue

        # add ions
        for i in pwf.replicates:
            ions = f'{pwf.hybPath}/{edge}/complex/crd/ions{i}.pdb'  # ion file
            topfile = f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/topol{i}.top'  # top file
            shutil.copy(topology, topfile)

            process = subprocess.Popen(['gmx', 'genion',
                                        '-s', tprfile,
                                        '-p', topfile,
                                        '-conc', str(0.15),
                                        '-neutral',
                                        '-nname', 'CL',
                                        '-pname', 'NA',
                                        '-o', ions],
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.stdin.write('SOL'.encode())

            out = process.communicate()

            if pwf.verbose:
                print('STDOUT{} '.format(out[0].decode("utf-8")))
                print('STDERR{} '.format(out[1].decode("utf-8")))

# clean unnecessary files
def cleanUp(pwf):
    pmxworkflow.printInfo(runtype='Cleaning up',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')
    for edge, item in pwf.edges.items():
        print(f'    - {edge}')

        # clean up
        toclean = glob.glob(f'{pwf.hybPath}/{edge}/water/tmp/*')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/complex/tmp/*')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/water/crd/*#')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/water/top/{pwf.forcefield}/*#')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/complex/crd/*#')
        toclean += glob.glob(f'{pwf.hybPath}/{edge}/complex/top/{pwf.forcefield}/*#')
        for clean in toclean:
            os.remove(clean)
        os.rmdir(f'{pwf.hybPath}/{edge}/water/tmp/')
        os.rmdir(f'{pwf.hybPath}/{edge}/complex/tmp/')



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
                        help='The force field used.')
    parser.add_argument('-p',
                        '--path',
                        metavar='PATH',
                        type=str,
                        default='./',
                        help='The path to the data.')
    parser.add_argument('-r',
                        '--replicates',
                        metavar='REPLICATES',
                        nargs='+',
                        type=str,
                        default='1',
                        help='Comma or space separated list of integers or range of integers denoted by hyphen')
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

    replicates = args.replicates
    replicates = [n for r in replicates for n in r.split(',') if n is not '']
    replicates = [n for r in replicates for n in ' - '.join(r.split('-')).split()]
    replicatesToUse = []
    for i, r in enumerate(replicates):
        if r == '-':
            replicatesToUse += range(int(replicates[i - 1]) + 1, int(replicates[i + 1]), 1)
        else:
            replicatesToUse.append(int(r))

    pwf = pmxworkflow.pmxvariables(target=args.target,
                                   forcefield=args.forcefield,
                                   path=args.path,
                                   replicates=replicates,
                                   verbose=args.verbose)

    if args.edges[0] != 'all':
        edges = dict()
        for edge in args.edges:
            edges[edge] = pwf.edges[edge]
        pwf.edges = edges

    solvateLigand(pwf)
    addIonsLigand(pwf)

    combineProteinLigand(pwf)

    solvateComplex(pwf)
    addIonsComplex(pwf)

    cleanUp(pwf)


