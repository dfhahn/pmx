#!/usr/bin/env python
# coding: utf-8

"""
pmxworkflow.py
A collection of notebooks and scripts for calculating free energies with pmx.

Handles the primary functions
"""

import os
import shutil

from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

# from openforcefield.utils import toolkits

# ### OpenEye version: uncomment the following if you have and if you want to use the OpenEye toolkit, then RDKit and Ambertools toolkits
# ### Attention: there could be problems in the workflow due to atom order changes, this is related to this issue: https://github.com/openforcefield/openforcefield/issues/475
# #toolkit_precedence = [toolkits.OpenEyeToolkitWrapper, toolkits.RDKitToolkitWrapper, toolkits.AmberToolsToolkitWrapper]

# ### Non-OpenEye version: uncomment the following if you want to use the rdkit and ambertools
# toolkit_precedence = [toolkits.RDKitToolkitWrapper, toolkits.AmberToolsToolkitWrapper]

# toolkits.GLOBAL_TOOLKIT_REGISTRY = toolkits.ToolkitRegistry(toolkit_precedence=toolkit_precedence)

import pmx
from plbenchmark import targets, ligands, edges


__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"



def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote

def infoString(runtype='em', run=1, target='xxx', edge='yyy_zzz', wc='water', state='stateA'):
    '''
    Generates an formatted string with information about which simulation is currently processsed.
    :param runtype: type of simulation
    :param run: number of repitition
    :param target: target name
    :param edge: edge name
    :param wc: water or complex simulation
    :param state: endstate of simulation
    :return: formatted information string
    '''
    outputString = '=' * 86 + '\n'
    outputString += f'=== {runtype:20}{str(run):<2s} {target:8s} {edge:29s} {wc:8s} {state:7s} ===\n'
    outputString += '=' * 86 + '\n'
    return outputString

def printInfo(runtype='em', run=1, target='xxx', edge='yyy_zzz', wc='water', state='stateA'):
    '''
    Prints information about which simulation is currently processsed to standard output.
    :param runtype: type of simulation
    :param run: number of repitition
    :param target: target name
    :param edge: edge name
    :param wc: water or protein simulation
    :param state: endstate of simulation
    :return:
    '''
    print(infoString(runtype=runtype, run=run, target=target, edge=edge, wc=wc, state=state), end='')

def read_ligands( target ):
    '''
    Reads in a dict of the ligands to be calculated for one target.
    :param target: target name
    :return: py:class`dict` with ligands
    '''
    try:
        return ligands.LigandSet(target).get_list()
    except Exception as e:
        print(e)
        return []

def read_edges( target ):
    '''
    Reads in a dict of the edges to be calculated for one target.
    :param target: target name
    :return: py:class`dict` with edges
    '''
    try:
        return edges.EdgeSet(target).get_dict()
    except Exception as e:
        print(e)
        return []

# a function to prepare a .top file
def create_top(fname='topol.top', ff='amber14sb_OL15.ff', water='tip3p',
               itp=['merged.itp'], mols=[['MOL', 1]],
               destination='', toppaths=[]):
    fp = open(fname, 'w')
    # ff itp
    fp.write('#include "%s/forcefield.itp"\n' % ff)
    # additional itp
    for i in itp:
        fp.write('#include "%s"\n' % i)
        # water itp
    fp.write('#include "%s/%s.itp"\n' % (ff, water))
    # ions
    fp.write('#include "%s/ions.itp"\n\n' % ff)
    # system
    fp.write('[ system ]\n')
    fp.write('simulation system\n\n')
    # molecules
    fp.write('[ molecules ]\n')
    for mol in mols:
        fp.write('%s %s\n' % (mol[0], mol[1]))
    fp.close()

    # also copy (if needed) topology files into the working folder
    if len(destination) > 0 and len(toppaths) > 0:
        for topfile in itp:
            for toppath in toppaths:
                if os.path.isfile(toppath + '/' + topfile):
                    try:
                        shutil.copy(toppath + '/' + topfile, destination)
                    except shutil.SameFileError:
                        pass
                    break


# assemble system structure
def assemble_system(system, molfiles=[]):
    fpout = open(system, 'w')
    for molfile in molfiles:
        fp = open(molfile, 'r')
        lines = fp.readlines()
        fp.close()

        for l in lines:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                fpout.write(l)
    fpout.close()

class pmxvariables:
    def __init__(self,
                 target='thrombin',
                 forcefield='smirnoff99Frosst.offxml',
                 path='./',
                 replicates=[1,2,3],
                 verbose=True):
        # path to pmx scripts
        self.scriptpath = pmx.__path__[0] + '/scripts/ligands/'
        # path to mdp
        self.mdpPath = os.path.abspath(pmx.__path__[0]) + '/data/mdppath/'
        # path to forcefield files
        self.ffPath = os.path.abspath(pmx.__path__[0]) + '/data/'
        # forcefield
        self.forcefield = forcefield
        # target
        self.target = target

        # path where to find edges, topologies, etc.
        self.basePath = f'{path}'

        targets.set_data_dir(path)
        for t in targets.target_dict:
           if t == self.target:
                self.targetDir = targets.target_dict[t]['dir']
                break
        else:
            self.targetDir = None
            print('Target not found')

        self.protPath = f'{self.basePath}/{self.targetDir}/01_protein/'
        self.ligPath = f'{self.basePath}/{self.targetDir}/02_ligands/'
        self.hybPath = f'{self.basePath}/{self.targetDir}/03_hybrid/'
        self.runPath = f'{self.basePath}/{self.targetDir}/06_pmx/'

        # all simulations will be done in 3 replicates
        self.replicates = replicates

        # make all output verbose
        self.verbose = verbose

    # # workpath/[water|protein]/edge* - every edge has its own folder
    # global waterProtein
    # waterProtein = ['water', 'protein']
    # # workpath/[water|protein]/edge*/state[A|B] - two states will be considered for every edge
    # global states
    # states = ['stateA', 'stateB']
    # # workpath/[water|protein]/edge*/state[A|B]/[em|eq_nvt|eq|morphes] - all simulations will be in separate folders
    # global simulations
    # simulations = ['em', 'eq_nvt', 'eq', 'morphes']

        self.ligands = read_ligands(target)

        # read pre-defined edges
        self.edges = read_edges(target)


def step1_parameterize_ligands(target, forcefield='openff-1.0.0', verbose=False):
    ligSet = ligands.LigandSet(target).get_dataframe(columns=['name', 'smiles', 'docked'])
    for index, lig in ligSet.iterrows():
        printInfo(runtype='para', run=0, target=target, edge=lig['name'].values[0], wc='', state='')
    return

def step2_create_hybrid_systems(target, forcefield='openff-1.0.0', repeats=3, verbose=False):
    edgesToUse = read_edges(target)
    for edge in edgesToUse.keys():
        lig1 = edgesToUse[edge][0]
        lig2 = edgesToUse[edge][1]
        for r in range(repeats):
            for wp in ['water', 'protein']:
                printInfo(runtype='hybrid', run=r, target=target, edge=edge, wc=wp, state='')
    return


def step3_simulation_setup(target='jnk1', forcefield='openff-1.0.0', simulationType=['em', 'nvt', 'eq', 'morphes'], queueType='sge', repeats=3, verbose=False):
    edgesToUse = read_edges(target)
    for edge in edgesToUse.keys():
        lig1 = edgesToUse[edge][0]
        lig2 = edgesToUse[edge][1]
        for r in range(repeats):
            for wp in ['water', 'protein']:
                for state in ['stateA', 'stateB']:
                    for sim in simulationType:
                        printInfo(runtype=sim, run=r, target=target, edge=edge, wc=wp, state=state)
    return

def step4_submit_simulations(target='jnk1', forcefield='openff-1.0.0', queueType='sge', repeats=3, verbose=False):
    edgesToUse = read_edges(target)
    for edge in edgesToUse.keys():
        lig1 = edgesToUse[edge][0]
        lig2 = edgesToUse[edge][1]
        for r in range(repeats):
            for wp in ['water', 'protein']:
                for state in ['stateA', 'stateB']:
                    printInfo(runtype='pmx', run=r, target=target, edge=edge, wc=wp, state=state)

    return

def step4_check_simulations(target='jnk1', forcefield='openff-1.0.0', queueType='sge', repeats=3, verbose=False):
    edgesToUse = read_edges(target)
    for edge in edgesToUse.keys():
        lig1 = edgesToUse[edge][0]
        lig2 = edgesToUse[edge][1]
        for r in range(repeats):
            for wp in ['water', 'protein']: 
                for state in ['stateA', 'stateB']:
                    printInfo(runtype='pmx', run=r, target=target, edge=edge, wc=wp, state=state)
                    print('running')
    return

def step5_analysis(target='jnk1', forcefield='openff-1.0.0', repeats=3, verbose=False):
    edgesToUse = read_edges(target)
    for edge in edgesToUse.keys():
        lig1 = edgesToUse[edge][0]
        lig2 = edgesToUse[edge][1]
        for r in range(repeats):
            printInfo(runtype='pmx', run=r, target=target, edge=edge, wc='', state='')
    return

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
