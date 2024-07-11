#!/usr/bin/env python
# coding: utf-8

"""workflow1_analyse.py: Parameterizes ligand molecules."""

import argparse
import os
import numpy as np
import warnings

import parmed as pmd

from pmx.workflow import pmxworkflow
from pmx import forcefield as pmxff

from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField

from multiprocessing import Pool

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

# function to convert openFF molecule (ligand) to a parmed structure
def ligandToPMD(ligand, pwf):
    ligand_positions = ligand.conformers[0]

    # Calculate am1bcc charges
    print(ligand.partial_charges)
    if ligand.partial_charges == None:
        try:
            ligand.compute_partial_charges_am1bcc()
        except Exception as e:
            raise Exception('Error in charge calculation for ligand {}: {}'.format(ligand.name, e))
    # Give all atoms unique names so we can export to GROMACS
    for idx, atom in enumerate(ligand.atoms):
        atom.name = f'{atom.element.symbol}{idx}'

    # initialize ForceField object
    openff_forcefield = ForceField(pwf.forcefield)
    # Do not assign H-bond constraints now, instead have ParmEd add them later
    del openff_forcefield._parameter_handlers['Constraints']

    ligand_topology = ligand.to_topology()
    try:
        ligand_system = openff_forcefield.create_openmm_system(ligand_topology, charge_from_molecules=[ligand])
    except Exception as e:
        raise Exception('Error in creating openmm system: {}'.format(e))
    # Create OpenMM Topology from OpenFF Topology
    omm_top = ligand_topology.to_openmm()

    # Convert OpenMM System to a ParmEd structure.
    pmd_structure = pmd.openmm.load_topology(omm_top, ligand_system, ligand_positions)

    return pmd_structure, ligand_topology, ligand_system, ligand_positions


# functions for Gromacs force field  manuplation/conversion
def set_charge_to_zero(itp):
    q = 0.0
    n = 0
    for a in itp.atoms:
        q += a.q
        n += 1
    intq = round(q)
    diffq = intq - q
    # round to 6 digit precision
    deltaq = np.around(diffq / float(n), decimals=6)

    newq = 0.0
    for a in itp.atoms:
        a.q += deltaq
        a.q = np.around(a.q, decimals=6)
        newq += a.q
    # add remainder to first atom
    intq = round(newq)
    diffq = intq - newq
    itp.atoms[0].q += diffq


def change_atomtypes(itp, suffix):
    for a in itp.atoms:
        #        newtype = str(a.atomtype)+'_'+str(a.name)+str(suffix)
        newtype = str(a.atomtype) + str(suffix)
        a.atomtype = newtype

    newdict = {}
    for atkey in itp.atomtypes.keys():
        at = itp.atomtypes[atkey]
        newtype = str(atkey) + str(suffix)
        newdict[newtype] = at
    itp.atomtypes = newdict


def write_ff(atypes, fname, ff='amber14sb'):
    fp = open(fname, 'w')
    print('[ atomtypes ]', file=fp)
    for atkey in atypes.keys():
        at = atypes[atkey]
        print('%8s %12.6f %12.6f %3s %12.6f %12.6f' % (atkey, at[1], at[2], at[3], at[4], at[5]), file=fp)


def write_posre(itp, fname, fc=1000):
    fp = open(fname, 'w')
    print('[ position_restraints ]', file=fp)
    for i, atom in enumerate(itp.atoms):
        print("%d   1    %d   %d    %d" % (i + 1, fc, fc, fc), file=fp)
    fp.close()

def calculate_parameters(lig, pwf):
    print(f'    - {lig}')

    # make topology directory
    os.makedirs(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/', exist_ok=True)
    
    # ligand sdf file
    ligFile = f'{pwf.ligPath}/{lig}/crd/{lig}.sdf'

    if os.path.isfile(f'{ligFile}'):
        ligand = Molecule.from_file(f'{ligFile}', allow_undefined_stereo=True)
    elif os.path.isfile(f'{pwf.ligPath}/{lig}/crd/{lig}.pdb'):
        # Try to read in PDB file instead of a SDF, only works with OpenEye
        warnings.warn('    SDF file not available. Trying to read in PDB file and automatically convert it to SDF. '
                      'This might lead to wrong bond orders.')
        ligand = Molecule.from_file(
            f'{pwf.ligPath}/{lig}/crd/{lig}.pdb',
            allow_undefined_stereo=True)
        # save as sdf file
        # ATTENTION: automatic conversion to SDF
        ligand.name = lig
        ligand.to_file(f'{ligFile}', 'sdf')
    else:
        warnings.warn(f'      File not found. Ligand {lig} cannot be read in. Continuing with next ligand.')
        return

    try:
        pmd_structure, ligand_topology, ligand_system, ligand_positions = ligandToPMD(ligand, pwf)
    except Exception as e:
        print('      ' + str(e))
        return

    # Export GROMACS files.
    pmd_structure.save(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/{lig}.top', overwrite=True)
    pmd_structure.save(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/{lig}.gro', overwrite=True, precision=3)

    # Create GROMACS ITP file
    itp = pmxff.read_gaff_top(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/{lig}.top')
    itp.set_name('MOL')
    change_atomtypes(itp, lig)
    set_charge_to_zero(itp)
    
    itp.write(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/{lig}.itp')
    write_ff(itp.atomtypes, f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/ff{lig}.itp')
    write_posre(itp, f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/posre_{lig}.itp')

    # Export AMBER files.
    pmd_structure.save(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/{lig}.prmtop', overwrite=True)
    pmd_structure.save(f'{pwf.ligPath}/{lig}/top/{pwf.forcefield}/{lig}.inpcrd', overwrite=True)

def coordToTopo(pwf):
    pmxworkflow.printInfo(runtype='Parameterize Ligands',
                          run='',
                          target=pwf.target,
                          edge='', wc='', state='')

    with Pool(4) as p:
        p.starmap(calculate_parameters, [(lig, pwf) for lig in pwf.ligands])


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
    parser.add_argument('-l',
                        '--ligands',
                        metavar='LIGANDS',
                        nargs='+',
                        type=str,
                        default=['all'],
                        help='Either "all" or a list of ligands to be calculated.')
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

    if args.ligands[0] != 'all':
        pwf.ligands = args.ligands

    coordToTopo(pwf)
