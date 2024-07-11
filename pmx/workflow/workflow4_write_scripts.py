#!/usr/bin/env python
# coding: utf-8

"""workflow4_write_scripts.py: Writes simulation scripts."""

import argparse
import os
import sys
import shutil
import glob
from pmx.workflow import pmxworkflow, sge, slurm

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def decide_on_resources( wc, simType ):
    # simulation time in days
    simtime = 1
    # number of cpus used
    simcpu = 4
    if simType=='eq':
        simtime = 2
        simcpu = 8
        if wc=='complex':
            simtime = 7
            simcpu = 8
    return(simtime,simcpu)


def query_yes_no(question, default='no'):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes": "yes", "y": "yes",
             "no": "no", "n": "no"}
    if default == None:
        prompt = " [yes/no] "
    elif default == "yes":
        prompt = " [yes/no] "
    elif default == "no":
        prompt = " [yes/no] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " \
                             "(or 'y' or 'n').\n")

def getRunCoord(pwf, runtype, run=1, target='xxx', edge='yyy_zzz', wc='water', state='stateA'):
    if runtype == 'em':
        return f'{pwf.hybPath}/{edge}/{wc}/crd/ions{run}.pdb'
    elif runtype == 'nvt':
        return f'{pwf.runPath}/{edge}/{wc}/{state}/em{run}/em{run}.gro'
    elif runtype == 'eq':
        return f'{pwf.runPath}/{edge}/{wc}/{state}/nvt{run}/nvt{run}.gro'
    elif runtype == 'morphes':
        return f'{pwf.runPath}/{edge}/{wc}/{state}/eq{run}/eq{run}.gro'
    else:
        print('runtype not known')
        return ''


def writeScript(pwf, edge, wc, state, runtype, run, queueType):
    jobname = f'pmx{run}_{pwf.target}_{edge}_{wc}_{state}'
    # make parent directory for runs
    os.makedirs(os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}'), exist_ok=True)
    jobscriptFile = f'{pwf.runPath}/{edge}/{wc}/{state}/pmx{run}.sh'
    if queueType == 'sge':
        sge.scriptHeader(jobscriptFile, jobname)
    elif queueType == 'slurm':
        slurm.scriptHeader(jobscriptFile, jobname, simtime=7, simcpu=8)

    for rt in ['em', 'nvt', 'eq', 'morphes']:
        if rt not in runtype:
            continue
        # set simulation path
        simPath = os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/')
        # create folder
        os.makedirs(simPath, exist_ok=True)

        if os.path.isfile(f'{pwf.runPath}/{edge}/{wc}/{state}/pmx{run}.log'):
            # run log exists, so nothing is deleted
            print(f'Simulation {rt} has been run already, nothing is deleted. Continue with next simulation.')
            continue

        # specify input files
        mdp = os.path.abspath(f'{pwf.mdpPath}/{rt}_{state}.mdp')
        topology = os.path.abspath(f'{pwf.hybPath}/{edge}/{wc}/top/{pwf.forcefield}/topol{run}.top')
        coord = os.path.abspath(getRunCoord(pwf, rt, run=run, target=pwf.target, edge=edge, wc=wc, state=state))
        # specify output files
        tprfile = f'{rt}{run}.tpr'  # temporary tpr file
        mdout = f'mdout.mdp'

        # decide about resources
        simtime, simcpu = decide_on_resources(wc, rt)

        if rt == 'morphes':
            tprfile = os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}/eq{run}/eq{run}.tpr') # tpr file  
            trjfile = os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}/eq{run}/eq{run}.trr') # trr file  
            mdout = os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}/mdout.mdp')

            framefile = os.path.abspath(f'{simPath}/frame.gro') # frames

            commands = f'echo 0 | '\
                       f'gmx trjconv -s {tprfile}\\\n'\
                       f'            -f {trjfile}\\\n'\
                       f'            -o {framefile}\\\n'\
                       f'            -sep \\\n'\
                       f'            -ur compact \\\n'\
                       f'            -pbc mol \\\n'\
                       f'            -b 2256\n\n'\
                       f'mv {simPath}/frame0.gro {simPath}/frame80.gro\n\n'

            if queueType == 'sge':
                sge.scriptAppend(jobscriptFile, commands)
            elif queueType == 'slurm':
                slurm.scriptAppend(jobscriptFile, commands)

            # create array jobscript
            arrayjobname = f'morphes{run}_{pwf.target}_{edge}_{wc}_{state}'
            arrayjobscriptFile = f'{pwf.runPath}/{edge}/{wc}/{state}/morphes{run}.sh'

            if queueType == 'sge':
                gromppline = f'gmx grompp -p  {topology}\\\n'\
                             f'           -c  {simPath}/frame$SGE_TASK_ID.gro\\\n'\
                             f'           -o  tpr.tpr\\\n'\
                             f'           -f  {mdp}\\\n'\
                             f'           -po mdout$SGE_TASK_ID.mdp\\\n'\
                             f'           -maxwarn 2\n'

                sge.scriptArrayjob(arrayjobscriptFile, gromppline, simPath, arrayjobname, rt, run)
                submit_exec = shutil.which('qsub')
            elif queueType == 'slurm':
                gromppline = f'gmx grompp -p  {topology}\\\n'\
                             f'           -c  {simPath}/frame$SLURM_ARRAY_TASK_ID.gro\\\n'\
                             f'           -o  tpr.tpr\\\n'\
                             f'           -f  {mdp}\\\n'\
                             f'           -po mdout$SLURM_ARRAY_TASK_ID.mdp\\\n'\
                             f'           -maxwarn 2\n'

                slurm.scriptArrayjob(arrayjobscriptFile, gromppline, simPath, arrayjobname, rt, run)
                submit_exec = shutil.which('sbatch')

            submitPath = os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}/')
            commands = f'cd {submitPath}\n'\
                       f'{submit_exec} morphes{run}.sh\n\n'
            if queueType == 'sge':
                sge.scriptAppend(jobscriptFile, commands)
            elif queueType == 'slurm':
                slurm.scriptAppend(jobscriptFile, commands)

        else: # em, nvt, eq
            gromppline = f'gmx grompp -p {topology} '\
                         f'-c {coord} '\
                         f'-o {tprfile} '\
                         f'-f {mdp} '\
                         f'-po {mdout} '\
                         f'-maxwarn 3'

            if queueType == 'sge':
                sge.scriptMain(jobscriptFile, gromppline, simPath, rt, run, simcpu=simcpu)
            elif queueType == 'slurm':
                slurm.scriptMain(jobscriptFile, gromppline, simPath, rt, run, simcpu=simcpu)
    if queueType == 'sge':
        sge.scriptFooter(jobscriptFile)
    elif queueType == 'slurm':
        slurm.scriptFooter(jobscriptFile)

def deleteRunFiles(pwf, runtype):
    if query_yes_no('WARNING: All files in simulation directories will be deleted '
                    '(specified by -d/--deleterunfiles).'
                    'Do you really want to continue?',
                    default='no') == 'no':
        exit()
    # workpath/[water|complex]/edge* - every edge has its own folder
    waterComplex = ['water', 'complex']
    # workpath/[water|complex]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA', 'stateB']

    for edge in pwf.edges.keys():
        pmxworkflow.printInfo(runtype='Delete Sim files', run='', target=pwf.target, edge=edge, wc='', state='')
        for wc in waterComplex:
            for state in states:
                for run in pwf.replicates:
                    for rt in ['em', 'nvt', 'eq', 'morphes']:
                        if rt not in runtype:
                            continue
                        # set simulation path
                        simPath = os.path.abspath(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/')
                        # create folder
                        os.makedirs(simPath, exist_ok=True)
                        # everything is deleted!
                        toclean = glob.glob(f'{simPath}/*.*')
                        for clean in toclean:
                            os.remove(clean)

def prepareSimulations(pwf, runtype, queueType):

    # workpath/edge*/[water|complex]/ - every edge has its own folder
    waterComplex = ['water', 'complex']
    # workpath/[water|complex]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA', 'stateB']

    for edge in pwf.edges.keys():
        pmxworkflow.printInfo(runtype='Prepare Sim Scripts', run='', target=pwf.target, edge=edge, wc='', state='')
        for wc in waterComplex:
            print(f'    - {wc}')
            for state in states:
                print(f'        - {state}')
                for run in pwf.replicates:
                    print(f'            - Replicate {run}')
                    writeScript(pwf, edge, wc, state, runtype, run, queueType)

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
    parser.add_argument('-s',
                        '--simulationtype',
                        metavar = 'SIMULATION_TYPE',
                        nargs='+',
                        type=str,
                        default='em',
                        choices = ['em', 'eq', 'nvt', 'morphes'],
                        help='The simulation type.')
    parser.add_argument('-q',
                        '--queuetype',
                        metavar = 'QUEUE',
                        type=str,
                        default='slurm',
                        choices = ['slurm', 'sge'],
                        help='The queue type of the job scripts.')
    parser.add_argument('-d',
                        '--deleterunfiles',
                        type=bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Deletes run files before writing new simulation scripts.')
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
                                   replicates=replicatesToUse,
                                   verbose=args.verbose)

    if args.edges[0] != 'all':
        edges = dict()
        for edge in args.edges:
            edges[edge] = pwf.edges[edge]
        pwf.edges = edges

    if args.deleterunfiles:
        deleteRunFiles(pwf, args.simulationtype)


    prepareSimulations(pwf, args.simulationtype, args.queuetype)
