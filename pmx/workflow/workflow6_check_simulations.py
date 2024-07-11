#!/usr/bin/env python
# coding: utf-8

"""workflow6_check_simulations.py: Checks state of pmx simulations."""

import os
import re
import subprocess
import argparse
from pmx.workflow import pmxworkflow
from pmx.workflow import workflow4_write_scripts
__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def checkSimulations(pwf, runtype, queueType='sge'):
    outputString=''

    # workpath/[water|complex]/edge* - every edge has its own folder
    waterComplex = ['water','complex']
    # workpath/[water|complex]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA','stateB']

    # queue output
    try:
        if queueType == 'slurm':
            process = subprocess.run(f'sacct --user dfhahn --format jobid,jobname%60,state'.split(),
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        elif queueType == 'sge':
            process = subprocess.run(f'qstat -u dhahn3 -xml -r'.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        queue = process.stdout.decode("utf-8")
        queueError = process.stderr.decode("utf-8")
    except Exception as e:
        queue = ''
        queueError = ''
        print('Queue could not be called.')

    # current working directory 
    cwd = os.getcwd()

    for edge in pwf.edges.keys():
        pmxworkflow.printInfo(runtype='Check Simulations', run='', target=pwf.target, edge=edge, wc='', state='')
        for wc in waterComplex:
            print(f'    - {wc}')
            for state in states:
                print(f'        - {state}')
                for run in pwf.replicates:
                    print(f'            - Replicate {run}')
                    runtypesForResub = []
                    for rt in ['em', 'nvt', 'eq', 'morphes']:
                        if rt not in runtype:
                            continue
                        print(f'                - {rt}: ', end='')
                        if rt == 'morphes':
                            jobname = f'morphes{run}_{pwf.target}_{wc}_{edge}_{state}'
                        else:
                            jobname = f'pmx{run}_{pwf.target}_{wc}_{edge}_{state}'
                        queuejob = ''
                        for line in queue.split('\n'):
                            if re.search(jobname, line) and (re.search('PENDING', line) or re.search('RUNNING', line)):
                                queuejob += line
                        if rt == 'morphes':
                            jobsresub = 0
                            jobsrun = 0
                            jobsfinished = 0
                            outputString += pmxworkflow.infoString(runtype=rt, run=run, target=pwf.target, edge=edge, wc=wc, state=state)
                            jobsforresub = []
                            for i in range(1, 81):
                                if re.search(f'[{i}]', queuejob) or re.search(f'_{i} ', queuejob):
                                    jobsrun += 1
#                                    outputString += f'\33[31mJob {i} in queue.\33[0m\n'
#                                    outputString += f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/\n'
                                    continue
                                if os.path.isfile(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/dgdl{i}.xvg'):
                                    count = len(open(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/dgdl{i}.xvg').readlines(  ))
                                    if count != 25019:
                                        print('\33[31mOutput file corrupted.\33[0m')
#                                        outputString += '\33[31mOutput file corrupted.\33[0m\n'
#                                        outputString += f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/\n'
                                        jobsforresub.append(i)
                                    else:
                                        jobsfinished += 1
                                else:
                                    #                               print(f'Simulation not yet done. File dgdl{i}.xvg does not exist.')                            
#                                    outputString += '\33[31mSimulation not (yet) finished.\33[0m\n'
                                    jobsforresub.append(i)
                      
                            if len(jobsforresub) > 0:
                                runtypesForResub.append(rt)
                      
                            print(f'Jobs in queue/failed/finished: {jobsrun}/{jobsresub}/{jobsfinished}')
                        else:
                            if os.path.isfile(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/{rt}{run}.log') and os.path.isfile(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/{rt}{run}.gro'):
                                print('log and final coordinate file exist:', end=' ')
                                if rt=='em':
                                    with open(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/{rt}{run}.log') as f:
                                        for line in f.readlines():
                                            if re.search('Steepest Descents converged', line):
                                                print(line.strip())
                                                break
                                        else:
                                            print('\33[31mSimulation not (yet) converged.\33[0m')
#                                            outputString += pmxworkflow.infoString(runtype=rt, run=run, target=pwf.target, edge=edge, wc=wc, state=state)
#                                            outputString += '\33[31mSimulation not (yet) converged.\33[0m\n'
#                                            outputString += queuejob + '\n'
                                            runtypesForResub.append(rt)
                                else:
                                    with open(f'{pwf.runPath}/{edge}/{wc}/{state}/{rt}{run}/{rt}{run}.log') as f:
                                        for line in f.readlines():
                                            if re.search('Finished mdrun', line):
                                                print(line.strip())
                                                break
                                        else:
                                            print('\33[31mSimulation not (yet) finished.\33[0m')
#                                            outputString += pmxworkflow.infoString(runtype=rt, run=run, target=pwf.target, edge=edge, wc=wc, state=state)
#                                            outputString += '\33[31mSimulation not (yet) finished.\33[0m'
#                                            outputString += queuejob + '\n'
                                            runtypesForResub.append(rt)

                            else:
                                print('\33[31mNo log file available, simulation probably not finished\33[0m')
#                                outputString += pmxworkflow.infoString(runtype=rt, run=run, target=pwf.target, edge=edge, wc=wc, state=state)
#                                outputString += '\33[31mNo log file available, simulation probably not finished\33[0m\n'
#                                outputString += queuejob + '\n'
                                runtypesForResub.append(rt)
                    if len(runtypesForResub) > 0:
                        print(f'Submitting job ... ')
                        exePath = f'{pwf.runPath}/{edge}/{wc}/{state}/'
                        workflow4_write_scripts.writeScript(pwf, edge, wc, state, runtypesForResub, run, queueType)
                        os.chdir(exePath)
                        if queueType == 'slurm':
                            process = subprocess.run(f'sbatch pmx{run}.sh'.split(),
                                                     stdout=subprocess.PIPE,
                                                     stderr=subprocess.PIPE)
                        elif queueType == 'sge':
                            print('qsub ... ')
                            process = subprocess.run(f'qsub pmx{run}.sh'.split(),
                                                     stdout=subprocess.PIPE,
                                                     stderr=subprocess.PIPE)
                            
                        if args.verbose:
                            print('STDOUT{} '.format(process.stdout.decode('utf8')))
                            print('STDERR{} '.format(process.stderr.decode('utf8')))
                                        
                        outputString += f'Submitted jobs\n'
                        os.chdir(cwd)

    return outputString


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', 
                        '--target', 
                        metavar = 'TARGET',
                        type=str,
                        default='jnk1',
                        help='The target protein.')
    parser.add_argument('-s',
                        '--simulationtype',
                        metavar = 'SIMULATION_TYPE',
                        nargs='+',
                        type=str,
                        default='em',
                        choices = ['em', 'eq', 'nvt', 'morphes'],
                        help='The simulation type.')
    parser.add_argument('-f', 
                        '--forcefield', 
                        metavar = 'FORCEFIELD',
                        type=str,
                        default='openff-2.0.0-rc.2.offxml',
                        choices = ['smirnoff99Frosst-1.1.0.offxml', 'openff-1.0.0.offxml', 'gaff2', 'openff-2.0.0-rc.2.offxml'],
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
    parser.add_argument('-q', 
                        '--queuetype', 
                        metavar = 'QUEUE',
                        type=str,
                        default='slurm',
                        choices = ['slurm', 'sge'],
                        help='The queue type of the job scripts.')
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

    outputString = checkSimulations(pwf, runtype=args.simulationtype, queueType=args.queuetype)

    
    print('=' * 80)
    print('=' * 80)
    print('=== SUMMARY OF UNSUCCESSFUL RUNS ===')
    print('=' * 80)
    print('=' * 80)
    print(outputString)



