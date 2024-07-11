#!/usr/bin/env python
# coding: utf-8

"""workflow5_submit_simulations.py: Submits simulation scripts."""

import os
import subprocess
import argparse
from pmx.workflow import pmxworkflow

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def submitSimulations(pwf, queueType='sge', numsim=None):
    numstarted = 1
    # workpath/[water|complex]/edge* - every edge has its own folder
    waterComplex = ['water','complex']
    # workpath/edge*/[water|complex]/state[A|B] - two states will be considered for every edge
    states = ['stateA','stateB']

    cwd = os.getcwd()

    for edge in pwf.edges.keys():
        pmxworkflow.printInfo(runtype='Prepare Sim Scripts', run='', target=pwf.target, edge=edge, wc='', state='')
        for wc in waterComplex:
            print(f'    - {wc}')
            for state in states:
                print(f'        - {state}')
                for run in pwf.replicates:
                    print(f'            - Replicate {run}')
                    jobname = f'pmx{run}_{pwf.target}_{wc}_{edge}_{state}'
                    exePath = f'{pwf.runPath}/{edge}/{wc}/{state}/'
                
                    os.chdir(exePath)
                
                    if os.path.isfile(f'eq{run}/eq{run}.log'):
                        print('Simulation has been run already.\nContinue with next simulation')
                        os.chdir(cwd)
                        continue

                    if not os.path.isfile(f'pmx{run}.sh'):
                        print('Submission script does not exist.\nContinue with next simulation.')
                        os.chdir(cwd)
                        continue

                    if queueType == 'sge':
                        process = subprocess.run(f'qsub pmx{run}.sh'.split(), 
                                                 stdout=subprocess.PIPE, 
                                                 stderr=subprocess.PIPE)
                    elif queueType == 'slurm':
                        process = subprocess.run(f'sbatch pmx{run}.sh'.split(), 
                                                 stdout=subprocess.PIPE, 
                                                 stderr=subprocess.PIPE)

                    
                    if args.verbose:
                        print('STDOUT{} '.format(process.stdout.decode('utf8')))
                        print('STDERR{} '.format(process.stderr.decode('utf8')))
                    
                    os.chdir(cwd)

                    numstarted += 1
                    if numsim != None and numstarted > numsim:
                        break
                    
                else:
                    continue
                break
            else:
                continue
            break
        else:
            continue
        break

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', 
                        '--target', 
                        metavar = 'TARGET',
                        type=str,
                        default='01_jnk1',
                        help='The target protein.')
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
    parser.add_argument('-n', 
                        '--numsim', 
                        metavar = 'NUMSIM',
                        type=int,
                        default=None,
                        help='Number of simulations to be started. Default results in normal behavior.')
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

    submitSimulations(pwf, queueType=args.queuetype, numsim=args.numsim)

                

    
