#!/usr/bin/env python
# coding: utf-8

"""workflow7_analyse.py: Analyses the pmx non-equilibrium trajectories."""

import argparse
import os
import subprocess

import numpy as np
import pandas as pd
from pmx.workflow import pmxworkflow

from multiprocessing import Pool

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"


def read_neq_results( fname ):
    if not os.path.exists(fname):
        print('File does not exist')
        return []
    fp = open(fname,'r')
    lines = fp.readlines()
    fp.close()
    out = []
    for l in lines:
        l = l.rstrip()
        foo = l.split()
        if 'BAR: dG' in l:
            out.append(float(foo[-2]))
        elif 'BAR: Std Err (bootstrap)' in l:
            out.append(float(foo[-2]))
        elif 'BAR: Std Err (analytical)' in l:
            out.append(float(foo[-2]))
        elif 'BAR: Conv' in l:
            out.append(float(foo[-1]))
    return(out)


def output_textfile(pwf, df, fname):
    fp = open(fname, 'w')
    fp.write('#%6s  %6s  %6s\n' % ('1_edge', '2_ddg', '3_ddg_err'))
    fp.write('#all values in kcal/mol\n')

    for edge, ligs in pwf.edges.items():
        val = df.loc[f'{edge}_ddg', 'val']
        err = df.loc[f'{edge}_ddg', 'err']

        fp.write('%10s  %4.2f  %4.2f\n' % (edge, val, err))

    fp.close()

def run_commands(cmdline):
    process = subprocess.Popen(cmdline,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    
    out = process.communicate()
    print('STDOUT{} '.format(out[0].decode("utf-8")))
    print('STDERR{} '.format(out[1].decode("utf-8")))
    return

def analyzeSimulations(pwf):
    bootnum = 1000
    ######## read into a data frame #########
    df = pd.DataFrame()
    arrays = [['water', 'complex', 'ddg'] * 3,
              [int(i / 3) + 1 for i in range(9)],
              [''] * 9]
    arrays = [[] * 9] * 3
    tuples = list(zip(*arrays))

    index = pd.MultiIndex.from_tuples(tuples, names=['leg', 'repeat', ''])
    newdf = pd.DataFrame(columns=index)

    # workpath/[water|complex]/edge* - every edge has its own folder
    waterComplex = ['water','complex']
    # workpath/[water|complex]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA','stateB']

    commands = []
    for edge in pwf.edges.keys():
        pmxworkflow.printInfo(runtype='Analyse', run='', target=pwf.target, edge=edge, wc='', state='')
        for wc in waterComplex:
            print(f'    - {wc}')
            for run in pwf.replicates:
                os.makedirs(f'{pwf.runPath}/{edge}/{wc}/analyse{run}', exist_ok=True)

                # input files
                fa = [f'{pwf.runPath}/{edge}/{wc}/stateA/morphes{run}/dgdl{i}.xvg' for i in range(1, 81)]
                fb = [f'{pwf.runPath}/{edge}/{wc}/stateB/morphes{run}/dgdl{i}.xvg' for i in range(1, 81)]
                # output files
                result = f'{pwf.runPath}/{edge}/{wc}/analyse{run}/results.txt'
                outWorkA = f'{pwf.runPath}/{edge}/{wc}/analyse{run}/integA.dat'
                outWorkB = f'{pwf.runPath}/{edge}/{wc}/analyse{run}/integB.dat'
                # plot = f'{pwf.runPath}/{edge}/{wc}/analyse{run}/plot.png'

                cmdline = f'pmx analyse -fA {" ".join(fa)} '\
                f'-fB {" ".join(fb)} '\
                f'-o {result} '\
                f'-oA {outWorkA} '\
                f'-oB {outWorkB} '\
                f'--work_plot none '\
                f'-b 100 '\
                f'-t 298 '\
                f'-m bar '\
                f'--no_ks '\
                f'--unit kcal'
                
                commands.append(cmdline.split())
#                print(cmdline)
    with Pool(16) as p:
        p.map(run_commands, commands)

def unite_results(pwf):
    bootnum = 1000
    ######## read into a data frame #########
    df = pd.DataFrame()
    arrays = [['water', 'complex', 'ddg'] * 3,
              [int(i / 3) + 1 for i in range(9)],
              [''] * 9]
    arrays = [[] * 9] * 3
    tuples = list(zip(*arrays))

    index = pd.MultiIndex.from_tuples(tuples, names=['leg', 'repeat', ''])
    newdf = pd.DataFrame(columns=index)

    # workpath/[water|complex]/edge* - every edge has its own folder
    waterComplex = ['water','complex']
    # workpath/[water|complex]/edge*/state[A|B] - two states will be considered for every edge
    states = ['stateA','stateB']

    for edge in pwf.edges.keys():
        pmxworkflow.printInfo(runtype='Analyse', run='', target=pwf.target, edge=edge, wc='', state='')
        for wc in waterComplex:
            print(f'    - {wc}')
            for run in pwf.replicates:
                # output files
                result = f'{pwf.runPath}/{edge}/{wc}/analyse{run}/results.txt'
                if not os.path.exists(result):
                    print(f"Results file {result} does not exist.")

                # --------------------------------
                # Processing of the results: calculate ddG, compute means, output
                # ---------------------------------

                print(result)
                foo = read_neq_results( result )
                if len(foo) > 1:
                    df.loc[f'{edge}_{wc}{run}','val'] = foo[0]
                    df.loc[f'{edge}_{wc}{run}','err'] = foo[2]
                    df.loc[f'{edge}_{wc}{run}','aerr'] = foo[1]
                    df.loc[f'{edge}_{wc}{run}','conv'] = foo[3]
                    for t in ['val', 'err', 'aerr', 'conv']:
                        newdf.loc[f'{edge}', (wc, run, t)] = df.loc[f'{edge}_{wc}{run}', t]
                else:
                    df.loc[f'{edge}_{wc}{run}','val'] = np.nan
                    df.loc[f'{edge}_{wc}{run}','err'] = np.nan
                    df.loc[f'{edge}_{wc}{run}','aerr'] = np.nan
                    df.loc[f'{edge}_{wc}{run}','conv'] = np.nan
                    for t in ['val', 'err', 'aerr', 'conv']:
                        newdf.loc[f'{edge}', (wc, run, t)] = df.loc[f'{edge}_{wc}{run}', t]
                    print('Results could not be read')


        vals = []
        errs = []
        aerrs = []
        for run in pwf.replicates:
            ##### calculate ddg #####
            ddg = df.loc[f'{edge}_complex{run}','val'] - df.loc[f'{edge}_water{run}','val']
            err = np.sqrt( np.power(df.loc[f'{edge}_complex{run}','err'],2.0) +
                           np.power(df.loc[f'{edge}_water{run}','err'],2.0) )
            aerr = np.sqrt( np.power(df.loc[f'{edge}_complex{run}','aerr'],2.0) +
                            np.power(df.loc[f'{edge}_water{run}','aerr'],2.0) )
            print(ddg, err, aerr)
            df.loc[f'{edge}_ddg{run}','val'] = ddg
            df.loc[f'{edge}_ddg{run}','err'] = err
            df.loc[f'{edge}_ddg{run}','aerr'] = aerr
            newdf.loc[f'{edge}', ('ddg', run, 'val')] = ddg
            newdf.loc[f'{edge}', ('ddg', run, 'err')] = err
            newdf.loc[f'{edge}', ('ddg', run, 'aerr')] = aerr
            vals.append(ddg)
            errs.append(err)
            aerrs.append(aerr)


        ###### calculate mean dg with err ######
        # mean ddg
        mean = np.average(vals)
        df.loc[f'{edge}_ddg', 'val'] = mean

        # error
        # 1) create three distributions
        distribs = []
        for v, e in zip(vals, errs):
            distribs.append(np.random.normal(v, e, size=bootnum))
        if len(distribs) > 1:
            distr = np.vstack(distribs)
            # 2) calculate stderrs
            stderr = np.mean(np.sqrt(np.var(distr, ddof=1, axis=0) / np.float(len(distribs))))
            df.loc[f'{edge}_ddg', 'err'] = stderr
            print(mean, stderr)
        else:
            stderr = errs[0]
            df.loc[f'{edge}_ddg','err'] = stderr

        for t in ['val', 'err', 'aerr']:
            newdf.loc[f'{edge}', ('ddg_mean', '-', t)] = df.loc[f'{edge}_ddg', t]

    ###### output ######
    os.makedirs(f'{pwf.runPath}/{pwf.forcefield}/results', exist_ok=True)
    output_textfile( pwf,  df, f'{pwf.runPath}/{pwf.forcefield}/results/{pwf.target}_{pwf.forcefield}.dat' )

    newdf.to_csv(f'{pwf.runPath}/{pwf.forcefield}/results/{pwf.target}_{pwf.forcefield}.csv', float_format='%10.2f')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',
                        '--target',
                        metavar='TARGET',
                        type=str,
                        default='jnk1',
                        help='The target protein.')
    parser.add_argument('-s',
                        '--simulationtype',
                        metavar='SIMULATION_TYPE',
                        nargs='+',
                        type=str,
                        default=['em'],
                        choices=['em', 'eq', 'nvt', 'morphes'],
                        help='The simulation type.')
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
    parser.add_argument('-c',
                        '--combine',
                        type=bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='Just combine data without calculation of work values.')
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

    if not args.combine:
        analyzeSimulations(pwf)
        
    unite_results(pwf)
