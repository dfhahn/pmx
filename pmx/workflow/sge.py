#!/usr/bin/env python
# coding: utf-8

"""sge.py: Handles SGE queue system commands."""

import os
import pmx

__author__ = "David Hahn and Vytas Gapsys"
__copyright__ = "Copyright (c) 2020 Open Force Field Consortium and de Groot Lab"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "David Hahn"
__email__ = "davidfriedrichhahn@gmail.com"
__status__ = "Development"

def scriptHeader( fname, jobname, simtime=1, simcpu=4):
    fp = open(fname,'w')
    gmxlibpath = os.path.abspath(pmx.__path__[0]) + 'data/mutff45/'
    jobline = f'#! /usr/bin/bash\n\
#$ -N {jobname}\n\
#$ -cwd\n\
#$ -q all.q\n\
#$ -j yes\n\
#$ -l slot_type=gromacs\n\n\
export LD_LIBRARY_PATH=/shared/app/cuda/10.1/toolkit/lib64/:$LD_LIBRARY_PATH\n\n\
source /usr/local/gromacs/2021/bin/GMXRC\n\n\
export GMXLIB={gmxlibpath}\n\n\
mkdir -p $TMPDIR/$JOB_ID/\n\
cd $TMPDIR/$JOB_ID/\n\n'
    fp.write(jobline)
    fp.close()


def scriptAppend( fname, commands ):
    fp = open(fname, 'a')
    fp.write(commands)
    fp.close()


def scriptMain( fname, gromppline, simpath, runtype, run, simtime=4, simcpu=1 ):
    fp = open(fname, 'a')

    runargs = ''
    print(runtype)
    if runtype == 'em':
        runargs = f'-ntmpi 1 -ntomp {simcpu} -pin on'
    else:
        runargs = f'-pme gpu -ntmpi 1 -ntomp {simcpu} -pin on'

    jobline = f'{gromppline}\n\
gmx mdrun {runargs} -s {runtype}{run}.tpr -deffnm {runtype}{run}\n\
rsync -avzIi ./ {simpath}/\n\
rm -v $TMPDIR/$JOB_ID/*.*\n\n'
    fp.write(jobline)
    fp.close()


def scriptFooter( fname ):
    fp = open(fname, 'a')

    jobline = f'\n\
rmdir -v $TMPDIR/$JOB_ID\n\
            '
    fp.write(jobline)
    fp.close()


def scriptArrayjob(fname, gromppline, simpath, jobname, runtype, run, simtime=4, simcpu=1):
    fp = open(fname,'w')
    
    runargs = f'-pme gpu -ntmpi 1 -ntomp {simcpu} -pin on'
    gmxlibpath = os.path.abspath(pmx.__path__[0]) + '/data/mutff45/'
    jobline = f'#! /usr/bin/bash\n\
#$ -N {jobname}\n\
#$ -cwd\n\
#$ -q all.q\n\
#$ -j yes\n\
#$ -t 1-80\n\
#$ -l slot_type=gromacs\n\n\
export LD_LIBRARY_PATH=/shared/app/cuda/10.1/toolkit/lib64/:$LD_LIBRARY_PATH\n\n\
source /usr/local/gromacs/2021/bin/GMXRC\n\n\
export GMXLIB={gmxlibpath}\n\n\
\n\
WORKDIR=$(pwd)\n\
mkdir -p $TMPDIR/$JOB_ID/$SGE_TASK_ID\n\
cd $TMPDIR/$JOB_ID/$SGE_TASK_ID\n\
\n\
{gromppline}\n\n\
gmx mdrun {runargs} -s tpr.tpr \n\n\
rsync -avzIi dhdl.xvg {simpath}/dgdl$SGE_TASK_ID.xvg\n\n\
rm -v $TMPDIR/$JOB_ID/$SGE_TASK_ID/*.*\n\n\
rmdir -v $TMPDIR/$JOB_ID/$SGE_TASK_ID'
    fp.write(jobline)    
    fp.close()

