#!/usr/bin/env python
# coding: utf-8


# ./_cluster.py

# Description:
# TODO: Make better documenation


from __future__ import division
from subprocess import Popen, PIPE


def pwrmsd(pdbFile, ftFile, rotPRM, num):

"""
sblu measure pwrmsd [OPTIONS] PDB_FILE FTFILE ROTPRM
Options:
  --sort-ftresults          Sort ftresults before using
  -n, --nftresults INTEGER  Number of ftresults to use
  --only-CA                 Only C-alpha atoms
  --only-backbone           Only backbone atoms
  --only-interface
  --interface-radius FLOAT
  --rec PATH                PDB to use for calculating interface
  -o, --output FILENAME
  --help
"""

cmd = 'sblu measure pwrmsd -n {num} {pdbFile} {ftFile} {rotPRM}'
print(cmd)
proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
so, se = proc.communicate()


def clusterResults(pwrmsd):

"""
Usage: sblu docking cluster [OPTIONS] PWRMSDS

Options:
  -r, --radius FLOAT
  -s, --min-cluster-size INTEGER
  -l, --max-clusters INTEGER
  -o, --output FILENAME
  --json / --no-json
  --help                          Show this message and exit.
"""

cmd = 'sblu docking cluster {pwrmsd}'
print(cmd)
proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
so, se = proc.communicate()
