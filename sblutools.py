#!/usr/bin/env python
# coding: utf-8

# ./sblutools.py

# Description:
# Uses sblu (see: https://bitbucket.org/bu-structure/sb-lab-utils/src/master/)

from __future__ import division
from subprocess import Popen, PIPE


def rmsd(lig, crys_lig, ft, rot, output, options=None):
    """Usage: sblu measure ftrmsd [OPTIONS] LIG_FILE LIG_CRYS FTFILE ROTPRM

    Options:
      --sort-ftresults / --no-sort-ftresults
      -n, --nftresults INTEGER
      --only-CA                       Only C-alpha atoms
      --only-backbone                 Only backbone atoms
      --only-interface                Only use inteface atoms. Requires --rec.
      --interface-radius FLOAT        Radius around receptor to consider.
      --rec PATH                      Receptor file if using interface mode.
      -o, --output FILENAME           Write output to file (default: stdout)
      --help                          Show this message and exit.
    """
    DEFAULT_RMSD_OPTIONS = ""

    if options is None:
        options = DEFAULT_RMSD_OPTIONS

    cmd = 'sblu measure ftrmsd -o {} {} {} {} {} {}'.format(output, options, lig, crys_lig, ft, rot)
    print cmd
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()


def pwrmsd(pdb_file, ft, rot, output, options=None):
    """Usage: sblu measure pwrmsd [OPTIONS] PDB_FILE FTFILE ROTPRM

    Options:
      --sort-ftresults          Sort ftresults before using
      -n, --nftresults INTEGER  Number of ftresults to use
      --only-CA                 Only C-alpha atoms
      --only-backbone           Only backbone atoms
      --only-interface
      --interface-radius FLOAT
      --rec PATH                PDB to use for calculating interface
      -o, --output FILENAME
      --help                    Show this message and exit.
    """
    DEFAULT_PWRMSD_OPTIONS = ""

    if options is None:
        options = DEFAULT_PWRMSD_OPTIONS

    cmd = 'sblu measure pwrmsd -o {} {} {} {} {}'.format(output, options, pdb_file, ft, rot)
    print(cmd)
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()


def cluster(pwrmsd_file, output, options=None):
    """Usage: sblu docking cluster [OPTIONS] PWRMSDS

    Options:
      -r, --radius FLOAT
      -s, --min-cluster-size INTEGER
      -l, --max-clusters INTEGER
      -o, --output FILENAME
      --json / --no-json
      --help                          Show this message and exit.

    """

    DEFAULT_CLUSTER_OPTIONS = ""

    if options is None:
        options = DEFAULT_CLUSTER_OPTIONS

    cmd = 'sblu docking cluster -o {} {} {}'.format(output, options, pwrmsd_file)
    print(cmd)
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()
