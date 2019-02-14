#!/usr/bin/env python
# coding: utf-8

# ./rmsd.py

# Description:
# Calculates an rmsd for based on a ligand 
# and a crystallized ligand


from __future__ import division
from subprocess import Popen,PIPE


def reportRMSD(ligfile, ligcrys, ftfile, rotprm, output):
    """Uses the sblu call:
    
    Usage: sblu measure ftrmsd [OPTIONS] LIG_FILE LIG_CRYS FTFILE ROTPRM

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
    ftname = os.path.basename(ftfile)
    cmd = 'sblu measure ftrmsd -o {o} {lf} {lc} {ft} {rp}'.format(o=output, lf=ligfile, 
                                                                  lc=ligcrys, ft=ftfile, rp=rotprm)
    print cmd
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so,se = proc.communicate()



def graphRMSD():
    #TODO: - Need to fix matplotlib depedency issues...

