#!/usr/bin/env python
# coding: utf-8


# ./_cluster.py

# Description:
# TODO: Make better documenation


from __future__ import division
from subprocess import Popen, PIPE


def pwrmsd(pdbFile, ftFile, rotPRM, num, out, rec):

    cmd = 'sblu measure pwrmsd -o {} --only-interface --rec {} -n {} {} {} {}'.format(
        out, rec, num, pdbFile, ftFile, rotPRM)
    print(cmd)
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()


def clusterResults(pwrmsd, outter):

    cmd = 'sblu docking cluster -r 6.0 -o {} {}'.format(outter, pwrmsd)
    print(cmd)
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()
