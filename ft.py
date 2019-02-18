#!/usr/bin/env python
# coding: utf-8

# ./ft.py

# Description:
# Returns a ft file based on a type of input


from __future__ import division
import os
import sys

from glob import glob
from subprocess import Popen, PIPE

import pandas as pd
from pymol import cmd


def fromList(ftList, ftOrig, ftPathOut):
    """A generalized script that generates an FT file from a list
    """
    for entry in ftList:
        with open(ftOrig, 'r') as rf:
            for num, ff in enumerate(rf):
                if num == int(entry):
                    with open(ftPathOut, 'a') as ftOut:
                        ftOut.write(ff)


def fromPymolSession(session, ftOrig, ftPathOut):
    """Writes one or more FT files based on ligand names that are saved to a
    given PyMOL session.

    Assumes that the ligand is saved in the form of "*.ft_entry.ft_file.00.pdb"

    """

    if os.path.isfile(session) is True:
        cmd.load(session)
        pymol.finish_launching()
        # otherwise stops @ error 'PyMOL not running, entering library mode (experimental)'
        cmd.load(session)

        ligands = cmd.get_object_list("*.00")

        lig_000 = []
        lig_002 = []
        lig_004 = []
        lig_006 = []

        for lig in ligands:
            ft_ligand = lig.split('.')[-2]
            name_ligand = lig.split('.')[-3]

            if ft_ligand == '000':
                lig_000.append(name_ligand)

            elif ft_ligand == '002':
                lig_002.append(name_ligand)

            elif ft_ligand == '004':
                lig_004.append(name_ligand)

            else:
                lig_006.append(name_ligand)

        if len(lig_000) > 0:
            ftPathOut = os.path.join(ftPathOut, "000")
            fromList(lig_000, ftOrig, ftPathOut)

        elif len(lig_002) > 0:
            ftPathOut = os.path.join(ftPathOut, "002")
            fromList(lig_002, ftOrig, ftPathOut)

        elif len(lig_004) > 0:
            ftPathOut = os.path.join(ftPathOut, "004")
            fromList(lig_004, ftOrig, ftPathOut)

        elif len(lig_006) > 0:
            ftPathOut = os.path.join(ftPathOut, "006")
            fromList(lig_006, ftOrig, ftPathOut)

        else:
            print("Check to see if the ligands in the PyMOL session are formatted as: *.ft_entry.ft_file.00.pdb ")


def fromReport(report, ftOrig, ftPathOut, cluster=None):
    """Writes a FT file based on a RMSD Report generated from ./rmsd.reportRMSD
    """
    entries = []
    if cluster is None:
        with open(report, 'r') as rf:
            lines = rf.readlines()[1:]
            for line in lines:
                bn = line.split('\t')[0]
                entries.append(bn)
        fromList(entries, ftOrig, ftPathOut)
