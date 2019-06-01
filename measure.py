#!/usr/bin/env python
# coding: utf-8

# ./measure.py

# Description:
# Measures things
from __future__ import division
from prody import *
import math
import numpy as np
from itertools import islice
import linecache


def getVDW(atomName):
    """Maps an atom name to its corresponding VdW radii
    """
    vdwSwitch = {
        "H": 1.10,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "P": 1.80,
        "S": 1.75,
        "CL": 1.75,
    }
    return vdwSwitch.get(atomName)


def pairwise_dist(sel1, sel2, sel1_index=None, sel2_index=None):
    """Calculates a pairwise distance between a first and second selection within a defined interface.

    Arguments:
    - sel1 - (str) first selection (from User Input)
    - sel2 - (str) second selection (from User Input)
    """

    if sel1_index is not None:
        sel1_coords = sel1.getCoordsets(sel1_index)
    else:
        sel1_coords = sel1.getCoords()

    if sel2_index is not None:
        sel2_coords = sel2.getCoordsets(sel2_index)
    else:
        sel2_coords = sel2.getCoords()

    clashes = []
    distances = []
    for c1 in range(len(sel1_coords)):
        for c2 in range(len(sel2_coords)):
            pw_distance = math.sqrt(
                sum(map(lambda f: (f[0]-f[1])**2, zip(sel1_coords[c1], sel2_coords[c2]))))
            element_1 = sel1.getElements()[c1]
            element_2 = sel2.getElements()[c2]
            vdw_1 = getVDW(element_1)
            vdw_2 = getVDW(element_2)
            clash_dist = abs(vdw_1 + vdw_2)
            if clash_dist > pw_distance:
                pw_distance = round(pw_distance, 2)
                clash_dist = round(clash_dist, 2)
                clashes.append(clash_dist)
                distances.append(pw_distance)
    clashing_pw = zip(clashes, distances)
    return clashing_pw


def count_hydrophobic_res(selection_res_names):
    """selection_res_names needs to be in selection.getResnames() format.
    Example:
    lig = parsePDB(ligand)
    lig_interface = lig.select('all within 10 of another selection')
    lig_res_names = lig_interface.getResnames()
    """
    HYDROPHOBIC_RESIDUES = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'MET', 'TRP']
    count = 0
    for residue in selection_res_names:
        if residue in HYDROPHOBIC_RESIDUES:
            count += 1
    num_res = len(selection_res_names)
    num_hyd = count
    percent_hydrophobic = float(num_hyd)/float(num_res)
    return count, len(selection_res_names), round((percent_hydrophobic)*100, 2)
