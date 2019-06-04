#!/usr/bin/env python
# coding: utf-8

# ./measure.py

# Description:
# Measures things
from __future__ import division
from prody import *
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances


def vdw_radii(sel):
    """
    Returns a VdW radii for an element name (e.g., AtomGroup.getElements())
    :param sel: list of str(element_names)
    :return: numpy array of VdW radii
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

    radii = [vdwSwitch[i] for i in sel]
    radii = np.asarray(radii)

    return radii


def pairwise_distance(sel1, sel2):
    """
    Calculates the pairwise euclidean distances between two sets of atoms
    :param sel1: first atom group with dimensions [[N, [xn, yn, zn]]
    :param sel2: second atom group with dimensions [[M, [xm, ym, zm]]
    :return: numpy array pairwise distances with dimensions [N, M]
    """
    sel1_coords = sel1.getCoords()
    sel2_coords = sel2.getCoords()

    pw_distances = euclidean_distances(sel1_coords, sel2_coords)

    return np.asarray(pw_distances)


def vdw_sum(sel1, sel2):
    """
    Calculates the VdW pairwise summation for two sets of atoms
    :param sel1: first atom group, with each atom having a VdW radii
    :param sel2: second atom group, with each atom having a VdW radii
    :return: numpy array pairwise distances with dimensions [N, M]
    """

    sel1_elements = sel1.getElements()
    sel2_elements = sel2.getElements()

    sel1_radii = vdw_radii(sel1_elements)
    sel2_radii = vdw_radii(sel2_elements)

    sel1_arr = np.array(sel1_radii)
    sel2_arr = np.array(sel2_radii)

    vdw_distances = sel1_arr[:, None] + sel2_arr

    return np.asarray(vdw_distances)


def clash_count(sel1, sel2, threshold=1.0):
    """
    Calculates the number of clashes for two sets of atoms
    :param sel1: first atom group
    :param sel2: second atom group
    :param threshold: VdW softener
    :return: number of clash counts
    """
    pwd = pairwise_distance(sel1, sel2)
    vdw_d = vdw_sum(sel1, sel2)

    pwd = pwd.flatten()
    vdw_d = vdw_d.flatten()

    clash_counts = sum(vdw_d * threshold > pwd)

    return clash_counts

def distance(sel1, sel2):
    """
    Calculates the euclidean distance between two sets of atoms
    :param sel1: first atom group
    :param sel2: second atom group
    :return: euclidean distance
    """

    sel1_coords = sel1.getCoords()
    sel2_coords = sel2.getCoords()

    dist = np.linalg.norm(sel1_coords-sel2_coords, axis=1)

    return dist
