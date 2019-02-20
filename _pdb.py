
# !/usr/bin/env python
# coding: utf-8


# ./_pdb.py

# Description:
# TODO: Make better documenation

from __future__ import division
import numpy as np
from pymol import cmd, stored, math
from subprocess import Popen, PIPE
from prody import *
import math
from itertools import islice
import os


FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3))])


def readRotations(filePath, limit=None):
    """
    Reads 3x3 rotation matrices from a file.

    Rotations may be a text file with 9 or 10 columns. If the text file has 10
    columns, the first column is assumed to be the line number, which will be discarded.

    Returns a numpy array with dimensions N x (3x3) (an array of 3x3 rotation
    matrices), where N is the smaller of number of rotations in the file and limit,
    if limit is provided.
    """
    with open(filePath, 'r') as stream:
        return readRotationsStream(stream, limit)


def readRotationsStream(stream, limit=None):
    """
    Read rotations from a stream.

    Rotations may be a text file with 9 or 10 columns. If the text file has 10
    columns, the first column is assumed to be the line number, which will be discarded.

    Returns a numpy array with dimensions N x (3x3) (an array of 3x3 rotation
    matrices), where N is the smaller of number of rotations in the file and limit,
    if limit is provided.
    """
    rotations = np.loadtxt(
        islice(iter(stream), 0, limit))
    if rotations.shape[-1] == 10:
        rotations = rotations[:, 1:]
    return rotations.reshape(-1, 3, 3)


def readFTResults(filePath, limit=None):
    """
    Reads ftResults from a file.

    See readFTResultsStream for details.
    """
    with open(filePath, "r") as f:
        return readFTResultsStream(f, limit)


def readFTResultsStream(stream, limit=None):
    """
    Read ftResults from a stream.

    ftResults are assumed to be in a text file with at least 5
    columns.  The first column will be the rotation index. The next
    three columns are the translation vector, and the last column is
    the total weighted energy.
    """
    stream = iter(stream)

    return np.loadtxt(
        islice(stream, 0, limit),
        dtype=FTRESULT_DTYPE,
        usecols=(0, 1, 2, 3, 4))


def getCenterandTV(atomGroup, ftResults, rotations, center=None):
    orig_coords = atomGroup.getCoords()
    if center is None:
        center = np.mean(orig_coords, axis=0)
    postTV = np.expand_dims(ftResults['tv'] + center, 1)
    center = center * -1
    return center, postTV


def getTTT(ft, rot, center, postTV):
    """
    Arguments include:
    ft-output fromread_ftresults_stream
    rot-output fromread_rotations_stream
    center-output fromgetCenterandTV
    postTV-output fromgetCeneterandTV

    Will return a PyMOL happy matrix (4x4) with the top-left most
    (3x3) is a rotation matrix, the right (3x1) is a post-rotation
    translation vector (accounting for the center) and the bottom
    (1x4) is the pre-rotational matrix ([-1]*center)
    """

    rot_number = ft['roti']
    rot_matrix = rot[rot_number]
    pymolMatrix = np.append(rot_matrix, postTV, axis=1)
    zeros = np.append(center, '1')
    pymolMatrix = np.append(pymolMatrix, zeros)
    return pymolMatrix


def transformInputs(rotationStream, entry, masterLig, ligName):
    master_lig_pdb = parsePDB(masterLig)
    centers, postTV = getCenterandTV(master_lig_pdb, entry, rotationStream)
    TTT = getTTT(entry, rotationStream, centers, postTV)
    target = cmd.get_object_list("{}".format(ligName))
    return target, TTT


def genPDBinPymol(rotationFile, ftentry, masterLig, name):
    rotationStream = readRotations(rotationFile)
    lig, matrix = transformInputs(rotationStream, ftentry, masterLig, name)

    # Formatting stuff
    lig = lig[0]
    matrix = list(matrix.flatten())
    matrix = map(float, matrix)
    cmd.transform_selection("{}".format(lig), matrix)


def genPDB(outPath, clusterFile, ft, rotprm, ligFile, ftGenTemp):
    """Usage: sblu docking gen_cluster_pdb [OPTIONS] CLUSTERFILE FTFILE ROTPRM
                                    LIG_FILE

    Options:
      -o, --output-prefix TEXT    Common prefix for output pdb files (default:
                                  lig)
      -l, --max-clusters INTEGER  Number of top clusters to build models for
                                  (default: all)
      -s, --symmetry TEXT         Type of symmetry
      --help                      Show this message and exit.

    """
    with open(ft, 'r') as rf:
        for num, ftentry in enumerate(rf):
            with open(ftGenTemp, 'x+') as ft_gen_pdb:
                ft_gen_pdb.write(ftentry)
                ft_gen_pdb.write(ftentry)
                ft_gen_pdb.flush()
                cmd = 'sblu docking gen_cluster_pdb -o {out}/{num} {clusterFile} {ft_gen_pdb} {rotprm} {ligFile}'
                print(cmd)
                proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
                so, se = proc.communicate()
