#!/usr/bin/env python
# coding: utf-8


# ./pdb.py

# Description:
# TODO: Make better documenation




from __future__ import division

import os
import linecache
from itertools import islice
import math
from prody import *

from glob import glob
from subprocess import Popen,PIPE

from pymol import cmd, stored, math
import numpy as np



def loadSession(session):
    """
    Loads a given PyMOL Session (e.g., an align-map)
    """
    cmd.load(session)
    pymol.finish_launching()
    cmd.load(session) # otherwise stops @ error 'PyMOL not running, entering library mode (experimental)'




def read_rotations(filepath, limit=None):
    """
    Reads 3x3 rotation matrices from a file.

    Rotations may be a text file with 9 or 10 columns. If the text file has 10
    columns, the first column is assumed to be the line number, which will be discarded.

    Returns a numpy array with dimensions N x (3x3) (an array of 3x3 rotation
    matrices), where N is the smaller of number of rotations in the file and limit,
    if limit is provided.
    """
    with open(filepath, 'r') as stream:
        return read_rotations_stream(stream, limit)

def read_rotations_stream(stream, limit=None):
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




def read_ftresults(filepath, limit=None):
    """
    Reads ftresults from a file.

    See read_ftresults_stream for details.
    """
    with open(filepath, "r") as f:
        return read_ftresults_stream(f, limit)
    
def read_ftresults_stream(stream, limit=None):
    """
    Read ftresults from a stream.

    Ftresults are assumed to be in a text file with at least 5
    columns.  The first column will be the rotation index. The next
    three columns are the translation vector, and the last column is
    the total weighted energy.
    """
    stream = iter(stream)

    return np.loadtxt(
        islice(stream, 0, limit),
        dtype=FTRESULT_DTYPE,
        usecols=(0, 1, 2, 3, 4))




def getCenterandTV(atom_group, ftresults, rotations, center=None):
    orig_coords = atom_group.getCoords()
    if center is None:
        center = np.mean(orig_coords, axis=0)  
    post_tv = np.expand_dims(ftresults['tv'] + center, 1)
    center = center * -1
    return center, post_tv




def getTTT(ft, rot, center, post_tv):
    """
    Arguments include: 
    ft --> output from read_ftresults_stream
    rot --> output from read_rotations_stream
    center --> output from getCenterandTV
    post_tv --> output from getCeneterandTV
    
    Will return a PyMOL happy matrix (4x4) with the top-left most 
    (3x3) is a rotation matrix, the right (3x1) is a post-rotation
    translation vector (accounting for the center) and the bottom
    (1x4) is the pre-rotational matrix ([-1]*center)
    """

    rot_number = ft['roti']
    rot_matrix = rot[rot_number]
    pymolMatrix = np.append(rot_matrix, post_tv, axis = 1)
    zeros = np.append(center, '1')
    pymolMatrix = np.append(pymolMatrix, zeros)
    return pymolMatrix



def transformInputs(rotation_stream, entry, master_lig, lig_name):
    master_lig_pdb = parsePDB(master_lig)
    centers, post_tv = getCenterandTV(master_lig_pdb, entry, rotation_stream)
    TTT = getTTT(entry, rotation_stream, centers, post_tv)
    target = cmd.get_object_list("{}".format(lig_name))
    return target, TTT




def genPDBinPymol(rotationFile, ftentry, masterLig):
    rotationStream = read_rotations(rotationFile)
    cmd.load(masterLig)
    name = os.path.basename(masterLig)
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
                print cmd
                proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
                so,se = proc.communicate()

