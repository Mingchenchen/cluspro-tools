#!/usr/bin/env python
# coding: utf-8

# ./_filter.py

# Description:

from __future__ import division
from pymol import cmd
import math


def doWeClash(atoms1, atoms2, distances, tolerance):
    """Determines whether two atoms clash based a softened VdW radii distance.
    Returns a percentage of clashing instances over number of considered atoms.

    Arguments:
    - atoms1 - (list) first set of atoms (from pairwiseDist)
    - atoms1 - (list) second set of atoms (from pariwiseDist)
    - distances - (list) pairwise distances corresponding to atoms1 and atoms2 (from pairwiseDist)
    - tolerance - (float) value indicating how much to soften the VdW radii distance (from User Input)
    """
    clash_count = []

    for i in range(0, len(atoms1)):
        atom1 = atoms1[i][0]
        atom2 = atoms2[i][0]
        distance = distances[i]
        a1 = getVDW(atom1)
        a2 = getVDW(atom2)
        vdwDist = a1 + a2
        clashDist = vdwDist*tolerance - float(distance)
        if clashDist > 0:
            clash_count.append(clashDist)
    try:
        percentClash = (len(clash_count) / len(atoms1))
        return percentClash
    except:
        percentClash = 0
        return percentClash


def pairwiseDist(sel1, sel2, interface):
    """Calculates a pairwise distance between a first and second selection within a defined interface.

    Arguments:
    - sel1 - (str) first PyMOL selection (from User Input)
    - sel2 - (str) second PyMOL selection (from User Input)
    - interface - (float) desired interface selection (from User Input)
    """

    m1 = cmd.get_model(sel2+" around "+str(interface)+" and "+sel1)
    m2 = cmd.get_model(sel1+" around "+str(interface)+" and "+sel2)

    a1List = []
    a2List = []
    dList = []
    for c1 in range(len(m1.atom)):
        for c2 in range(len(m2.atom)):
            distance = math.sqrt(
                sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord, m2.atom[c2].coord))))
            if distance < float(interface):
                a1List.append(m1.atom[c1].name)
                a2List.append(m2.atom[c2].name)
                dList.append(distance)
    return a1List, a2List, dList


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
    }
    return vdwSwitch.get(atomName)


def reportFilter(atoms1, atoms2, distances, tolerance, percentThreshold, linkerLength, atom1, atom2, out, num, somefile, dic, cutoff):
    """Writes to disk the FT entry number and correponding linker length that passes filter parameters

    Arguments:
    - atoms1 - (list) first set of atoms (from pairwiseDist)
    - atoms1 - (list) second set of atoms (from pariwiseDist)
    - tolerance - (float) value indicating how much to soften the VdW radii distance (from User Input)
    - percentThreshold - (float) value indicating how many clashes (%) are allowed to pass (from User Input)
    - linkerLength - (float) value to filter a feasible linker length (from User Input)
    - atom1 - (str) PyMOL formatted atom name for distance calculation between atom2 (from User Input)
    - atom2 - (str) PyMOL formatted atom name for distance calculation between atoms1 (from User Input)
    - out - (path) enumerated file path describing the combination of parameters used for filtering (from User Input)
    - num - (int) value indicating FT entry of interest
    """
    if doWeClash(atoms1, atoms2, distances, tolerance) < percentThreshold:

        dst = cmd.get_distance("{}".format(atom1), "{}".format(atom2))
        if dst < linkerLength:
            if out in dic:
                if dic[out] == int(cutoff):
                    print "wtf"
                    print sum(dic.values())//18
                    if len(dic.keys()) == 18:
                        if sum(dic.values())//18 == int(cutoff):
                            return True
                        else:
                            pass
                    else:
                        pass
                else:
                    dic[out] += 1
                    with open(somefile, "a") as sf:
                        entry = (str(out) + "\t" + str(num) + "\t" + str(dst) + "\n")
                        sf.write(entry)
                    sf.close()
            else:
                dic[out] = 1
                with open(somefile, "a") as sf:
                    entry = (str(out) + "\t" + str(num) + "\t" + str(dst) + "\n")
                    sf.write(entry)
                sf.close()
        else:
            pass
    else:
        pass
