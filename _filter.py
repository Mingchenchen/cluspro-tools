#!/usr/bin/env python
# coding: utf-8

# ./_filter.py

# Description:
# TODO: Make better documenation

from __future__ import division
import numpy as np
from pymol import cmd, stored, math
from subprocess import Popen, PIPE
import os
import linecache
from itertools import islice
import math
from prody import *
import pdb
from glob import glob

# gets the atom1, atom2 and max from output_path and runs func(vdwClash)


def doWeClash(output_path, threshold):
    with open(output_path, 'r') as rf:
        line = rf.readline()
        cnt = 1
        clash_count = []
        while line:
            atom1 = (line.split(' ')[0]).split('/')[-1]
            atom2 = (line.split(' ')[2]).split('/')[-1]
            distance = (line.split(' ')[3]).split('\n')[0]
            clash_d = vdwClash(atom1, atom2, distance, threshold)
            if clash_d > threshold:
                clash_count.append(clash_d)
            line = rf.readline()
            cnt += 1
    return (len(clash_count) / cnt)


def pairwise_dist(sel1, sel2, max_dist, output_path, output="N", sidechain="N", show="N"):
    """
    usage: pairwise_dist sel1, sel2, max_dist, [output=S/P/N, [sidechain=N/Y, [show=Y/N]]]
    sel1 and sel2 can be any to pre-existing or newly defined selections
    max_dist: maximum distance in Angstrom between atoms in the two selections

    --optional settings:
    output: accepts Screen/Print/None (default N)
    sidechain: limits (Y) results to sidechain atoms (default N)
    show: shows (Y) individual distances in pymol menu (default=N)
    """

    cmd.delete("dist*")
    extra = ""
    if sidechain == "Y":
        extra = " and not name c+o+n"

    # builds models
    m1 = cmd.get_model(sel2+" around "+str(max_dist)+" and "+sel1+extra)
    m1o = cmd.get_object_list(sel1)
    m2 = cmd.get_model(sel1+" around "+str(max_dist)+" and "+sel2+extra)
    m2o = cmd.get_object_list(sel2)

    # defines selections
    cmd.select("__tsel1a", sel1+" around "+str(max_dist)+" and "+sel2+extra)
    cmd.select("__tsel1", "__tsel1a and "+sel2+extra)
    cmd.select("__tsel2a", sel2+" around "+str(max_dist)+" and "+sel1+extra)
    cmd.select("__tsel2", "__tsel2a and "+sel1+extra)
    cmd.select("IntAtoms_"+max_dist, "__tsel1 or __tsel2")
    cmd.select("IntRes_"+max_dist, "byres IntAtoms_"+max_dist)

    # controlers-1
    if len(m1o) == 0:
        print "warning, '"+sel1+extra+"' does not contain any atoms."
        return
    if len(m2o) == 0:
        print "warning, '"+sel2+extra+"' does not contain any atoms."
        return

    # measures distances
    s = ""
    counter = 0
    for c1 in range(len(m1.atom)):
        for c2 in range(len(m2.atom)):
            distance = math.sqrt(
                sum(map(lambda f: (f[0]-f[1])**2, zip(m1.atom[c1].coord, m2.atom[c2].coord))))
            if distance < float(max_dist):
                s += "%s/%s/%s/%s/%s to %s/%s/%s/%s/%s: %.3f\n" % (m1o[0], m1.atom[c1].chain, m1.atom[c1].resn, m1.atom[c1].resi,
                                                                   m1.atom[c1].name, m2o[0], m2.atom[c2].chain, m2.atom[c2].resn, m2.atom[c2].resi, m2.atom[c2].name, distance)
                counter += 1
                if show == "Y":
                    cmd.distance(m1o[0]+" and "+m1.atom[c1].chain+"/"+m1.atom[c1].resi+"/"+m1.atom[c1].name,
                                 m2o[0]+" and "+m2.atom[c2].chain+"/"+m2.atom[c2].resi+"/"+m2.atom[c2].name)

    # controler-2
    if counter == 0:
        print "warning, no distances were measured! Check your selections/max_dist value"
        return

    # outputs
    if output == "S":
        print s
    if output == "P":
        f = open(output_path, 'w')
        f.write(s)
        f.close()
    # print "Number of distances calculated: %s" % (counter)
    cmd.hide("lines", "IntRes_*")
    if show == "Y":
        cmd.show("lines", "IntRes_"+max_dist)
    cmd.deselect()


cmd.extend("pairwise_dist", pairwise_dist)


def vdwClash(atom1, atom2, dist, threshold):

    h_vdw = 1.1
    c_vdw = 1.7
    n_vdw = 1.55
    o_vdw = 1.52
    f_vdw = 1.47
    p_vdw = 1.8
    s_vdw = 1.8
    cl_vdw = 1.75

    atomone = atom1[0]
    atomtwo = atom2[0]

    if atomone == 'H':
        atom1 = h_vdw
    if atomone == 'C':
        if atom1[0:1] == 'Cl':
            atom1 = cl_vdw
        else:
            atom1 = c_vdw
    if atomone == 'N':
        atom1 = n_vdw
    if atomone == 'O':
        atom1 = o_vdw
    if atomone == 'F':
        atom1 = g_vdw
    if atomone == 'P':
        atom1 = p_vdw
    if atomone == 'S':
        atom1 = s_vdw

    if atomtwo == 'H':
        atom2 = h_vdw
    if atomtwo == 'C':
        if atom2[0:1] == 'Cl':
            atom2 = cl_vdw
        else:
            atom2 = c_vdw
    if atomtwo == 'N':
        atom2 = n_vdw
    if atomtwo == 'O':
        atom2 = o_vdw
    if atomtwo == 'F':
        atom2 = g_vdw
    if atomtwo == 'P':
        atom2 = p_vdw
    if atomtwo == 'S':
        atom2 = s_vdw

    vdw_dist = atom1 + atom2
    clash_dist = vdw_dist*threshold - float(dist)
    return clash_dist


def getInterface(count, chain, pwd1, pwd2):
    MOL2 = cmd.get_object_list("{}".format(count))
    MOL2 = ''.join(MOL2)  # Needed to format so PyMOL is happy
    MOL1 = cmd.get_object_list("rec-mol1")
    MOL1 = ''.join(MOL1)  # Needed to format so PyMOL is happy
    MOBILE = cmd.get_object_list("mobile-mol2")
    MOBILE = ''.join(MOBILE)  # Needed to format so PyMOL is happy
    cmd.align(MOBILE, MOL2)
    cmd.copy("{}.mol2".format(MOL2), MOBILE)
    cmd.delete("{}".format(MOL2))
    cmd.set_name("{}.mol2".format(MOL2), "{}".format(MOL2))
    cmd.select("mol2", "{} and not chain {}".format(MOL2, chain))
    cmd.select("lig_interface1", "rec-mol1 within 10.0 of mol2")
    cmd.select("lig_interface2", "{} within 10.0 of mol1".format(MOL2))
    atomcnt1 = cmd.count_atoms("lig_interface1")
    atomcnt2 = cmd.count_atoms("lig_interface2")
    pairwise_dist("mol2", "lig_interface1", "5",
                  output_path=pwd1, output="P", sidechain="Y")
    pairwise_dist("mol1", "lig_interface2", "5",
                  output_path=pwd2, output="P", sidechain="Y")
    return atomcnt1, atomcnt2, MOL2


def reportFilter(linkerLength, threshold, percent_threshold, count, pwd1, pwd2, atom1, atom2, outputSummary, cnt1, cnt2, MOL2):
    good_linker = []
    good_Clash = []
    linker_distances = []
    atomcnt1 = cnt1
    atomcnt2 = cnt2

    if doWeClash(pwd1, threshold) > percent_threshold:
        pass
    elif doWeClash(pwd2, threshold) > percent_threshold:
        pass
    else:
        dst = cmd.get_distance("{}".format(atom1), "{}".format(atom2))
        if dst < linkerLength:
            with open(outputSummary, "a") as sf:
                outSumEntry = (str(count) + "\t" + str(dst) + "\n")
                sf.write(outSumEntry)
                sf.close()
            good_linker.append(count)
            good_Clash.append(count)
            linker_distances.append(dst)
        else:
            good_Clash.append(count)
