#!/usr/bin/env python
# coding: utf-8

# ./filter.py

# Description:
# todo: should probably delete this

import pandas as pd

from measure import *
from pdbtools import *


def get_selections(job, mol1_endpoint, mol2_endpoint, ft_file):

    # Sessionize:
    sesh = session(job)

    # Get translation and rotation files:
    rotfile = sesh.rot

    if ft_file == "000":
        ft = sesh.ft0
    if ft_file == "002":
        ft = sesh.ft2
    if ft_file == "004":
        ft = sesh.ft4
    if ft_file == "006":
        ft = sesh.ft6

    # Get PDBs:
    ligand = os.path.join(sesh.session_path, "mobile-mol2.pdb")
    receptor = os.path.join(sesh.session_path, "rec_mol2.pdb")
    orig_ligand = os.path.join(sesh.session_path, "lig.pdb")

    # Define prody objects
    lig = parsePDB(ligand)
    rec = parsePDB(receptor)
    org_lig = parsePDB(orig_ligand)

    # Define selections:
    mol1_endpoint = lig.select(mol1_endpoint)
    mol2_endpoint = rec.select(mol2_endpoint)

    # Get interfaces:
    rec.get_interface()



def clash_count(job, small_molecule, mol1_endpoint, mol2_endpoint, ft_file, clash_threshold=.5, rmsd_threshold=3,
                print_out="y"):
    """Only for PROTAC project (05/20/19)
    """
    sesh = session(job)

    # Get low RMSD poses:
    summary = sesh.near_native(threshold=70000, rmsd=rmsd_threshold)
    num_hits = summary[2]
    hit_index = summary[3]

    # Get files:
    rotfile = sesh.rot

    if ft_file == "000":
        ft = sesh.ft0
    if ft_file == "002":
        ft = sesh.ft2
    if ft_file == "004":
        ft = sesh.ft4
    if ft_file == "006":
        ft = sesh.ft6

    # Define objects:
    ligand = os.path.join(sesh.session_path, "mobile-mol2.pdb")
    receptor = os.path.join(sesh.session_path, "rec_mol2.pdb")
    orig_ligand = os.path.join(sesh.session_path, "lig.pdb")

    lig = parsePDB(ligand)
    rec = parsePDB(receptor)
    org_lig = parsePDB(orig_ligand)

    # define selections:
    mol1_endpoint = lig.select(mol1_endpoint)
    mol2_endpoint = rec.select(mol2_endpoint)
    rec_interface = rec.select('all within 10 of center', center=mol2_endpoint)

    # define centers:
    orig_coords = org_lig.getCoords()
    orig_center = np.mean(orig_coords, axis=0)

    # parse files:
    ftresults = read_ftresults(ft)
    rotresults = read_rotations(rotfile)

    # Apply transformation:
    transformed = apply_ftresults_atom_group(lig, ftresults, rotresults, center=orig_center)
    transformed_interface = transformed.select('all within 20 of center', center=mol1_endpoint)

    # getting the interface for all of the transformed ligands
    if print_out == "y":
        print("--\tJob: {}\t|\tRMSD Threshold: {}A\t|\tTotal RMSD Hits: {}\t--".format(
            summary[0], rmsd_threshold, num_hits))
        print("ft entry:\t # of clashes:\t # of severe clashes ({})".format(clash_threshold))

    # output lists:
    ft_entries = []
    clash_counts = []
    severe_counts = []

    # main:
    for i in hit_index:
        flag = True
        clash_count = 0
        severe_clash = 0
        pw = pairwise_dist(rec_interface, transformed_interface, sel2_index=i)
        for clash_dist, pw_dis in pw:
            if clash_dist > pw_dis:
                flag = False
                clash_count += 1
                if clash_dist * float(clash_threshold) > pw_dis:
                    severe_clash += 1

        if print_out == "y":
            if flag:
                clash_count = 0
                severe_clash = 0
            print("{}\t\t\t{}\t\t{}".format(i, clash_count, severe_clash))
            ft_entries.append(i)
            clash_counts.append(clash_count)
            severe_counts.append(severe_clash)
        else:
            if flag:
                pass
            else:
                ft_entries.append(i)
                clash_counts.append(clash_count)
                severe_counts.append(severe_clash)

    cc_out = ({'Job': [summary[0]], '<{}A Hits'.format(rmsd_threshold): [num_hits],
               'Avg Clash': [round(np.mean(clash_counts), 2)],
               'Avg Severe Clash ({})'.format(clash_threshold): [round(np.mean(severe_counts), 2)]})
    df = pd.DataFrame(cc_out)
    return df