#!/usr/bin/env python
# coding: utf-8

# ./report/filtering_parameters.py

# The purpose of this script is to obtain data to representing:
# clashes (and their severity), end point distances, and the corresponding rmsd

import os
import time

import pandas as pd
from glob import glob

from clusprotools.measure import distance, id_binding_interface, clash_count
from clusprotools.session import Session, Component



def main(job, param, rmsd_threshold, ft_threshold, mol_paths, clash_severities):
    """
    :param job: path for job name
    :param param: list of param  [small molecule name, ligand-endpoint, receptor-endpoint, ft file] (e.g., ["resname RN3", "serial 26", "serial 9", "004"])
    :param rmsd_threshold: float of maximum rmsd to consider for a given ft
    :param ft_threshold: integer of maximum ft entries to consider
    :param mol_paths: string of component base names (e.g., "mobile-mol2.pdb")
    :param clash_severities: list of clash severeties to report
    :return: a dictionary of clashes for each clash severity, the end point distances, and irmsd for a given pose
    """

    # Unpack parameters:

    lig_ep = param[1]
    rec_ep = param[2]
    ft_type = param[3]

    # Sessionize:
    sesh = Session(job)

    # Get low RMSD poses:
    summary = sesh.near_native(threshold=ft_threshold, rmsd=rmsd_threshold)
    hit_index = summary[3]

    # Get Atom Groups:
    ligand = Component(os.path.join(sesh.session_path, mol_paths[0]))

    receptor = Component(os.path.join(sesh.session_path, mol_paths[1]))
    receptor = receptor.atom_group

    og_ligand = Component(os.path.join(sesh.session_path, mol_paths[2]))
    og_ligand = og_ligand.atom_group

    # Transform ligand
    transformed = ligand.xform(ft_type, center=og_ligand)

    # Iterate through FT file with hit_index
    hit_d = {}
    for i in hit_index:

        # Define Dictionary Key:
        basename = os.path.basename(job)
        k = basename + "-" + str(i)

        # Define interfaces:
        transformed.setACSIndex(i)
        lig_interface, rec_interface = id_binding_interface(transformed, receptor)

        # Calculate number of atoms at interface:
        num_atoms_ilig = lig_interface.numAtoms()
        num_atoms_irec = rec_interface.numAtoms()
        num_atoms_interface = num_atoms_ilig + num_atoms_irec

        # Get clash info:
        clash_info = []
        clash_info.append(num_atoms_interface)
        for cc in clash_severities:
            clashes = clash_count(lig_interface, rec_interface, threshold=cc)
            clash_info.append(clashes)

        irmsd = ligand.pose_rmsd(i)
        clash_info.append(irmsd[0])

        # Select End Points:
        ep_1 = transformed.select(lig_ep)
        ep_2 = receptor.select(rec_ep)

        # Calculate end point distance:
        ep_d = distance(ep_1, ep_2)

        clash_info.append(round(ep_d[0], 2))

        hit_d[k] = clash_info

    return hit_d


def filtering_parameters(clash_analysis_inputs, labels, rmsd_threshold, ft_threshold, mol_paths, clash_severities):
    """
    Generates csv in each session folder with filter parameter data
   :param job: path for job name
   :param param: list of param  [small molecule name, ligand-endpoint, receptor-endpoint, ft file]
   (e.g., ["resname RN3", "serial 26", "serial 9", "004"])
   :param rmsd_threshold: float of maximum rmsd to consider for a given ft
   :param ft_threshold: integer of maximum ft entries to consider
   :param mol_paths: string of component base names (e.g., "mobile-mol2.pdb")
   :param clash_severities: list of clash severeties to report
    """
    for dir_path, params in clash_analysis_inputs:
        glob_dir = (os.path.join(dir_path, "*"))

        for job in glob(glob_dir):
            start = time.time()

            print("Session Path: {} is being processed".format(job))

            # Execute main function:
            data = main(job, params, rmsd_threshold, ft_threshold, mol_paths, clash_severities)
            df = pd.DataFrame.from_dict(data, orient='index')

            if len(df.columns) == len(labels):
                df.columns = labels

            # Export Data Frame to csv:
            csv_path = os.path.join(job, "clash_info.csv")
            df.to_csv(csv_path, sep="\t")

            end = time.time()
            total_time = end - start
            print("Collecting data for this session took {} seconds".format(round(total_time, 2)))
