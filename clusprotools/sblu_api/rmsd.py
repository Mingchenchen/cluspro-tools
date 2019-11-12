from __future__ import print_function

import numpy as np
from scipy.spatial.distance import cdist

from prody import parsePDB

from sblu.ft import (read_rotations_stream, read_ftresults_stream, apply_ftresults_atom_group)
from sblu.rmsd import calc_rmsd


def rmsd(lig_file, lig_crys, ftfile, rotprm, output, sort_ftresults=False, nftresults=None, only_ca=False,
         only_backbone=False, only_interface=False, interface_radius=10.0, rec=None, center_pdb=None):

    if only_interface and rec is None:
        print('Receptor is required')
        return False

    lig = parsePDB(lig_file)
    lig_crys = parsePDB(lig_crys)

    if sort_ftresults:
        ftresults = read_ftresults_stream(ftfile)
        ftresults.sort(order='E', kind='mergesort')  # only mergesort is stable
        ftresults = ftresults[:nftresults]
    else:
        ftresults = read_ftresults_stream(ftfile, limit=nftresults)

    rotations = read_rotations_stream(rotprm)

    if center_pdb:
        center_lig = parsePDB(center_pdb)
        center_lig_coords = center_lig.getCoords()
        center = np.mean(center_lig_coords, axis=0)
        transformed = apply_ftresults_atom_group(lig, ftresults, rotations, center=center)

    else:
        transformed = apply_ftresults_atom_group(lig, ftresults, rotations)

    if only_ca:
        transformed = transformed.select("name CA")
        lig_crys = lig_crys.select("name CA")
    elif only_backbone:
        transformed = transformed.backbone
        lig_crys = lig_crys.backbone

    lig_crys_coords = lig_crys.getCoords()

    interface = None
    if rec and only_interface:
        rec = parsePDB(rec)
        rec_coords = rec.getCoords()

        r_sq = interface_radius ** 2
        dists = cdist(rec_coords, lig_crys_coords, 'sqeuclidean')
        interface = np.any(dists < r_sq, axis=0).nonzero()[0]

    rmsds = calc_rmsd(transformed, lig_crys, interface)
    for rms in rmsds:
        print("{:.2f}".format(rms), file=output)
