from __future__ import print_function

import numpy as np
from prody import parsePDB

from sblu.rmsd import interface_pwrmsd
from sblu.rmsd import pwrmsd as pwrmsd_func
from sblu.ft import (read_rotations_stream, read_ftresults_stream, apply_ftresults_atom_group)


def pwrmsd(pdb_file, ftfile, rotprm, output, sort_ftresults=False, nftresults=None, only_ca=False,
         only_backbone=False, only_interface=False, interface_radius=10.0, rec=None, center_pdb=None):

    if only_interface and rec is None:
        print('Receptor is required')
        return False

    lig = parsePDB(pdb_file)

    if sort_ftresults:
        ftresults = read_ftresults_stream(ftfile, limit=None)
        ftresults.sort(order='E', kind='mergesort')
        if nftresults != -1: ftresults = ftresults[:nftresults]
    else:
        ftresults = read_ftresults_stream(ftfile, limit=nftresults)

    rotations = read_rotations_stream(rotprm)

    if center_pdb:
        center_lig = parsePDB(center_pdb)
        center_lig_coords = center_lig.getCoords()
        center = np.mean(center_lig_coords, axis=0)
    else:
        center = np.mean(lig._getCoords(), axis=0)
    if only_ca:
        lig = lig.calpha
    elif only_backbone:
        lig = lig.backbone

    transformed = apply_ftresults_atom_group(lig, ftresults, rotations, center=center)

    if only_interface:
        rec = parsePDB(rec)

        pairwise_dists = interface_pwrmsd(rec, transformed, interface_d=interface_radius)
    else:
        pairwise_dists = pwrmsd_func(transformed)

    for rms in pairwise_dists.flat:
        output.write("{:.2f}\n".format(rms))