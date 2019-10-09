from prody import *
import numpy as np


class SeleProp:
    def __init__(self, sele):
        self.sele = sele
        self.atom_names = sele.getElements()
        self.atom_nums = sele.numAtoms()
        self.coords = sele.getCoords()
        self.surface = sele.select('surface')

    def vdw_radii(self):
        """
        :return: (array) VdW radii for element names
        """
        vdw_switch = {
            "H": 1.10,
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "F": 1.47,
            "P": 1.80,
            "S": 1.75,
            "CL": 1.75,
            "AS": 1.7
        }
        radii = [vdw_switch[i] for i in self.atom_names]
        return np.asarray(radii)

    def resi_names(self):
        """
        :return: (array) residue names for a selection; each residue is only accounted for once
        """
        res_names = []
        for resi_name in self.sele.iterResidues():
            names.append(resi_name.getResname())
        return res_names

    def hyphb_residue_score(self):
        """
        :return:(array) hydrophobic scores for residue names
        """
        res_names = self.resi_names()
        kd_hyd_switch = {
            "ILE": 4.5,
            "VAL": 4.2,
            "LEU": 3.8,
            "PHE": 2.8,
            "CYS": 2.5,
            "MET": 1.9,
            "ALA": 1.8,
            "GLY": -0.4,
            "THR": -0.7,
            "SER": -0.8,
            "TRP": -0.9,
            "TYR": -1.3,
            "PRO": -1.6,
            "HIS": -3.2,
            "GLU": -3.5,
            "GLN": -3.5,
            "ASP": -3.5,
            "ASN": -3.5,
            "LYS": -3.9,
            "ARG": -4.5
        }

        hyphb_scores = [kd_hyd_switch[i] for i in res_names]
        return np.asarray(hyphb_scores)
