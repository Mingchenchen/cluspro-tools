#!/usr/bin/env python
# coding: utf-8

# ./session.py


from __future__ import division
import os
from glob import glob

import numpy as np
import pandas as pd
from prody import *

from .sblu import pwrmsd, rmsd
from .io import extract
from .transform import read_ftresults, read_rotations, apply_ftresults_atom_group


class RawSession:

    def __init__(self, raw_path, crys_lig=None):

        # Raw Files:
        self.raw_path = raw_path
        self.ft0_raw = os.path.join(self.raw_path, 'ft.000.00.gz')
        self.ft2_raw = os.path.join(self.raw_path, 'ft.002.00.gz')
        self.ft4_raw = os.path.join(self.raw_path, 'ft.004.00.gz')
        self.ft6_raw = os.path.join(self.raw_path, 'ft.006.00.gz')
        self.rot_raw = os.path.join(self.raw_path, 'prms/rot70k.0.0.4.prm')
        self.lig_raw = os.path.join(self.raw_path, 'lig.pdb.gz')
        self.rec_raw = os.path.join(self.raw_path, 'rec.pdb.gz')

        if crys_lig is None:
            self.crys_lig = self.lig

    def convert(self, session_path, contents=None):
        """
        Converts a raw session comprising a raw ClusPro output into a "working session"
        :param session_path: path of destination for working session directory
        :param contents: list of contents to convert, default = None
        """
        DEFAULT_CONTENTS = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                            self.ft6_raw, self.lig_raw, self.rec_raw, self.rot_raw]

        EXTRACTABLE_CONTENTS = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                                self.ft6_raw, self.lig_raw, self.rec_raw]
        if contents is None:
            contents = DEFAULT_CONTENTS

        for item in contents:

            if item in EXTRACTABLE_CONTENTS:
                extract(item, session_path)
            elif item == self.rot_raw:
                shutil.copy(self.rot_raw, session_path)


class Session:

    def __init__(self, session_path, crys_lig=None):

        # Standard:
        self.session_path = session_path
        self.ft0 = os.path.join(self.session_path, 'ft.000.00')
        self.ft2 = os.path.join(self.session_path, 'ft.002.00')
        self.ft4 = os.path.join(self.session_path, 'ft.004.00')
        self.ft6 = os.path.join(self.session_path, 'ft.006.00')
        self.rot = os.path.join(self.session_path, 'rot70k.0.0.4.prm')
        self.lig = os.path.join(self.session_path, 'lig.pdb')
        self.rec = os.path.join(self.session_path, 'rec.pdb')
        self.session_name = os.path.basename(self.session_path)

        # Optional:
        self.mobile_lig = os.path.join(self.session_path, 'mobile_lig.pdb')
        self.align_map = glob(os.path.join(self.session_path, '*.pse'))
        self.irmsd = glob(os.path.join(self.session_path, '*.irmsd'))
        self.ipwrmsd = glob(os.path.join(self.session_path, '*.ipwrmsd'))
        self.pwrmsd = glob(os.path.join(self.session_path, '*.pwrmsd'))
        self.rmsd = glob(os.path.join(self.session_path, '*.rmsd'))

        if crys_lig is None:
            self.crys_lig = self.lig

    def interface_rmsd(self, ft_type, output=None):
        """
        Generates the interface rmsd for a given ligand/receptor complex and ft_type
        :param ft_type: string list for corresponding ft type (e.g., ["0"] or ["0", "6"]
        :param output: default is in the same session folder
        """

        option = "--only-interface --rec {}".format(self.rec)

        for ft in ft_type:
            DEFAULT_OUTPUT = os.path.join(
                self.session_path, "{}.{}.irmsd".format(self.session_name, ft))

            if output is None:
                output = DEFAULT_OUTPUT

            if ft == "0":
                rmsd(self.lig, self.crys_lig, self.ft0, self.rot, output, option)
            elif ft == "2":
                rmsd(self.lig, self.crys_lig, self.ft2, self.rot, output, option)
            elif ft == "4":
                rmsd(self.lig, self.crys_lig, self.ft4, self.rot, output, option)
            elif ft == "6":
                rmsd(self.lig, self.crys_lig, self.ft6, self.rot, output, option)
            else:
                print("ft_type should be in a list (e.g., ['0'] or ['2', '6'])")

    def interface_pwrmsd(self, ft_type, output=None, num_entries=None, ft_special=None):

        """
        Generates the interface pwrmsd for a given ligand/receptor complex and ft_type
        :param ft_type: string list for corresponding ft type (e.g., ["0"] or ["0", "special"]
        :param output: default is in the same session folder
        :param num_entries: default is 1000 ft entries used for interface_pwrmsd calculation
        :param ft_special: file path for ft_type="special" flag
        """

        DEFAULT_NUM_ENTRIES = 1000
        if num_entries is None:
            num_entries = DEFAULT_NUM_ENTRIES

        option = "--only-interface --rec {} --nftresults {}".format(self.rec, num_entries)

        for ft in ft_type:
            DEFAULT_OUTPUT = os.path.join(
                self.session_path, "{}.{}.ipwrmsd".format(self.session_name, ft))

            if output is None:
                output = DEFAULT_OUTPUT

            if ft == "0":
                pwrmsd(self.lig, self.ft0, self.rot, output, option)
            elif ft == "2":
                pwrmsd(self.lig, self.ft2, self.rot, output, option)
            elif ft == "4":
                pwrmsd(self.lig, self.ft4, self.rot, output, option)
            elif ft == "6":
                pwrmsd(self.lig, self.ft6, self.rot, output, option)
            elif ft == "special":
                pwrmsd(self.lig, ft_special, self.rot, output, option)
            else:
                print("ft_type should be in a list (e.g., ['0'] or ['2', 'special'])")

    def near_native(self, threshold, rmsd=10.0, rmsd_file=None):
        """
        Generates information for near native poses
        :param threshold: number of ft entries to consider
        :param rmsd: number representing maximum rmsd to consider
        :param rmsd_file: default is the interface_rmsd, can be file path to other rmsd file
        :return: session name, threshold, number of near native poses, near native pose entries
        """
        if rmsd_file is None:
            rmsd_file = self.irmsd

        for irmsd_file in rmsd_file:
            df = pd.read_csv(irmsd_file, names=["RMSD"])
            df = df[:int(threshold)]
            bad_poses = df[df["RMSD"] > rmsd].index
            df.drop(bad_poses, inplace=True)
            shape = df.shape
            num_hits = shape[0]
            hits = df.index.values
            base_name = os.path.basename(irmsd_file)
            base_name = base_name.split(".")[0]
            return base_name, threshold, num_hits, hits

    def get_ft(self, ft_type):
        """
        Returns the ft file based on an entered ft_type
        :param ft_type: string representing ft-type (e.g., "000")
        :return: requested ft file path
        """

        if ft_type == "000":
            ft = self.ft0
        if ft_type == "002":
            ft = self.ft2
        if ft_type == "004":
            ft = self.ft4
        if ft_type == "006":
            ft = self.ft6

        return ft

    def from_component(self, component_path):

        return Session(os.path.split(component_path)[0])


class Component(Session):

    def __init__(self, component_path):
        self.component_path = component_path
        self.session_path = os.path.split(component_path)[0]
        self.atom_group = parsePDB(self.component_path)

    # def atom_group(self):
    #     """
    #     Returns an atom group for a specified molecule
    #     :return: the ProDy atom group
    #     """
    #     ag = parsePDB(self.component_path)
    #     return ag

    def pose_rmsd(self, num_entry, rmsd_file=None):
        """
        Returns the rmsd for a specified pose
        :param num_entry: int of ft entry
        :param rmsd_file: default is the interface_rmsd, can be file path to other rmsd file
        :return: rmsd for the pose @ the specified num_entry
        """
        if rmsd_file is None:
            session = super().from_component(self.component_path)
            rmsd_file = session.irmsd

        for rfile in rmsd_file:
            df = pd.read_csv(rfile, names=["RMSD"])
            return df.iloc[num_entry]

    def xform(self, ft_type, center=None):
        """
        Transform a molecule based on translation and rotation matricies
        :param ft_type: string representing ft-type (e.g., "000")
        :return: transformed coordinate of component
        """
        session = super().from_component(self.component_path)
        ft_path = session.get_ft(ft_type)
        rot_path = session.rot
        protein = self.atom_group

        if center is not None:
            center_coords = center.getCoords()
            center = np.mean(center_coords, axis=0)

        ft_results = read_ftresults(ft_path)
        rot_results = read_rotations(rot_path)

        transformed = apply_ftresults_atom_group(protein, ft_results, rot_results, center)

        return transformed
