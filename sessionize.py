#!/usr/bin/env python
# coding: utf-8

# ./sessionize.py

# Description:
# Job refers to a raw ClusPro output folders
# Session refers to a folder usable for subsequent analysis

from __future__ import division
import os
from glob import glob
from subprocess import Popen, PIPE
from sblutools import *
import pandas as pd


class session:
    """Do cool stuff with sessions
    """

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

        if crys_lig is None:
            self.crys_lig = self.lig

    def interface_rmsd(self, ft_type, output=None):
        """Generates the interface rmsd for a given ligand/receptor complex and ft ft_type
        ft_type (list): (e.g., 0, 2, 4, 6)
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
        """Generates the interface pwrmsd for a given ligand/receptor complex and ft ft_type
        ft_type (list): (e.g., 0, 2, 4, 6, special)
        num_entries (int)
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

    def near_native(self, threshold, rmsd=10.0):
        """Returns the number of near native structures
        """
        for irmsd_file in self.irmsd:
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
