#!/usr/bin/env python
# coding: utf-8

# ./pdb.py

import os

import pymol
from pymol import cmd, math, stored


class PDB:

    def __init__(self, pdb, chain=None, small_molecule=None):

        # Defined Features:
        if os.path.exists(pdb):
            self.pdb_path = pdb
            self.pdb = os.path.basename(pdb).split(".")[0]
            self.get = cmd.load('{}'.format(self.pdb_path))
        else:
            self.pdb = pdb
            self.get = cmd.fetch('{}'.format(self.pdb))

        self.chain = chain
        self.small_molecule = small_molecule
        self.base_name = self.pdb + '.pdb'
        self.selection = 'sele'
        self.select = cmd.select('{}'.format(self.selection), '{}'.format(self.pdb))
        self.sele_macro = '/{}//{}'.format(self.pdb, self.chain)

        # Initialize PyMOL:
        pymol.finish_launching()

    def save(self, save_path, selection_name=None):

        if selection_name is None:
            selection_name = self.selection
            self.select
        else:
            cmd.select('{}'.format(selection_name))

        save_path = os.path.join(save_path, self.base_name)
        cmd.save('{}'.format(save_path), '{}'.format(selection_name))

    def clean(self, selection_name=None, small_molecule=None):

        if selection_name is None:
            selection_name = self.selection

        if small_molecule is None:
            extract_string = '{}'.format(self.sele_macro)

        if small_molecule is not None:

            if small_molecule == self.small_molecule:
                small_molecule = self.small_molecule

            extract_string = '{} or (resn {} w. 10 of {})'.format(
                self.sele_macro, small_molecule, self.sele_macro)

        self.get

        cmd.remove('solvent')
        cmd.extract('{}'.format(selection_name), extract_string)
        print('{}'.format(selection_name), extract_string)

    def align(self, target, pse=None):

        if pse is not None:
            cmd.load(pse)
            self.load
            cmd.align('{}'.format(self.pdb), '{}'.format(target))

        elif os.path.exists(target):
            cmd.load('{}'.format(target))
            self.load
            target_pdb = os.path.basename(target).split(".")[0]
            cmd.align('{}'.format(self.pdb), '{}'.format(target_pdb))

        else:
            cmd.fetch('{}'.format(target))
            self.load
            cmd.align('{}'.format(self.pdb), '{}'.format(target))