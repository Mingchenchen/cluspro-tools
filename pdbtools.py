#!/usr/bin/env python
# coding: utf-8

# ./pdbtools.py

# Description:
# come back to this

from __future__ import division
from __future__ import print_function
from prody import *


def get_ligand(cci, ligand_path):
    ligand = fetchPDBLigand("{}".format(cci))
    ligand_model = ligand.get('model')
    writePDB(ligand_path, ligand_model)
