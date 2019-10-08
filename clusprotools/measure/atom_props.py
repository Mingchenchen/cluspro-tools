from prody import *

def vdw_radius(atom_name):
    """
    gets van der waal radius for a particular atom
    :param atom_sel: (string) atom name
    :return: (float) VdW radius for atom name
    """
    vdwSwitch = {
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

    return vdwSwitch(atom_name, 2.0)