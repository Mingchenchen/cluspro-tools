
def hydrophobicity(residue_name):
    """
    gets the hydrophobicity attribute of a particular residue according to values set forth by Kyte and Doolittle
    :param residue_name: (string) three letter residue name
    :return: (float) hydrophobicity attribute
    """
    kdHydSwitch = {
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

    return kdHydSwitch(residue_name, 0.0)