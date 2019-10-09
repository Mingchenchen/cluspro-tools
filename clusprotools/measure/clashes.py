from .pairwise import pw_distance, pw_vdw_sum

def clashes(sele1, sele2, threshold=1.0):
    """
    Calculates the number of clashes for two sets of atoms
    :param sele1: first atom group
    :param sele2: second atom group
    :param threshold: VdW softener
    :return: number of clash counts
    """
    pwd = pw_distance(sele1, sele2)
    vdw_d = pw_vdw_sum(sele1, sele2)

    pwd = pwd.flatten()
    vdw_d = vdw_d.flatten()

    return sum(vdw_d * threshold > pwd)
