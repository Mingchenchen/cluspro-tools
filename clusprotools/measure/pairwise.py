import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from clusprotools.properties import SeleProp

def pw_distance(sele1, sele2):
    """
    Calculates the pairwise euclidean distances between two sets of atoms
    :param sele1: (ProDy AtomGroup) with dimensions [[N, [xn, yn, zn]]
    :param sele2: (ProDy AtomGroup) with dimensions [[N, [xn, yn, zn]]
    :return: numpy array pairwise distances with dimensions [N, M]
    """
    sele1_coords = SeleProp(sele1).coords
    sele2_coords = SeleProp(sele2).coords

    pw_distances = euclidean_distances(sele1_coords, sele2_coords)
    return np.asarray(pw_distances)


def pw_vdw_sum(sele1, sele2):
    """
    Calculates the VdW pairwise summation for two sets of atoms
    :param sele1: (ProDy AtomGroup) with dimensions [[N, [xn, yn, zn]]
    :param sele2: (ProDy AtomGroup) with dimensions [[N, [xn, yn, zn]]
    :return: numpy array pairwise distances with dimensions [N, M]
    """

    sele1_radii = SeleProp(sele1).vdw_radii()
    sel1_arr = np.array(sele1_radii)

    sele2_radii = SeleProp(sele2).vdw_radii()
    sel2_arr = np.array(sele2_radii)

    vdw_distances = sel1_arr[:, None] + sel2_arr

    return np.asarray(vdw_distances)
