#!/usr/bin/env python
# coding: utf-8

# ./transform.py

# Description:
# Transforms proteins using rotational and translational matricies


from itertools import islice
import linecache

import numpy as np

FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3)), ('E', 'f8')])


def read_rotations(filepath, limit=None):
    """Reads 3x3 rotation matrices from a file.

    Rotations may be a text file with 9 or 10 columns. If the text file has 10
    columns, the first column is assumed to be the line number, which will be discarded.

    Returns a numpy array with dimensions N x (3x3) (an array of 3x3 rotation
    matrices), where N is the smaller of number of rotations in the file and limit,
    if limit is provided.
    """
    with open(filepath, 'r') as stream:
        return read_rotations_stream(stream, limit)


def read_rotations_stream(stream, limit=None):
    """Read rotations from a stream.

    Rotations may be a text file with 9 or 10 columns. If the text file has 10
    columns, the first column is assumed to be the line number, which will be discarded.

    Returns a numpy array with dimensions N x (3x3) (an array of 3x3 rotation
    matrices), where N is the smaller of number of rotations in the file and limit,
    if limit is provided.
    """
    rotations = np.loadtxt(
        islice(iter(stream), 0, limit))
    if rotations.shape[-1] == 10:
        rotations = rotations[:, 1:]
    return rotations.reshape(-1, 3, 3)


def read_ftresults(filepath, limit=None):
    """Reads ftresults from a file.

    See read_ftresults_stream for details."""
    with open(filepath, "r") as f:
        return read_ftresults_stream(f, limit)


def read_ftresults_stream(stream, limit=None):
    """Read ftresults from a stream.

    Ftresults are assumed to be in a text file with at least 5
    columns.  The first column will be the rotation index. The next
    three columns are the translation vector, and the last column is
    the total weighted energy.
    """
    stream = iter(stream)

    return np.loadtxt(
        islice(stream, 0, limit),
        dtype=FTRESULT_DTYPE,
        usecols=(0, 1, 2, 3, 4))


def get_ftresult(filepath, index):
    """Get ftresult at index from file.

    index should be zero offset.
    """
    line = linecache.getline(filepath, index + 1)
    if not line:
        return None

    tokens = line.strip().split()
    return np.array(
        (int(tokens[0]), [float(c) for c in tokens[1:4]], float(tokens[4])),
        dtype=FTRESULT_DTYPE)


def apply_ftresult(coords, ftresult, rotations, center=None, out=None):
    """Apply the ftresult to coords.

    `coords` and `out` cannot point to the same numpy array.
    """
    if center is None:
        center = np.mean(coords, axis=0)

    if out is None:
        out = np.empty_like(coords)

    out = np.dot(coords - center, rotations[ftresult['roti']].T)
    np.add(out, ftresult['tv'] + center, out)

    return out


def symmetrize_ftresults(ftresults, rotations, subunit):
    try:
        if len(ftresults.shape) == 0:
            ftresults = np.expand_dims(ftresults, 0)
    except:
        raise ValueError("ftresults does not seem to be an ndarray or void")

    loc_ft = ftresults.copy()
    mov_ft = loc_ft.copy()
    rot = rotations[ftresults['roti']].copy()

    # To construct subunits for C-n symmetry, translations are calculated by
    # rotating the translation vector and adding it to the previous
    # translation.
    # A way to visualize this is drawing a polygon. We have the first side of
    # the polygon as our initial translation. We then rotate that vector and
    # add it to the initial vector to draw the next side and get to the
    # position of the next subunit. And so on.
    # How do we rotate it? We use the same rotation matrix we would apply to
    # the initial subunit. E.g. for a trimer, that rotation is 120, so we are
    # rotating the vector # according to the interior angle of an equilateral
    # polygon.
    for _ in range(subunit-1):
        # we have a vector for each input transformation, rotate that vector by
        # the corresponding rotation
        for i, vec in enumerate(mov_ft['tv']):
            mov_ft['tv'][i] = np.dot(rot[i], vec.T).transpose()
        loc_ft['tv'] += mov_ft['tv']
    # Rotations are constructed by multiplying the rotation matrix by itself
    # the number of times necessary for a subunit.
    for i, mat in enumerate(rot):
        rot[i] = np.linalg.matrix_power(mat, subunit)
        loc_ft['roti'][i] = i
    return loc_ft, rot


def apply_ftresults_atom_group(atom_group,
                               ftresults,
                               rotations,
                               center=None,
                               out=None):
    """Apply ftresult(s) to an atomgroup, returning a new atomgroup.

    The new atomgroup will have one coordinate set for each ftresult passed.
    ftresult can either be a single ftresult object, or an array of ftresult
    objects."""
    orig_coords = atom_group.getCoords()  # This returns a copy so we can mutate it
    if center is None:
        center = np.mean(orig_coords, axis=0)
    np.subtract(orig_coords, center, orig_coords)

    try:
        if len(ftresults.shape) == 0:
            ftresults = np.expand_dims(ftresults, 0)
    except:
        raise ValueError("ftresults does not seem to be an ndarray or void")

    if out is None:
        out = atom_group.copy()

    new_coords = np.dot(rotations[ftresults['roti']],
                        orig_coords.T).transpose(0, 2, 1)
    np.add(new_coords, np.expand_dims(ftresults['tv'] + center, 1), new_coords)
    out._setCoords(new_coords, overwrite=True)

    return out
