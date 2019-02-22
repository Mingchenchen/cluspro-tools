from __future__ import division
import os
import sys
from _pdb import *
from _filter3 import *
from inputMassager import *
import pymol
from glob import glob
import time

FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3))])  # Do not change this


def loadSession(session):
    """
    Loads a given PyMOL Session (e.g., an align-map)
    """
    cmd.load(session)
    pymol.finish_launching()
    cmd.load(session)


folders = "/home/ajindal/Desktop/JOBS/6BOY/SESSIONS"

makeSession("/home/ajindal/Desktop/JOBS/6BOY", folders, multiSession="Yes")

sessions = idSessionFolder(folders)
linkerLength = [15, 20]
threshold = [0.4, 0.6, 0.8]
percentAllow = [0.025, 0.05, 0.075]
outputFolder = "/home/ajindal/Desktop/JOBS/6BOY/OUTPUT"
somefile = "/home/ajindal/Desktop/JOBS/6BOY/OUTPUT/test2.txt"

for k in sessions:
    cwd = os.path.basename(os.path.dirname(k[0]))
    ftFiles = k[:4]
    lig = "/home/ajindal/Desktop/JOBS/6BOY/mobile-mol2.pdb"
    rotprm = k[4]
    for ftFile in ftFiles:
        dic = {}
        ftFileStream = readFTResults(ftFile)
        for num, entry in enumerate(ftFileStream[:35000]):
            flag = False
            loadSession("/home/ajindal/Dropbox/PROTAC/6BOY/6BOY-align-map.pse")
            cmd.load(lig, "{}".format(num))
            genPDBinPymol(rotprm, entry, lig, num)
            selection1 = "mol1"
            selection2 = "{} and not chain C".format(num)
            a1list, a2list, dlist = pairwiseDist(selection1, selection2, 10)
            for l in linkerLength:
                for t in threshold:
                    for p in percentAllow:
                        atom1 = "/mol1/E/B/502/C16"
                        atom2 = "/{}/E/B/502/C3".format(num)
                        ft = os.path.basename(ftFile)
                        outputSummary = "{}.{}.{}.{}.{}.txt".format(cwd, ft, l, t, p)
                        flag = reportFilter(a1list, a2list, dlist, t, p, l, atom1,
                                            atom2, outputSummary, num, somefile, dic, 3)
                        if flag:
                            break
                    if flag:
                        break
                if flag:
                    break
            if flag:
                break
