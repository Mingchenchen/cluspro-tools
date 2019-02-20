#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ./inputMassager.py

# Description:
# Takes a ClusPro protein-docking output folder and
# returns useful files for subsequent analysis


from __future__ import division

import os
from glob import glob
import shutil

from subprocess import Popen, PIPE


def extractThings(things, newThingsDir=None):
    """Extracts things

    If a directory is specified (newThingsDir),then
    newThings will be stored there.
    If the directory does not exist, then it will
    be created.
    """
    newThings = os.path.join(os.path.basename(things).replace('.gz', ''))

    if newThingsDir is None:
        cmd = 'gzip -dc {}'.format(things)
        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd

    elif not os.path.isdir(newThingsDir):
        os.makedirs(newThingsDir)
        newThingsPath = os.path.abspath(newThingsDir)
        cmd = 'gzip -dc {} > {}'.format(things, newThingsPath)
        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd
        print "test"

    else:
        newThingsPath = os.path.join(newThingsDir, newThings)
        cmd = 'gzip -dc {} > {}'.format(things, newThingsPath)
        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd


def idInputs(jobFolder, ligCrys=None):
    """Returns the ft files (ft*), ligand (lig),
    rotational matrix (rot_prm).

    Returns a crystallized ligand (ligCrys), which is commonly
    used for RMSD calculations, if specified.
    """
    if os.path.exists(jobFolder) is True:
        rotprm = glob(os.path.join(jobFolder, 'prms/rot70k.0.0.4.prm'))
        lig = glob(os.path.join(jobFolder, 'lig.pdb.gz'))
        ft000 = glob(os.path.join(jobFolder, 'ft.000.*.gz'))
        ft002 = glob(os.path.join(jobFolder, 'ft.002.*.gz'))
        ft004 = glob(os.path.join(jobFolder, 'ft.004.*.gz'))
        ft006 = glob(os.path.join(jobFolder, 'ft.006.*.gz'))
        ft = ft000 + ft002 + ft004 + ft006

        if ligCrys is None:
            inputs = lig + ft + rotprm
            return inputs

        else:
            ligCrys = glob(ligCrys)
            inputs = lig + ft + rotprm + ligCrys
            return inputs

    else:
        print("Input folder does not have a valid path")


def makeSession(jobFolder, sessionFolder, ligCrys=None, multiSession=None):
    """Extracts, via extractThings, inputs returned from idInputs
    to a destination folder (sessionFolder)

    default: multiSession=None, but can set to "Yes" if jobFolder
    contains more than one ClusPro docking session folders

    see idInputs for information on ligCrys
    """
    if multiSession is None:
        if os.path.exists(os.path.join(jobFolder, "prms")) is True:
            sessionInputs = idInputs(jobFolder, ligCrys)

            if ligCrys is None:
                for item in sessionInputs[:-1]:
                    extractThings(item, sessionFolder)
                    shutil.copy(sessionInputs[-1], sessionFolder)

            else:
                for item in sessionInputs[:-2]:
                    extractThings(item, sessionFolder)
                    shutil.copy(sessionInputs[-2:], sessionFolder)
        else:
            print("Double check to see if your attempting to unpack"
                  "a multi-session folder (do: multiSession='Yes')")

    else:
        jobs = []
        folderPath = os.path.abspath(jobFolder)
        jobsName = os.listdir(jobFolder)

        for job in jobsName:
            for i in (glob(os.path.join(folderPath, job))):
                out = (os.path.join(sessionFolder, job))
                sessionInputs = idInputs(i, ligCrys)

                if ligCrys is None:
                    for item in sessionInputs[:-1]:
                        extractThings(item, out)
                        shutil.copy(sessionInputs[-1], out)

                else:
                    for item in sessionInputs[:-2]:
                        extractThings(item, out)
                        shutil.copy(sessionInputs[-2], out)
                        shutil.copy(sessionInputs[-1], out)


def getInputs(job):
    rotprm = glob(os.path.join(job[0], 'rot70k.*'))
    lig = glob(os.path.join(job[0], 'lig.pdb'))
    ft000 = glob(os.path.join(job[0], 'ft.000.*'))
    ft002 = glob(os.path.join(job[0], 'ft.002.*'))
    ft004 = glob(os.path.join(job[0], 'ft.004.*'))
    ft006 = glob(os.path.join(job[0], 'ft.006.*'))
    inputs = ft000 + ft002 + ft004 + ft006 + lig + rotprm
    return inputs


def idSessionFolder(complexFolder):
    jobs = []
    for i in os.listdir(complexFolder):
        jobs.append(os.path.join(complexFolder, i))
    x = []
    for job in jobs:
        y = getInputs(glob(job))
        x.append(y)
    return x
