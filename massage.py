#!/usr/bin/env python
# coding: utf-8

# ./massage.py

from __future__ import division

import os
from glob import glob
import shutil

from subprocess import Popen, PIPE


def extract(compressed_file, extracted_dir=None):
    """Extracts compressed_file

    If a directory is specified (extracted_dir),then
    extracted_file will be stored there.
    If the directory does not exist, then it will
    be created.
    """
    extracted_file = os.path.join(os.path.basename(compressed_file).replace('.gz', ''))

    if extracted_dir is None:
        cmd = 'gzip -dc {}'.format(compressed_file)
        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd

    elif not os.path.isdir(extracted_dir):
        os.makedirs(extracted_dir)
        extracted_path = os.path.abspath(extracted_dir)
        cmd = 'gzip -dc {} > {}'.format(compressed_file, extracted_path)
        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd

    else:
        extracted_path = os.path.join(extracted_dir, extracted_file)
        cmd = 'gzip -dc {} > {}'.format(compressed_file, extracted_path)
        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd


def idInputs(job_dir, lig_crys=None):
    """Returns the ft files (ft*), ligand (lig),
    rotational matrix (rot_prm).

    Returns a crystallized ligand (lig_crys), which is commonly
    used for RMSD calculations, if specified.
    """
    if os.path.exists(job_dir) is True:
        rec = glob(os.path.join(job_dir, 'rec.pdb.gz'))
        rotprm = glob(os.path.join(job_dir, 'prms/rot70k.0.0.4.prm'))
        lig = glob(os.path.join(job_dir, 'lig.pdb.gz'))
        ft000 = glob(os.path.join(job_dir, 'ft.000.*.gz'))
        ft002 = glob(os.path.join(job_dir, 'ft.002.*.gz'))
        ft004 = glob(os.path.join(job_dir, 'ft.004.*.gz'))
        ft006 = glob(os.path.join(job_dir, 'ft.006.*.gz'))
        ft = ft000 + ft002 + ft004 + ft006

        if lig_crys is None:
            inputs = lig + ft + rec + rotprm
            return inputs

        else:
            lig_crys = glob(lig_crys)
            inputs = lig + ft + rotprm + rec + lig_crys
            return inputs

    else:
        print("Input folder does not have a valid path")


def genSession(job_dir, session_dir, lig_crys=None, multisession=None):
    """Extracts, via extract, inputs returned from idInputs
    to a destination folder (session_dir)

    default: multisession=None, but can set to "Yes" if job_dir
    contains more than one ClusPro docking session folders

    see idInputs for information on lig_crys
    """
    if multisession is None:
        if os.path.exists(os.path.join(job_dir, "prms")) is True:
            session_inputs = idInputs(job_dir, lig_crys)

            if lig_crys is None:
                for item in session_inputs[:-1]:
                    extract(item, session_dir)
                    shutil.copy(session_inputs[-1], session_dir)

            else:
                for item in session_inputs[:-2]:
                    extract(item, session_dir)
                    shutil.copy(session_inputs[-2:], session_dir)
        else:
            print("Double check to see if your attempting to unpack"
                  "a multi-session folder (do: multisession='Yes')")

    else:

        for job in glob(job_dir):
            for i in (glob(os.path.join(folderPath, job))):
                out = (os.path.join(session_dir, job))
                session_inputs = idInputs(i, lig_crys)

                if lig_crys is None:
                    for item in session_inputs[:-1]:
                        extract(item, out)
                        shutil.copy(session_inputs[-1], out)

                else:
                    for item in session_inputs[:-2]:
                        extract(item, out)
                        shutil.copy(session_inputs[-2], out)
                        shutil.copy(session_inputs[-1], out)
