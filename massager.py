#!/usr/bin/env python
# coding: utf-8

# ./massager.py

# Description:
# Job refers to a raw ClusPro output folders
# Session refers to a folder usable for subsequent analysis

from __future__ import division
import os
from subprocess import Popen, PIPE
import shutil


class job:
    """A typical ClusPro output directory will look like:
    |-- /output_directory/
        |-- job_1/
            |-- .fts
            |-- .pdbs
            |-- prms/
                |-- .prms
        |-- job_2/
        |-- job_3/
    """

    def __init__(self, job_path):
        self.job_path = job_path
        self.ft0_raw = os.path.join(self.job_path, 'ft.000.00.gz')
        self.ft2_raw = os.path.join(self.job_path, 'ft.002.00.gz')
        self.ft4_raw = os.path.join(self.job_path, 'ft.004.00.gz')
        self.ft6_raw = os.path.join(self.job_path, 'ft.006.00.gz')
        self.rot_raw = os.path.join(self.job_path, 'prms/rot70k.0.0.4.prm')
        self.lig_raw = os.path.join(self.job_path, 'lig.pdb.gz')
        self.rec_raw = os.path.join(self.job_path, 'rec.pdb.gz')

    def extract(self, file_path, destination_path=None, gzip_options=None):
        DEFAULT_GZIP = 'dc'
        """Extracts a file in 'file_path' to a 'destination_path'
        file_path (str): file for extraction
        destination_path (str): extraction destination path
        gzip_options (str): Default set at 'dc'
        """

        DEFAULT_GZIP = 'dc'

        if gzip_options is None:
            gzip_options = DEFAULT_GZIP

        if destination_path is None:
            cmd = 'gzip -{} {}'.format(gzip_options, file_path)

        else:
            if not os.path.isdir(destination_path):
                os.makedirs(destination_path)

            file_base_name = os.path.join(os.path.basename(file_path).replace('.gz', ''))
            new_file_path = os.path.join(destination_path, file_base_name)

            cmd = 'gzip -{} {} > {}'.format(gzip_options, file_path, new_file_path)

        proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
        proc.communicate()
        print cmd

    def convert(self, session_path, contents=None):
        """Makes a session with desired contents
        session_path (str): session destination_path
        contents (list)
        """
        DEFAULT_CONTENTS = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                            self.ft6_raw, self.lig_raw, self.rec_raw, self.rot_raw]

        EXTRACTABLE_CONTENTS = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                                self.ft6_raw, self.lig_raw, self.rec_raw]
        if contents is None:
            contents = DEFAULT_CONTENTS

        for item in contents:

            if item in EXTRACTABLE_CONTENTS:
                self.extract(item, session_path)
            elif item == self.rot_raw:
                shutil.copy(self.rot_raw, session_path)
