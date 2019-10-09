from __future__ import division
import os
from subprocess import Popen, PIPE


def extract(file_path, destination_path=None, gzip_options=None):
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
    print(cmd)