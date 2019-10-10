from subprocess import Popen, PIPE

def pwrmsd(pdb_file, ft, rot, output, options=None):
    """Usage: sblu measure pwrmsd [OPTIONS] PDB_FILE FTFILE ROTPRM
    Options:
      --sort-ftresults          Sort ftresults before using
      -n, --nftresults INTEGER  Number of ftresults to use
      --only-CA                 Only C-alpha atoms
      --only-backbone           Only backbone atoms
      --only-interface
      --interface-radius FLOAT
      --rec PATH                PDB to use for calculating interface
      -o, --output FILENAME
      --help                    Show this message and exit.
    """
    DEFAULT_PWRMSD_OPTIONS = ""

    if options is None:
        options = DEFAULT_PWRMSD_OPTIONS

    cmd = 'sblu measure pwrmsd -o {} {} {} {} {}'.format(output, options, pdb_file, ft, rot)
    print(cmd)
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()