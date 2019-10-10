from subprocess import Popen, PIPE

def rmsd(lig, crys_lig, ft, rot, output, options=None):
    """Usage: sblu measure ftrmsd [OPTIONS] LIG_FILE LIG_CRYS FTFILE ROTPRM
    Options:
      --sort-ftresults / --no-sort-ftresults
      -n, --nftresults INTEGER
      --only-CA                       Only C-alpha atoms
      --only-backbone                 Only backbone atoms
      --only-interface                Only use inteface atoms. Requires --rec.
      --interface-radius FLOAT        Radius around receptor to consider.
      --rec PATH                      Receptor file if using interface mode.
      -o, --output FILENAME           Write output to file (default: stdout)
      --help                          Show this message and exit.
    """
    DEFAULT_RMSD_OPTIONS = ""

    if options is None:
        options = DEFAULT_RMSD_OPTIONS

    cmd = 'sblu measure ftrmsd -o {} {} {} {} {} {}'.format(output, options, lig, crys_lig, ft, rot)
    print(cmd)

    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()