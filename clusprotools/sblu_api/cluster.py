from subprocess import Popen, PIPE

def cluster(pwrmsd_file, output, options=None):
    """Usage: sblu docking cluster [OPTIONS] PWRMSDS
    Options:
      -r, --radius FLOAT
      -s, --min-cluster-size INTEGER
      -l, --max-clusters INTEGER
      -o, --output FILENAME
      --json / --no-json
      --help                          Show this message and exit.
    """

    DEFAULT_CLUSTER_OPTIONS = ""

    if options is None:
        options = DEFAULT_CLUSTER_OPTIONS

    cmd = 'sblu docking cluster -o {} {} {}'.format(output, options, pwrmsd_file)
    print(cmd)
    proc = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    so, se = proc.communicate()