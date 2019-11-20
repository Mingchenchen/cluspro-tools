import os
import shutil
import pandas as pd

from glob import glob

from clusprotools.sblu_api import rmsd, pwrmsd, cluster


class Session(object):

    def __init__(self, session_path, raw=False, ft_type=None, ft_file=None):
        """
        gets the files, directories for either a raw or new session
        :param session_path: file path for a session
        :param raw: Bool flag designating raw session
        :param ft_type: string indicating ft type of interest (e.g., '0', '2' or 'custom').
                        will probably be used later to find rmsd/cluster files
        :param ft_file: file path for ft file of interest
        """

        self.session_path = session_path
        self.session_name = os.path.basename(self.session_path)
        os.chdir(self.session_path)

        self._ft_type = ft_type
        self._ft_file = ft_file

        self.raw = raw
        if self.raw:
            # Raw ClusPro result files of interest:
            self.ft0_raw = 'ft.000.00.gz'
            self.ft2_raw = 'ft.002.00.gz'
            self.ft4_raw = 'ft.004.00.gz'
            self.ft6_raw = 'ft.006.00.gz'
            self.lig_raw = 'lig.pdb.gz'
            self.rec_raw = 'rec.pdb.gz'
            self.rot_raw = os.path.join('prms', 'rot70k.0.0.4.prm')
            self._initialized = False
        else:
            # Standard session files of interest:
            self.ft0 = 'ft.000.00'
            self.ft2 = 'ft.002.00'
            self.ft4 = 'ft.004.00'
            self.ft6 = 'ft.006.00'
            self.rot = 'rot70k.0.0.4.prm'
            self.lig = 'lig.pdb'
            self.rec = 'rec.pdb'

            # Output directories
            self.cluster_dir = 'cluster'
            self.rmsd_dir = 'rmsd'
            self.pdb_dir = 'pdb'
            self.ft_dir = 'ft'
            self.result_dir = 'result'
            self.log_dir = 'log'
            self.initialized_dirs = [self.cluster_dir, self.rmsd_dir, self.ft_dir,
                                     self.pdb_dir, self.result_dir, self.log_dir]
            self.initialized = self.initialized_dirs

    @property
    def initialized(self):
        return self._initialized

    @initialized.setter
    def initialized(self, initialized_dirs):
        if len(initialized_dirs) == 6 or len(self.initialized_dirs) == 6:
            self._initialized = True

    #####
    # FT:
    #####

    @property
    def ft_type(self):
        return self._ft_type

    @ft_type.setter
    def ft_type(self, ft_type):
        if ft_type in ['0', '2', '4', '6', 'custom']:
            self._ft_type = ft_type
        else:
            self._ft_type = None

    @property
    def ft_file(self):
        return self._ft_file

    @ft_file.setter
    def ft_file(self, custom_ft=None):
        """
        requires ft_type thru either:
        1.) Session creation (e.g., Session('/session/path', ft_type='0')
        2.) setter ft_type (e.g., SessionObject.ft_type = 'custom')
        :param custom_ft: basename in /ft/ OR absolute file path to custom ft
        :return: ft file path based on either ft_type or custom_ft
        """
        if self._ft_type is not None:
            if self._ft_type == '0':
                self._ft_file = self.ft0
            elif self._ft_type == '2':
                self._ft_file = self.ft2
            elif self._ft_type == '4':
                self._ft_file = self.ft4
            elif self._ft_type == '6':
                self._ft_file = self.ft6
            elif self._ft_type == 'custom' and custom_ft is not None:
                custom_ft_file = self.custom_ft_file(custom_ft)
                self._ft_file = custom_ft_file
        else:
            self._ft_file = None

    def custom_ft_file(self, ft_base_name):
        if self._ft_type == 'custom' and self._initialized:
            self._ft_file = os.path.join(self.ft_dir, ft_base_name)

        # in case we want to use an ft file outside of initialized dir /ft/
        else:
            if os.path.isfile(ft_base_name):
                self._ft_file = ft_base_name
            else:
                return 'ERROR: Enter valid ft file path selected'

    def convert_raw(self, dest_parent_dir, initialize=True, dest_session_name=None, contents=None, delete_src=False):
        """
        Takes specific contents from a raw ClusPro job and stores them in a new session directory
        :param dest_parent_dir: path for destination parent directory
        :param initialize: flag indicating whether to initialize the destination session
        :param dest_session_name: string for destination path basename
        :param contents: list of compressed files to extract to new_session_path
        :param delete_src: flag indicating whether to delete source directory
        """

        from .utilities import extract

        if self.raw:
            default_contents = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                                self.ft6_raw, self.lig_raw, self.rec_raw, self.rot_raw]

            if contents is None:
                contents = default_contents

            if dest_session_name is None:
                dest_session_name = self.session_name

            dest_session_path = os.path.join(dest_parent_dir, dest_session_name)
            try:
                print dest_session_path
                os.mkdir(dest_session_path)
                for item in contents:
                    print item
                    if item.endswith('.gz'):
                        extract(item, dest_session_path)
                    else:
                        shutil.copy(item, dest_session_path)

                if initialize:
                    os.chdir(dest_session_path)
                    cluster_dir = 'cluster'
                    rmsd_dir = 'rmsd'
                    pdb_dir = 'pdb'
                    ft_dir = 'ft'
                    result_dir = 'result'
                    initialized_dirs = [cluster_dir, rmsd_dir, ft_dir, pdb_dir, result_dir]

                    for dir in initialized_dirs:
                        os.mkdir(dir)

                    self.initialized = glob(initialized_dirs)

                if delete_src:
                    shutil.rmtree(self.session_path)

            except:
                print 'WARNING: {} Raw conversion failed'.format(self.session_name)

        else:
            print '{} is not set as raw session. Flag with \'raw=True\' '.format(self.session_path)
            return False

    def gen_rmsd(self, lig_file=None, lig_crys=None, ft_file=None, rotprm=None, output=None, sort_ftresults=False,
                 nftresults=None, only_ca=False, only_backbone=False, only_interface=False, interface_radius=10.0,
                 rec=None, center_pdb=None):

        if lig_file is None:
            lig_file = os.path.abspath(self.lig)

        if lig_crys is None:
            lig_crys = os.path.abspath(self.lig)

        if ft_file is None:
            ft_file = os.path.abspath(self._ft_file)

        if rotprm is None:
            rotprm = os.path.abspath(self.rot)

        if output is None:
            output = os.path.abspath(os.path.join(self.rmsd_dir, 'irmsd.{}'.format(self._ft_type)))

        rmsd(lig_file, lig_crys, ft_file, rotprm, output, sort_ftresults, nftresults, only_ca, only_backbone,
             only_interface, interface_radius, rec, center_pdb)

    def near_native(self, threshold, rmsd=10.0, rmsd_file=None):
        """
        Generates information for near native poses
        :param threshold: number of ft entries to consider
        :param rmsd: number representing maximum rmsd to consider
        :param rmsd_file: default is the interface_rmsd, can be file path to other rmsd file
        :return: session name, threshold, number of near native poses, near native pose entries
        """
        if rmsd_file is None or rmsd_file is '[]':
            rmsd_file = glob(os.path.join(self.rmsd_dir, '*'))

        df = pd.read_csv(rmsd_file, names=["RMSD"])
        df = df[:int(threshold)]
        bad_poses = df[df["RMSD"] > rmsd].index
        df.drop(bad_poses, inplace=True)
        shape = df.shape
        num_hits = shape[0]
        hits = df.index.values

        return self.session_name, num_hits, hits