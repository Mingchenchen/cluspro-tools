import os
import shutil

from glob import glob

from clusprotools.sblu_api import rmsd, pwrmsd, cluster


class Session:

    def __init__(self, session_path, raw=False, custom=False):

        # Path info
        self.session_path = session_path
        self.session_name = os.path.basename(self.session_path)
        os.chdir(self.session_path)

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

        else:
            # Standard session files of interest:
            self.ft0 = 'ft.000.00'
            self.ft2 = 'ft.002.00'
            self.ft4 = 'ft.004.00'
            self.ft6 = 'ft.006.00'
            self.rot = 'rot70k.0.0.4.prm'
            self.lig = 'lig.pdb'
            self.rec = 'rec.pdb'

        # Custom Session Information:
        self.custom = custom
        self.ft_file = None
        self.lig_crys = None
        self.lig_custom = None
        self.rec_custom = None

        if custom:
            allowed_keys = {'lig_crys', 'lig_custom', 'rec_custom', 'ft_custom', 'rmsd_custom',
                            'pwrmsd_custom', 'coord_custom', 'pdb_custom', 'cluster_dir', 'ft_dir', 'pdb_dir',
                            'result_dir'}
            for key, value in kwargs.items():
                setattr(self, key, value)

                if key not in allowed_keys:
                    print 'WARNING: Session attribute {} is not allowed.  Proceed with caution.'.format(key)

    @property
    def ft_file(self):
        return self.ft_file

    @property
    def lig_file(self):
        return self.lig

    @ft_file.setter
    def ft_file(self, ft_type, ft_file_path=None, ft_dir_path=None):

        if os.path.isfile(ft_file_path):
            self.ft_file = ft_file_path

        if os.path.isdir(ft_dir_path):
            self.ft_file = glob(os.path.join(ft_dir_path, '*.{}.ft'.format(ft_type)))[0]

        else:
            if ft_type == '0':
                self.ft_file = self.ft0
            if ft_type == '2':
                self.ft_file = self.ft2
            if ft_type == '4':
                self.ft_file = self.ft4
            if ft_type == '6':
                self.ft_file = self.ft6

    @lig_file.setter
    def lig_file(self, lig_path):
        self.lig = lig_path

    def convert_raw(self, dest_parent_dir, initialize=True, dest_session_name=None, contents=None, delete_src=False):
        """
        Takes specific contents from a raw ClusPro job and stores them in a new session directory
        :param dest_parent_dir: path for destination parent directory
        :param initialize: flag indicating whether to intialize the destination session
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
            convert_complete = False

            try:
                os.mkdir(dest_session_path)

                for item in contents:
                    if item.endswith('.gz'):
                        extract(item, dest_session_path)

                    else:
                        shutil.copy(item, dest_session_path)

                convert_complete = True

            except:
                print 'WARNING: {} Raw conversion failed'.format(self.session_name)

            if initialize:
                os.chdir(dest_session_path)
                cluster_dir = 'cluster'
                rmsd_dir = 'rmsd'
                pdb_dir = 'pdb'
                ft_dir = 'ft'
                result_dir = 'result'
                standard_dirs = [cluster_dir, rmsd_dir, ft_dir, pdb_dir, result_dir]

                for dir in standard_dirs:
                    os.mkdir(dir)

            if delete_src and convert_complete:
                shutil.rmtree(self.session_path)

        else:
            print '{} is not set as raw session. Flag with \'raw=True\' '.format(self.session_path)
            return False

    # def gen_rmsd(self, lig_file=None, lig_crys=None, ft_file=None, rotprm=None, output=None, sort_ftresults=False,
    #              nftresults=None, only_ca=False, only_backbone=False, only_interface=False, interface_radius=10.0,
    #              rec=None, center_pdb=None):
    #
    #     input_check = False
    #
    #     if self.ft_file is None and ft_file is None:
    #         print 'ERROR: No ft file selected'
    #
    #     if self.lig is None and lig_file is None:
    #         print 'ERROR: No lig_file selected'

