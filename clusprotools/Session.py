import os
from glob import glob
import sblu_api.rmsd as _rmsd
import sblu_api.pwrmsd as _pwrmsd
import sblu_api.cluster as _cluster


class Session:

    def __init__(self, session_path, crys_lig=None, type_ft=None, custom_lig=None, custom_rec=None, custom_ft=None,
                 custom_rmsd=None, custom_pwrmsd=None, dir_custom_ft=None, dir_custom_pwrmsd=None, dir_custom_rmsd=None,
                 file_ep_coords=None):

        # Standard Files:
        self.session_path = session_path
        self.ft0 = os.path.join(self.session_path, 'ft.000.00')
        self.ft2 = os.path.join(self.session_path, 'ft.002.00')
        self.ft4 = os.path.join(self.session_path, 'ft.004.00')
        self.ft6 = os.path.join(self.session_path, 'ft.006.00')
        self.rot = os.path.join(self.session_path, 'rot70k.0.0.4.prm')
        self.lig = os.path.join(self.session_path, 'lig.pdb')
        self.rec = os.path.join(self.session_path, 'rec.pdb')

        if crys_lig is None:
            self.crys_lig = self.lig

        # Standard Directories:
        self.dir_cluster = os.path.join(self.session_path, 'cluster')
        self.dir_pwrmsd = os.path.join(self.session_path, 'pwrmsd')
        self.dir_custom_ft = os.path.join(self.session_path, 'custom_ft')
        self.dir_rmsd = os.path.join(self.session_path, 'rmsd')

        standard_directories = [self.dir_cluster, self.dir_rmsd, self.dir_custom_ft, self.dir_pwrmsd]

        for directory in standard_directories:
            if not os.path.exists(directory):
                os.mkdir(directory)

        # Custom Information:
        self.type_ft = type_ft

        # Custom Files:
        self.custom_ft = custom_ft
        self.custom_lig = custom_lig
        self.custom_rec = custom_rec
        self.custom_rmsd = custom_rmsd
        self.custom_pwrmsd = custom_pwrmsd
        self.file_ep_coords = file_ep_coords

        # Custom Directories:
        self.dir_custom_ft = dir_custom_ft
        self.dir_custom_pwrmsd = dir_custom_pwrmsd
        self.dir_custom_rmsd = dir_custom_rmsd

    def id_file(self, file_type, custom_file=False):
        if self.type_ft is None:
            print 'No ft type is selected'
            return None

        type_ft = self.type_ft
        if len(self.type_ft) == 1:
            type_ft = '00{}'.format(self.type_ft)

        if file_type == 'ft':
            if custom_file:
                if self.dir_custom_ft is None:
                    print 'No custom ft directory selected'
                    return None

                file_ft = glob(os.path.join(self.dir_custom_ft, '*.{}.*'.format(type_ft)))

                if len(file_ft) == 1:
                    return file_ft[0]

                else:
                    print 'Multiple ft files selected: {}'.format(file_ft)
                    return None

            else:
                if type_ft == '000':
                    return self.ft0
                if type_ft == '002':
                    return self.ft2
                if type_ft == '004':
                    return self.ft4
                if type_ft == '006':
                    return self.ft6

        if file_type == 'rmsd':
            if self.dir_custom_rmsd is None:
                print 'No custom rmsd directory selected'
                return None

            file_rmsd = glob(os.path.join(self.dir_custom_rmsd, '*.{}.*'.format(type_ft)))

            if len(file_rmsd) == 1:
                return file_rmsd[0]

            else:
                print 'Multiple rmsd files selected: {}'.format(file_rmsd)
                return None

        if file_type == 'pwrmsd':
            if self.dir_custom_pwrmsd is None:
                print 'No custom pwrmsd directory selected'
                return None

            file_pwrmsd = glob(os.path.join(self.dir_custom_pwrmsd, '*.{}.*'.format(type_ft)))

            if len(file_pwrmsd) == 1:
                return file_pwrmsd[0]

            else:
                print 'Multiple pwrmsd files selected: {}'.format(file_pwrmsd)
                return None

    def rmsd(self, custom_ft=False, custom_lig=False):

        if self.type_ft is None and self.custom_ft is None:
            print 'No ft type or ft file was selected'
            return None

        if custom_ft and self.custom_ft is None:
            file_ft = self.id_file('ft', custom_file=True)

        elif custom_ft:
            file_ft = self.id_file('ft')

        if custom_ft:
            if self.dir_custom_ft is None:

        _rmsd(self.lig, self.crys_lig, file_ft, self.rot)



