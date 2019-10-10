import os
import shutil

class CompressedSession:

    def __init__(self, raw_path, crys_lig=None):

        # Raw Files:
        self.raw_path = raw_path
        self.ft0_raw = os.path.join(self.raw_path, 'ft.000.00.gz')
        self.ft2_raw = os.path.join(self.raw_path, 'ft.002.00.gz')
        self.ft4_raw = os.path.join(self.raw_path, 'ft.004.00.gz')
        self.ft6_raw = os.path.join(self.raw_path, 'ft.006.00.gz')
        self.rot_raw = os.path.join(self.raw_path, 'prms/rot70k.0.0.4.prm')
        self.lig_raw = os.path.join(self.raw_path, 'lig.pdb.gz')
        self.rec_raw = os.path.join(self.raw_path, 'rec.pdb.gz')

        if crys_lig is None:
            self.crys_lig = self.lig_raw

    def convert(self, session_path, contents=None):
        """
        Converts a raw session comprising a raw ClusPro output into a "working session"
        :param session_path: path of destination for working session directory
        :param contents: list of contents to convert, default = None
        """
        DEFAULT_CONTENTS = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                            self.ft6_raw, self.lig_raw, self.rec_raw, self.rot_raw]

        EXTRACTABLE_CONTENTS = [self.ft0_raw, self.ft2_raw, self.ft4_raw,
                                self.ft6_raw, self.lig_raw, self.rec_raw]
        if contents is None:
            contents = DEFAULT_CONTENTS

        for item in contents:

            if item in EXTRACTABLE_CONTENTS:
                extract(item, session_path)
            elif item == self.rot_raw:
                shultil.copy(self.rot_raw, session_path)