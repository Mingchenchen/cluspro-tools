import os
import tempfile
import shutil
import click

import Bio.Seq
import Bio.SeqRecord
import Bio.SeqIO
import Bio.Alphabet


class ModellerError(Exception):
    pass


import modeller
import modeller.automodel


class fixedautomodel(modeller.automodel.automodel):
    def select_atoms(self):
        aln = self.read_alignment()
        l = self.loops(aln, minlength=1, maxlength=len(aln.positions), insertion_ext=0, deletion_ext=0)
        return modeller.selection(l) | modeller.selection(self.atoms).only_sidechain()


def build_model(model_pdb_file, model_seq, ref_pdb_file, ref_chain):
    ref_pdb_file_abs = os.path.abspath(ref_pdb_file)
    ref_pdb_file_abs_dir = os.path.dirname(ref_pdb_file_abs)
    model_pdb_file_abs = os.path.abspath(model_pdb_file)

    current_dir = os.getcwd()
    tmp_dir = tempfile.mkdtemp()

    os.chdir(tmp_dir)
    try:
        model_code = 'UKNP'
        ref_code = '{}{}'.format(os.path.splitext(os.path.basename(ref_pdb_file_abs))[0], ref_chain)

        model_seq_pir_file = 'target_seq.pir'
        with open(model_seq_pir_file, 'w') as f:
            model_seq_obj = Bio.Seq.Seq(model_seq, Bio.Alphabet.generic_protein)
            model_seq_record = Bio.SeqRecord.SeqRecord(model_seq_obj, id=model_code, description='', name='sequence:{}::A::A::::'.format(model_code))
            Bio.SeqIO.write(model_seq_record, f, 'pir')

        modeller.log.minimal()  # do NOT use log.none() - will suppress some errors
        env = modeller.environ()
        env.io.atom_files_directory = [ref_pdb_file_abs_dir]
        aln = modeller.alignment(env)
        mdl = modeller.model(env, file=ref_pdb_file_abs, model_segment=('FIRST:{}'.format(ref_chain), 'LAST:{}'.format(ref_chain)))
        aln.append_model(mdl, atom_files=ref_pdb_file_abs, align_codes=ref_code)
        aln.append(file=model_seq_pir_file, align_codes=model_code)
        aln.align2d()

        alignment_pir_file = 'alignment.pir'
        aln.write(file=alignment_pir_file, alignment_format='PIR')

        a = fixedautomodel(env=env, alnfile=alignment_pir_file, knowns=ref_code, sequence=model_code)
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        shutil.move('{}.B99990001.pdb'.format(model_code), model_pdb_file_abs)
    except:
        raise ModellerError('Modeller error')
    finally:
        os.chdir(current_dir)
        shutil.rmtree(tmp_dir)


@click.command()
@click.argument('model_pdb_file', type=click.Path())
@click.argument('model_fasta_file', type=click.Path(exists=True))
@click.argument('ref_pdb_file', type=click.Path(exists=True))
@click.argument('ref_chain')
def _build_model_click(model_pdb_file, model_fasta_file, ref_pdb_file, ref_chain):
    '''
        Modeller-based script for homology modelling:\n
        1) Creates Modeller alignment between the chain REF_CHAIN of the reference PDB structure REF_PDB_FILE and the
        sequence MODEL_SEQ.\n
        2) Build structure model MODEL_PDB_FILE of the sequence MODEL_SEQ based on the alignment from 1).\n
        Note: Models only the unaligned portion of the sequence, the aligned portions are built with exactly the same
        backbone as the template, and the side chain positions are placed according to some unknown (by me) algorithm.
    '''
    model_seq_record = next(Bio.SeqIO.parse(model_fasta_file, 'fasta'))
    model_seq = str(model_seq_record.seq)
    build_model(model_pdb_file, model_seq, ref_pdb_file, ref_chain)


if __name__ == '__main__':
    _build_model_click()