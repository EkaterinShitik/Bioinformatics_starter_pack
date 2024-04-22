import os

import pytest

from bio_files_processor import OpenFasta, convert_multiline_fasta_to_oneline
from bioinf_starter_pack import (AminoAcidSequence, DNASequence, RNASequence,
                                 filter_fastq)


class TestGroupBiolSequence:
    @pytest.fixture
    def input_data(self):
        dna_seq = 'ATGC'
        return dna_seq

    def test_transcribed_dna_type(self, input_data):
        inp = input_data
        target_type = RNASequence
        result = DNASequence(inp).transcribe()
        assert isinstance(result, target_type)

    def test_rna_content(self, input_data):
        inp = input_data
        with pytest.raises(ValueError):
            RNASequence(inp)

    def test_prot_search4alt_frames(self):
        inp = 'MYRHHWWMYYYYYYYM'
        target = 'MYYYYYYYM'
        sample = AminoAcidSequence(inp)
        result = sample.search_for_alt_frames()[0].seq
        assert target == result


class InputFastaFile:
    @pytest.fixture
    def input_file_path(self):
        file_path = 'data/example_fasta.fasta'
        return file_path


class TestGroupOpenFasta(InputFastaFile):
    def test_number_records(self, input_file_path):
        inp = input_file_path
        target_number = 5
        with OpenFasta(inp, 'r') as fasta_file:
            result_number = len(fasta_file.read_records())
        assert target_number == result_number

    def test_last_record_id(self, input_file_path):
        inp = input_file_path
        target_id = 'GTD129563'
        with OpenFasta(inp, 'r') as fasta_file:
            last_record = fasta_file.read_records()[-1]
            result_id = last_record.id
        assert target_id == result_id


class TestGroupConvMultFasta(InputFastaFile):
    def test_rewrite_files(self, input_file_path):
        output_path = input_file_path.split('.fasta')[0]
        with pytest.raises(ValueError):
            convert_multiline_fasta_to_oneline(input_file_path, output_path)


class TestGroupFilterFastq:
    @pytest.fixture
    def input_file_path(self):
        file_path = 'data/example_fastq.fastq'
        return file_path

    @pytest.fixture
    def tmp_file(self):
        file_name = 'tmp.fastq'
        yield file_name
        obligatory_dir = 'fastq_filtrator_results'
        file_path = os.path.join(obligatory_dir, file_name)
        if os.path.exists(file_path):
            os.remove(file_path)
        if not os.listdir(obligatory_dir):
            os.rmdir(obligatory_dir)

    def test_output_file_exists(self, input_file_path, tmp_file):
        filter_fastq(input_file_path, tmp_file)
        target_file_path = os.path.join('fastq_filtrator_results', tmp_file)
        assert os.path.exists(target_file_path)

    def test_quality_threshold_output(self, input_file_path, tmp_file):
        target_seqs = ['GACCTTTCCGCAAGCTGTCGC', 'CATGGTGGCG', 'C']
        filter_fastq(input_file_path, tmp_file, quality_threshold=38)
        file_path = os.path.join('fastq_filtrator_results', tmp_file)
        result_seqs = []
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('@'):
                    seq = file.readline().strip()
                    result_seqs.append(seq)
        assert target_seqs == result_seqs
