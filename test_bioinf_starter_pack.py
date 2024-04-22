import os

import pytest

from bio_files_processor import OpenFasta
from bioinf_starter_pack import (AminoAcidSequence,
                                 DNASequence,
                                 RNASequence,
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


class TestGroupOpenFasta:
    @pytest.fixture
    def input_file_path(self):
        file_path = 'data/example_fasta.fasta'
        return file_path

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


class TestGroupFilterFastq:
    @pytest.fixture
    def input_file_path(self):
        file_path = 'data/example_fastq.fastq'
        return file_path

    @pytest.fixture
    def tmp_file(self):
        file_name = 'tmp.fastq'
        yield file_name
        file_path = f'fastq_filtrator_results/{file_name}'
        if os.path.exists(file_path):
            os.remove(file_path)
        if os.listdir('fastq_filtrator_results') == []:
            os.rmdir('fastq_filtrator_results')

    def test_output_file_exists(self, input_file_path, tmp_file):
        filter_fastq(input_file_path, tmp_file)
        assert os.path.exists(f'fastq_filtrator_results/{tmp_file}')
