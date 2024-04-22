import pytest

from bioinf_starter_pack import DNASequence, AminoAcidSequence, RNASequence
from bio_files_processor import OpenFasta


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
