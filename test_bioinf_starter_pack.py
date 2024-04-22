import pytest

from bioinf_starter_pack import DNASequence, AminoAcidSequence, RNASequence


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
