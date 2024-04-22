import pytest

from bioinf_starter_pack import DNASequence, AminoAcidSequence, RNASequence


class TestGroupBiolSequence:
    @pytest.fixture
    def input_data(self):
        inp = 'ATGC'
        return inp

    def test_transcribed_dna_type(self, input_data):
        inp = input_data
        target_type = RNASequence
        result = DNASequence(inp).transcribe()
        assert isinstance(result, target_type)

    def test_protein_content(self):
        inp = 'BGARB'
        with pytest.raises(ValueError):
            AminoAcidSequence(inp)
