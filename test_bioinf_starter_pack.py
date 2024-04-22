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

    def test_prot_search4alt_frames(self):
        inp = 'MYRHHWWMYYYYYYYM'
        target = 'MYYYYYYYM'
        sample = AminoAcidSequence(inp)
        result = sample.search_for_alt_frames()[0].seq
        assert target == result
