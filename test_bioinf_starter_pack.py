import pytest

from bioinf_starter_pack import DNASequence, AminoAcidSequence


def test_dna_type():
    inp = 'ATGC'
    target_type = DNASequence
    result = DNASequence(inp).complement()
    assert isinstance(result, target_type)


def test_protein_content():
    inp = 'BGARB'
    with pytest.raises(ValueError):
        AminoAcidSequence(inp)
