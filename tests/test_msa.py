import pytest
from aligner.msa import multiple_sequence_alignment, MSAError
from aligner.scoring import load_scoring_matrix

def test_multiple_sequence_alignment_basic():
    seqs = ["AGC", "ACGC", "AGGC"]
    aligned = multiple_sequence_alignment(seqs)
    assert len(aligned) == 3
    assert all(len(a) == len(aligned[0]) for a in aligned)
    assert "A-GC" in aligned or "AGC-" in aligned

def test_multiple_sequence_alignment_with_matrix():
    seqs = ["ILK", "IMK", "ILR"]
    matrix = load_scoring_matrix("BLOSUM62")
    aligned = multiple_sequence_alignment(seqs, scoring_matrix=matrix)
    assert len(aligned) == 3
    assert all(len(a) == 3 for a in aligned)

@pytest.mark.parametrize("seqs", [
    (["A"]),
     ([]),
])
def test_msa_error(seqs):
    with pytest.raises(MSAError):
        multiple_sequence_alignment(seqs)

def test_msa_identical():
    seqs = ["AAA", "AAA", "AAA"]
    aligned = multiple_sequence_alignment(seqs)
    assert all(a == "AAA" for a in aligned)