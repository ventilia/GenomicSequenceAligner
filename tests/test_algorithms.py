import pytest
from aligner.algorithms import needleman_wunsch, smith_waterman
from aligner.scoring import load_scoring_matrix

@pytest.mark.parametrize("seq1, seq2, expected_align1, expected_align2, expected_score", [
    ("AGC", "ACGC", "A-GC", "ACGC", 1),
    ("AA", "AA", "AA", "AA", 2),
    ("A", "", "A", "-", -2),
])
def test_needleman_wunsch(seq1, seq2, expected_align1, expected_align2, expected_score):
    align1, align2, score = needleman_wunsch(seq1, seq2)
    assert align1 == expected_align1
    assert align2 == expected_align2
    assert score == expected_score

@pytest.mark.parametrize("seq1, seq2, expected_align1, expected_align2, expected_score", [
    ("AGC", "ACGC", "GC", "GC", 2),
    ("ATGC", "TGCA", "TGC", "TGC", 3),
])
def test_smith_waterman(seq1, seq2, expected_align1, expected_align2, expected_score):
    align1, align2, score = smith_waterman(seq1, seq2)
    assert align1 == expected_align1
    assert align2 == expected_align2
    assert score == expected_score

def test_needleman_wunsch_with_matrix():
    matrix = load_scoring_matrix("BLOSUM62")
    align1, align2, score = needleman_wunsch("IL", "IM", match_score=0, mismatch_score=0, gap_penalty=0, scoring_matrix=matrix)
    assert score > 0