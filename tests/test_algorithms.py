import pytest
from unittest.mock import patch, MagicMock
from aligner.algorithms import needleman_wunsch, smith_waterman
from aligner.scoring import load_scoring_matrix
from aligner.msa import multiple_sequence_alignment, MSAError
from subprocess import run, CalledProcessError
import os
import sys
import timeit


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


@pytest.mark.parametrize("seq1, seq2, bandwidth, expected_score", [
    ("AGCTAGCT", "AGCTGGCT", 4, 6),
    ("AGCTAGCT", "AGCTGGCT", None, 6),
])
def test_needleman_wunsch_banded(seq1, seq2, bandwidth, expected_score):
    _, _, score = needleman_wunsch(seq1, seq2, bandwidth=bandwidth)
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


def test_needleman_wunsch_affine():
    align1, align2, score = needleman_wunsch("AGCT", "A--T", gap_penalty=-2, gap_open=-5, gap_extend=-1)
    assert score == 0


def test_multiple_sequence_alignment_basic():
    seqs = ["AGC", "ACGC", "AGGC"]
    aligned = multiple_sequence_alignment(seqs)
    assert len(aligned) == 3
    assert all(len(a) == len(aligned[0]) for a in aligned)
    assert "A-GC" in aligned or "AGC-" in aligned


def test_msa_parallel():
    seqs = ["AGC", "ACGC", "AGGC"]
    with patch('aligner.msa.Pool') as mock_pool:
        mock_instance = MagicMock()
        mock_pool.return_value.__enter__.return_value = mock_instance

        # Side effect для map: вычисляем serially, чтобы dist не zero
        from aligner.msa import pairwise_distance  # Импорт для real func
        def real_map(f, tasks):
            return [pairwise_distance(arg) for arg in tasks]

        mock_instance.map.side_effect = real_map

        multiple_sequence_alignment(seqs, threads=2)
        mock_pool.assert_called_once_with(2)


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


@pytest.fixture
def dummy_fasta():
    with open("test1.fasta", "w") as f:
        f.write(">seq1\nAGC")
    with open("test2.fasta", "w") as f:
        f.write(">seq2\nACGC")
    with open("test_multi.fasta", "w") as f:
        f.write(">seq1\nAGC\n>seq2\nACGC\n>seq3\nAGGC")
    yield
    os.remove("test1.fasta")
    os.remove("test2.fasta")
    os.remove("test_multi.fasta")


def test_cli_pairwise(dummy_fasta):
    try:
        run([sys.executable, "-m", "aligner.cli", "global", "--input1", "test1.fasta", "--input2", "test2.fasta", "--output", "test_out.txt", "--lang", "ru"], check=True)
        assert os.path.exists("test_out.txt")
        with open("test_out.txt", "r") as f:
            content = f.read()
            assert "Alignment Score" in content or "Счет выравнивания" in content  # check lang
        os.remove("test_out.txt")
    except CalledProcessError as e:
        pytest.fail(f"CLI failed: {e}")


def test_cli_msa(dummy_fasta):
    try:
        run([sys.executable, "-m", "aligner.cli", "msa", "--input1", "test_multi.fasta", "--output", "test_msa.txt", "--lang", "en"], check=True)
        assert os.path.exists("test_msa.txt")
        with open("test_msa.txt", "r") as f:
            content = f.read()
            assert "Seq1" in content
        os.remove("test_msa.txt")
    except CalledProcessError as e:
        pytest.fail(f"CLI failed: {e}")


def test_cli_benchmark():

    seq1 = "A" * 1000
    seq2 = "A" * 1000
    time_full = timeit.timeit(lambda: needleman_wunsch(seq1, seq2), number=1)
    time_banded = timeit.timeit(lambda: needleman_wunsch(seq1, seq2, bandwidth=100), number=1)
    assert time_banded < time_full * 1.5