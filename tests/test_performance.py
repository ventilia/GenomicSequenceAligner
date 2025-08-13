import pytest
import os
import sys
import timeit
from aligner.algorithms import needleman_wunsch, smith_waterman

def test_cli_benchmark():

    seq1 = "A" * 1000
    seq2 = "A" * 1000
    time_full = timeit.timeit(lambda: needleman_wunsch(seq1, seq2), number=1)
    time_banded = timeit.timeit(lambda: needleman_wunsch(seq1, seq2, bandwidth=100), number=1)
    assert time_banded < time_full * 1.5

def test_smith_waterman_performance():
    """Тест производительности локального выравнивания"""
    seq1 = "ATGC" * 100
    seq2 = "TGCA" * 100
    time_taken = timeit.timeit(lambda: smith_waterman(seq1, seq2), number=10)
    assert time_taken < 5  # Должно выполниться менее чем за 5 секунд

def test_msa_performance():
    """Тест производительности множественного выравнивания"""
    from aligner.msa import multiple_sequence_alignment
    seqs = ["AGCT" * 50, "ACGT" * 50, "AGGT" * 50, "ATGT" * 50]
    time_taken = timeit.timeit(lambda: multiple_sequence_alignment(seqs), number=1)
    assert time_taken < 10  # Должно выполниться менее чем за 10 секунд