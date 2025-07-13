import numpy as np
from typing import Tuple, Optional, Dict
from collections import defaultdict
from aligner.algorithms import _get_pair_score, smith_waterman, needleman_wunsch


def hirschberg_needleman_wunsch(
        seq1: str,
        seq2: str,
        match_score: int = 1,
        mismatch_score: int = -1,
        gap_penalty: int = -2,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None
) -> Tuple[str, str, int]:
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n, m = len(seq1), len(seq2)

    if n == 0:
        return '-' * m, seq2, gap_penalty * m
    if m == 0:
        return seq1, '-' * n, gap_penalty * n
    if n == 1 or m == 1:
        return needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty, scoring_matrix)

    mid = n // 2
    score_left = _compute_nw_row(seq1[:mid], seq2, match_score, mismatch_score, gap_penalty, scoring_matrix)
    score_right = _compute_nw_row(seq1[mid:][::-1], seq2[::-1], match_score, mismatch_score, gap_penalty,
                                  scoring_matrix)[::-1]

    split = np.argmax(score_left + score_right)

    align1_left, align2_left, _ = hirschberg_needleman_wunsch(seq1[:mid], seq2[:split], match_score, mismatch_score,
                                                              gap_penalty, scoring_matrix)
    align1_right, align2_right, _ = hirschberg_needleman_wunsch(seq1[mid:], seq2[split:], match_score, mismatch_score,
                                                                gap_penalty, scoring_matrix)

    full_align1 = align1_left + align1_right
    full_align2 = align2_left + align2_right
    full_score = score_left[split] + score_right[split]

    return full_align1, full_align2, full_score


def _compute_nw_row(
        seq1: str,
        seq2: str,
        match_score: int,
        mismatch_score: int,
        gap_penalty: int,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]]
) -> np.ndarray:
    m = len(seq2)
    prev = np.arange(0, gap_penalty * (m + 1), gap_penalty)
    for char1 in seq1:
        curr = np.zeros(m + 1, dtype=int)
        curr[0] = prev[0] + gap_penalty
        for j in range(1, m + 1):
            pair_score = _get_pair_score(char1, seq2[j - 1], match_score, mismatch_score, scoring_matrix)
            curr[j] = max(
                prev[j - 1] + pair_score,
                prev[j] + gap_penalty,
                curr[j - 1] + gap_penalty
            )
        prev = curr
    return prev


def heuristic_local_align(
        seq1: str,
        seq2: str,
        k: int = 3,
        match_score: int = 1,
        mismatch_score: int = -1,
        gap_penalty: int = -2,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None
) -> Tuple[str, str, int]:
    index = defaultdict(list)
    for i in range(len(seq2) - k + 1):
        kmer = seq2[i:i + k]
        index[kmer].append(i)

    best_score = 0
    best_align = ("", "", 0)
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i:i + k]
        for j in index.get(kmer, []):
            left1 = seq1[:i]
            left2 = seq2[:j]
            right1 = seq1[i + k:]
            right2 = seq2[j + k:]

            _, _, left_score = smith_waterman(left1[::-1], left2[::-1], match_score, mismatch_score, gap_penalty,
                                              scoring_matrix)
            _, _, right_score = smith_waterman(right1, right2, match_score, mismatch_score, gap_penalty, scoring_matrix)
            hit_score = k * match_score + left_score + right_score

            if hit_score > best_score:
                full_align1, full_align2, full_score = smith_waterman(seq1[max(0, i - 50):i + k + 50],
                                                                      seq2[max(0, j - 50):j + k + 50], match_score,
                                                                      mismatch_score, gap_penalty, scoring_matrix)
                if full_score > best_score:
                    best_score = full_score
                    best_align = (full_align1, full_align2, full_score)

    return best_align