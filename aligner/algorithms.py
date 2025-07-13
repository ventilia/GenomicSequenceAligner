import numpy as np
from typing import Tuple, Optional, Dict


def _get_pair_score(
        char1: str,
        char2: str,
        match_score: int,
        mismatch_score: int,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]]
) -> int:

    if scoring_matrix:
        return scoring_matrix.get((char1.upper(), char2.upper()), mismatch_score)
    return match_score if char1 == char2 else mismatch_score


def _backtrace(
        dp: np.ndarray,
        seq1: str,
        seq2: str,
        start_i: int,
        start_j: int,
        gap_penalty: int,
        pair_score_func,
        stop_condition
) -> Tuple[str, str]:

    """
    backtrace для выравнивания

     dp: матрица dynamic programming
    seq1: 1 последовательность.
   seq2: 2 последовательность
     start_i: start_j: начальный индекс i и j для backtrace

     gap_penalty: штраф за gap
     pair_score_func: функция для вычисления pair_score
     stop_condition: функция, определяющая, когда остановить backtrace

    """
    align1, align2 = [], []
    i, j = start_i, start_j
    while not stop_condition(i, j, dp):
        if i > 0 and j > 0:
            pair_score = pair_score_func(seq1[i - 1], seq2[j - 1])
            if dp[i][j] == dp[i - 1][j - 1] + pair_score:
                align1.append(seq1[i - 1])
                align2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
            align1.append(seq1[i - 1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j - 1])
            j -= 1
    return ''.join(reversed(align1)), ''.join(reversed(align2))


def needleman_wunsch(
        seq1: str,
        seq2: str,
        match_score: int = 1,
        mismatch_score: int = -1,
        gap_penalty: int = -2,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None
) -> Tuple[str, str, int]:
    """
    Needleman-Wunsch

     seq1: 1 последовательность (дНК/белок)
     seq2: 2 последовательность.
     match_score: score за совпадение (если нет матрицы)
    mismatch_score: score за несовпадение
     gap_penalty: шраф за gap
     scoring_matrix: рпциональная матрица scoring

    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n, m = len(seq1), len(seq2)
    if n > 10000 or m > 10000:
        raise ValueError("последовательности слишком длинные для базовой версии. используйте оптимизированную")

    dp = np.zeros((n + 1, m + 1), dtype=int)


    for i in range(1, n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(1, m + 1):
        dp[0][j] = j * gap_penalty

    #  матрица  судьбы
    pair_score_func = lambda c1, c2: _get_pair_score(c1, c2, match_score, mismatch_score, scoring_matrix)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair_score = pair_score_func(seq1[i - 1], seq2[j - 1])
            dp[i][j] = max(
                dp[i - 1][j - 1] + pair_score,
                dp[i - 1][j] + gap_penalty,
                dp[i][j - 1] + gap_penalty
            )

    # backtrace
    stop_condition = lambda i, j, dp: i == 0 and j == 0
    align1, align2 = _backtrace(dp, seq1, seq2, n, m, gap_penalty, pair_score_func, stop_condition)
    return align1, align2, dp[n][m]


def smith_waterman(
        seq1: str,
        seq2: str,
        match_score: int = 1,
        mismatch_score: int = -1,
        gap_penalty: int = -2,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None
) -> Tuple[str, str, int]:
    """
 Smith-Waterman для локального

    seq1: первая последовательность
     seq2: вторая последовательность
     match_score: score за совпадение
     mismatch_score: score за несовпадение
    gap_penalty: iтраф за gap
     scoring_matrix: опциональная матрица.

    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n, m = len(seq1), len(seq2)
    if n > 10000 or m > 10000:
        raise ValueError("последовательности  длинные для этой версии используйте оптимизированную.")

    dp = np.zeros((n + 1, m + 1), dtype=int)

    #  матрица
    pair_score_func = lambda c1, c2: _get_pair_score(c1, c2, match_score, mismatch_score, scoring_matrix)
    max_score = 0
    max_i, max_j = 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair_score = pair_score_func(seq1[i - 1], seq2[j - 1])
            dp[i][j] = max(
                0,
                dp[i - 1][j - 1] + pair_score,
                dp[i - 1][j] + gap_penalty,
                dp[i][j - 1] + gap_penalty
            )
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_i, max_j = i, j

    # Backtrace
    stop_condition = lambda i, j, dp: i == 0 or j == 0 or dp[i][j] <= 0
    align1, align2 = _backtrace(dp, seq1, seq2, max_i, max_j, gap_penalty, pair_score_func, stop_condition)
    return align1, align2, max_score