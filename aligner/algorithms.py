import numpy as np
from typing import Tuple, Optional, Dict
import numba


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


def _backtrace_linear(
        dp: np.ndarray,
        seq1: str,
        seq2: str,
        start_i: int,
        start_j: int,
        gap_penalty: int,
        pair_score_func,
        stop_condition,
        bandwidth: Optional[int] = None
) -> Tuple[str, str]:
    align1, align2 = [], []
    i, j = start_i, start_j
    if bandwidth is not None:
        width = 2 * bandwidth + 1
        while not stop_condition(i, j, dp):
            col = j - (i - bandwidth)
            if col < 0 or col >= width:
                # Выход за band
                if col < 0:
                    align1.append(seq1[i - 1])
                    align2.append('-')
                    i -= 1
                else:
                    align1.append('-')
                    align2.append(seq2[j - 1])
                    j -= 1
                continue
            if i > 0 and j > 0 and col > 0:  # col-1 >=0
                pair_score = pair_score_func(seq1[i - 1], seq2[j - 1])
                if dp[i][col] == dp[i - 1][col - 1] + pair_score:
                    align1.append(seq1[i - 1])
                    align2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                    continue
            if i > 0 and dp[i][col] == dp[i - 1][col] + gap_penalty:
                align1.append(seq1[i - 1])
                align2.append('-')
                i -= 1
                continue
            if j > 0:
                align1.append('-')
                align2.append(seq2[j - 1])
                j -= 1
    else:
        while not stop_condition(i, j, dp):
            if i > 0 and j > 0:
                pair_score = pair_score_func(seq1[i - 1], seq2[j - 1])
                if dp[i][j] == dp[i - 1][j - 1] + pair_score:
                    align1.append(seq1[i - 1])
                    align2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                    continue
            if i > 0 and (j == 0 or dp[i][j] == dp[i - 1][j] + gap_penalty):
                align1.append(seq1[i - 1])
                align2.append('-')
                i -= 1
                continue
            if j > 0:
                align1.append('-')
                align2.append(seq2[j - 1])
                j -= 1
    while i > 0:
        align1.append(seq1[i - 1])
        align2.append('-')
        i -= 1
    while j > 0:
        align1.append('-')
        align2.append(seq2[j - 1])
        j -= 1
    return ''.join(reversed(align1)), ''.join(reversed(align2))


def _backtrace_affine(
        M: np.ndarray,
        Ix: np.ndarray,
        Iy: np.ndarray,
        seq1: str,
        seq2: str,
        start_i: int,
        start_j: int,
        gap_open: int,
        gap_extend: int,
        pair_score_func,
        stop_condition
) -> Tuple[str, str]:
    align1, align2 = [], []
    i, j = start_i, start_j
    current_matrix = np.argmax([M[i][j], Ix[i][j], Iy[i][j]])  # 0: M, 1: Ix, 2: Iy
    while i > 0 or j > 0:
        if current_matrix == 0:  # M
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1
            current_matrix = np.argmax([M[i][j], Ix[i][j], Iy[i][j]])
        elif current_matrix == 1:  # Ix (gap in seq2)
            align1.append(seq1[i - 1])
            align2.append('-')
            i -= 1
            if i >= 0 and Ix[i + 1][j] == M[i][j] + gap_open:
                current_matrix = 0
            else:
                current_matrix = 1
        else:  # Iy (gap in seq1)
            align1.append('-')
            align2.append(seq2[j - 1])
            j -= 1
            if j >= 0 and Iy[i][j + 1] == M[i][j] + gap_open:
                current_matrix = 0
            else:
                current_matrix = 2
    return ''.join(reversed(align1)), ''.join(reversed(align2))


@numba.jit(nopython=True)
def _compute_nw_row_vectorized(
        seq1: str,
        seq2: str,
        match_score: int,
        mismatch_score: int,
        gap_penalty: int
) -> np.ndarray:
    # векторизованная версия заполнения строки NW (для Hirschberg)
    n = len(seq1)
    m = len(seq2)
    prev = np.arange(0, gap_penalty * (m + 1), gap_penalty, dtype=np.float64)
    for i in range(1, n + 1):
        curr = np.zeros(m + 1, dtype=np.float64)
        curr[0] = prev[0] + gap_penalty
        char1 = seq1[i - 1]
        pair_scores = np.array([match_score if char1 == seq2[j - 1] else mismatch_score for j in range(1, m + 1)], dtype=np.float64)
        diag = prev[:-1] + pair_scores
        up = prev[1:] + gap_penalty
        left = curr[:-1] + gap_penalty
        curr[1:] = np.maximum.reduce([diag, up, left])
        prev = curr
    return prev


def needleman_wunsch(
        seq1: str,
        seq2: str,
        match_score: int = 1,
        mismatch_score: int = -1,
        gap_penalty: int = -2,
        gap_open: Optional[int] = None,
        gap_extend: Optional[int] = None,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None,
        bandwidth: Optional[int] = None
) -> Tuple[str, str, int]:
    # Needleman-Wunsch с векторизацией, banded и affine gaps
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n, m = len(seq1), len(seq2)
    if n > 10000 or m > 10000:
        raise ValueError("Последовательности слишком длинные для базовой версии. Используйте оптимизированную")

    affine = gap_open is not None and gap_extend is not None
    if affine and bandwidth is not None:
        raise ValueError("Banded not supported for affine gaps yet")
    pair_score_func = lambda c1, c2: _get_pair_score(c1, c2, match_score, mismatch_score, scoring_matrix)

    if affine:
        # Affine gaps (full DP only)
        M = np.zeros((n + 1, m + 1), dtype=np.float64)
        Ix = np.full((n + 1, m + 1), -np.inf, dtype=np.float64)
        Iy = np.full((n + 1, m + 1), -np.inf, dtype=np.float64)
        for i in range(1, n + 1):
            Ix[i, 0] = gap_open + (i - 1) * gap_extend
            M[i, 0] = Ix[i, 0]
        for j in range(1, m + 1):
            Iy[0, j] = gap_open + (j - 1) * gap_extend
            M[0, j] = Iy[0, j]
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                pair = pair_score_func(seq1[i-1], seq2[j-1])
                M[i, j] = pair + max(M[i-1, j-1], Ix[i-1, j-1], Iy[i-1, j-1])
                Ix[i, j] = max(M[i-1, j] + gap_open, Ix[i-1, j] + gap_extend)
                Iy[i, j] = max(M[i, j-1] + gap_open, Iy[i, j-1] + gap_extend)
        stop_condition = lambda i, j, M: i == 0 or j == 0
        align1, align2 = _backtrace_affine(M, Ix, Iy, seq1, seq2, n, m, gap_open, gap_extend, pair_score_func, stop_condition)
        score = max(M[n, m], Ix[n, m], Iy[n, m])
    else:
        if bandwidth is not None:
            width = 2 * bandwidth + 1
            dp = np.full((n + 1, width), -np.inf, dtype=np.float64)

            for j in range(m + 1):
                col = j + bandwidth  # col = j - (0 - bandwidth)
                if 0 <= col < width:
                    dp[0, col] = j * gap_penalty
            for i in range(1, n + 1):
                # left column (gaps в seq2)
                col_left = bandwidth - i  # col для j=0: 0 - (i - bandwidth) = bandwidth - i
                if 0 <= col_left < width:
                    dp[i, col_left] = i * gap_penalty
                for col in range(width):
                    j = col + (i - bandwidth)
                    if j < 0 or j > m:
                        continue
                    if j == 0:
                        continue
                    pair = pair_score_func(seq1[i - 1], seq2[j - 1])

                    diag = (dp[i - 1, col] + pair) if j > 0 else -np.inf
                    up = (dp[i - 1, col + 1] + gap_penalty) if col + 1 < width else -np.inf
                    left = (dp[i, col - 1] + gap_penalty) if col - 1 >= 0 else -np.inf
                    dp[i, col] = max(diag, up, left)
            stop_condition = lambda i, j, dp: i == 0 or j == 0
            align1, align2 = _backtrace_linear(dp, seq1, seq2, n, m, gap_penalty, pair_score_func, stop_condition, bandwidth)
            col = m - (n - bandwidth)
            score = dp[n, col] if 0 <= col < width else np.max(dp[n])
        else:
            dp = np.zeros((n + 1, m + 1), dtype=np.float64)
            if gap_penalty != 0:
                dp[:, 0] = np.arange(0, gap_penalty * (n + 1), gap_penalty)
                dp[0, :] = np.arange(0, gap_penalty * (m + 1), gap_penalty)
            else:
                dp[:, 0] = 0
                dp[0, :] = 0
            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    pair = pair_score_func(seq1[i - 1], seq2[j - 1])
                    dp[i][j] = max(
                        dp[i - 1][j - 1] + pair,
                        dp[i - 1][j] + gap_penalty,
                        dp[i][j - 1] + gap_penalty
                    )
            stop_condition = lambda i, j, dp: i == 0 or j == 0
            align1, align2 = _backtrace_linear(dp, seq1, seq2, n, m, gap_penalty, pair_score_func, stop_condition, bandwidth)
            score = dp[n][m]
    if np.isinf(score) or np.isnan(score):
        score = -10000000000
    return align1, align2, int(score)


def smith_waterman(
        seq1: str,
        seq2: str,
        match_score: int = 1,
        mismatch_score: int = -1,
        gap_penalty: int = -2,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None,
        bandwidth: Optional[int] = None
) -> Tuple[str, str, int]:
    # Smith-Waterman
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n, m = len(seq1), len(seq2)
    if n > 10000 or m > 10000:
        raise ValueError("Последовательности слишком длинные для этой версии. Используйте оптимизированную.")

    dp = np.zeros((n + 1, m + 1), dtype=np.float64)
    pair_score_func = lambda c1, c2: _get_pair_score(c1, c2, match_score, mismatch_score, scoring_matrix)
    max_score = 0.0
    max_i, max_j = 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = pair_score_func(seq1[i - 1], seq2[j - 1])
            dp[i][j] = max(
                0.0,
                dp[i - 1][j - 1] + pair,
                dp[i - 1][j] + gap_penalty,
                dp[i][j - 1] + gap_penalty
            )
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_i, max_j = i, j

    stop_condition = lambda i, j, dp: i == 0 or j == 0 or dp[i][j] == 0
    align1, align2 = _backtrace_linear(dp, seq1, seq2, max_i, max_j, gap_penalty, pair_score_func, stop_condition, bandwidth)
    min_len = min(len(align1), len(align2))
    start = 0
    while start < min_len and (align1[start] == '-' or align2[start] == '-'):
        start += 1
    end = min_len
    while end > start and (align1[end-1] == '-' or align2[end-1] == '-'):
        end -= 1
    align1 = align1[start:end]
    align2 = align2[start:end]
    return align1, align2, int(max_score)


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
        return needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty, None, None, scoring_matrix)

    mid = n // 2
    score_left = _compute_nw_row_vectorized(seq1[:mid], seq2, match_score, mismatch_score, gap_penalty)
    score_right = _compute_nw_row_vectorized(seq1[mid:][::-1], seq2[::-1], match_score, mismatch_score, gap_penalty)[::-1]

    split = np.argmax(score_left + score_right)

    align1_left, align2_left, _ = hirschberg_needleman_wunsch(seq1[:mid], seq2[:split], match_score, mismatch_score, gap_penalty, scoring_matrix)
    align1_right, align2_right, _ = hirschberg_needleman_wunsch(seq1[mid:], seq2[split:], match_score, mismatch_score, gap_penalty, scoring_matrix)

    full_align1 = align1_left + align1_right
    full_align2 = align2_left + align2_right
    full_score = score_left[split] + score_right[split]

    return full_align1, full_align2, int(full_score)