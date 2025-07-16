import numpy as np
from typing import List, Tuple, Optional, Dict
from aligner.algorithms import needleman_wunsch, hirschberg_needleman_wunsch
import logging
import random
from multiprocessing import Pool, cpu_count


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class MSAError(Exception):

    pass


def pairwise_distance(args):
    # Функция для parallel вычисления расстояний
    i, j, sequences, match, mismatch, gap, gap_open, gap_extend, scoring_matrix = args
    seq_i, seq_j = sequences[i], sequences[j]
    if seq_i == seq_j:
        score = len(seq_i) * match  # For identical, max score
    elif max(len(seq_i), len(seq_j)) > 5000:
        _, _, score = hirschberg_needleman_wunsch(seq_i, seq_j, match, mismatch, gap, scoring_matrix)
    else:
        _, _, score = needleman_wunsch(seq_i, seq_j, match, mismatch, gap, gap_open, gap_extend, scoring_matrix)
    max_len = max(len(seq_i), len(seq_j))
    normalized = -score / max_len if max_len > 0 else 0
    return i, j, normalized


def compute_distance_matrix(
        sequences: List[str],
        match: int = 1,
        mismatch: int = -1,
        gap: int = -2,
        gap_open: Optional[int] = None,
        gap_extend: Optional[int] = None,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None,
        threads: int = cpu_count()
) -> np.ndarray:
    # дистанционная матрица с параллелизацией
    n = len(sequences)
    if n > 100:
        raise MSAError("Слишком много последовательностей для MSA (max 100)")
    dist = np.zeros((n, n))
    tasks = [(i, j, sequences, match, mismatch, gap, gap_open, gap_extend, scoring_matrix) for i in range(n) for j in range(i + 1, n)]

    with Pool(threads) as pool:
        results = pool.map(pairwise_distance, tasks)

    for i, j, normalized in results:
        dist[i][j] = dist[j][i] = normalized
    return dist


def build_guide_tree(dist: np.ndarray) -> List[Tuple[int, int]]:
    # guide tree
    n = dist.shape[0]
    clusters = [[i] for i in range(n)]
    tree = []
    current_dist = dist.copy()

    while len(clusters) > 1:
        nonzero = np.nonzero(current_dist)
        if len(nonzero[0]) == 0 or np.all(current_dist[nonzero] == 0):
            for k in range(1, len(clusters)):
                tree.append((0, k))
            break
        min_val = np.min(current_dist[nonzero])
        i, j = np.where(current_dist == min_val)
        i, j = i[0], j[0]
        if i > j:
            i, j = j, i
        tree.append((i, j))
        new_cluster = clusters[i] + clusters[j]
        clusters = [c for k, c in enumerate(clusters) if k not in (i, j)] + [new_cluster]
        new_dist = np.full((len(clusters), len(clusters)), np.inf)
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                scores = [current_dist[p, q] for p in clusters[a] for q in clusters[b] if p != q and 0 <= p < current_dist.shape[0] and 0 <= q < current_dist.shape[0]]
                if scores:
                    avg = np.mean(scores)
                else:
                    avg = 0.0
                new_dist[a, b] = new_dist[b, a] = avg if not np.isnan(avg) else 0.0
        current_dist = new_dist

    return tree


def get_consensus_columnwise(group: List[str]) -> str:
    consensus = ''
    for chars in zip(*group):
        non_gap_chars = [c for c in chars if c != '-']
        if non_gap_chars:
            counts = {c: non_gap_chars.count(c) for c in set(non_gap_chars)}
            max_count = max(counts.values())
            candidates = [c for c, count in counts.items() if count == max_count]
            consensus += random.choice(candidates) if len(candidates) > 1 else candidates[0]
        else:
            consensus += '-'
    return consensus


def progressive_align(
        sequences: List[str],
        tree: List[Tuple[int, int]],
        match: int = 1,
        mismatch: int = -1,
        gap: int = -2,
        gap_open: Optional[int] = None,
        gap_extend: Optional[int] = None,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None
) -> List[str]:
    alignments = [[s] for s in sequences]

    for merge in tree:
        logging.debug(f"Merging clusters {merge}")
        i, j = min(merge), max(merge)  # Ensure i < j
        align_i = alignments[i]
        align_j = alignments[j]

        cons_i = get_consensus_columnwise(align_i)
        cons_j = get_consensus_columnwise(align_j)

        prof_align1, prof_align2, _ = needleman_wunsch(cons_i, cons_j, match, mismatch, gap, gap_open, gap_extend, scoring_matrix)

        def expand_alignment(original_group, profile_alignment):
            expanded = []
            for seq in original_group:
                new_seq = []
                pos = 0
                for char in profile_alignment:
                    if char == '-':
                        new_seq.append('-')
                    else:
                        new_seq.append(seq[pos])
                        pos += 1
                expanded.append(''.join(new_seq))
            return expanded

        new_align_i = expand_alignment(align_i, prof_align1)
        new_align_j = expand_alignment(align_j, prof_align2)
        new_align = new_align_i + new_align_j

        alignments.pop(j)
        alignments.pop(i)
        alignments.append(new_align)

    final_align = alignments[0]
    consensus = get_consensus_columnwise(final_align)
    logging.info(f"Consensus sequence: {consensus}")
    return final_align


def multiple_sequence_alignment(
        sequences: List[str],
        match: int = 1,
        mismatch: int = -1,
        gap: int = -2,
        gap_open: Optional[int] = None,
        gap_extend: Optional[int] = None,
        scoring_matrix: Optional[Dict[Tuple[str, str], int]] = None,
        threads: int = cpu_count()
) -> List[str]:

    if len(sequences) < 2:
        raise MSAError("Нужны хотя бы 2 последовательности для MSA")
    dist = compute_distance_matrix(sequences, match, mismatch, gap, gap_open, gap_extend, scoring_matrix, threads)
    tree = build_guide_tree(dist)
    return progressive_align(sequences, tree, match, mismatch, gap, gap_open, gap_extend, scoring_matrix)