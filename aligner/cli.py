import argparse
import time
import os
import sys
try:
    import psutil
except ImportError:
    psutil = None
from aligner.algorithms import needleman_wunsch, smith_waterman
from aligner.io_utils import load_sequences, format_alignment, format_msa
from aligner.msa import multiple_sequence_alignment
from aligner.scoring import load_scoring_matrix
import logging


def main():
    parser = argparse.ArgumentParser(
        description="Инструмент для выравнивания геномных последовательностей. Поддерживает pairwise (global/local) и multiple alignments (MSA). Используйте --mode для выбора, --matrix для scoring (e.g., BLOSUM62). Для MSA - один FASTA с multiple seq."
    )
    parser.add_argument("--input1", required=True, help="Путь к FASTA-файлу (для MSA — один файл с multiple seq).")
    parser.add_argument("--input2", help="Путь к второму FASTA (только для pairwise).")
    parser.add_argument("--mode", choices=["global", "local", "msa"], default="global",
                        help="Режим: global (NW), local (SW), msa (progressive).")
    parser.add_argument("--output", default="alignment.txt", help="Файл для вывода (txt with alignment and stats).")
    parser.add_argument("--match", type=int, default=1, help="Score за совпадение (default 1).")
    parser.add_argument("--mismatch", type=int, default=-1, help="Score за несовпадение (default -1).")
    parser.add_argument("--gap", type=int, default=-2, help="Штраф за gap (default -2).")
    parser.add_argument("--matrix", default=None, help="Scoring matrix (e.g., BLOSUM62 for proteins).")
    parser.add_argument("--clustal", action="store_true", help="Вывод MSA в Clustal format (для --mode msa).")
    parser.add_argument("--verbose", action="store_true", help="Включить детальный logging (debug level).")
    parser.add_argument("--subsample", type=int, default=0,
                        help="Subsample first N bases (for large genomes, default 0 - full).")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    scoring_matrix = load_scoring_matrix(args.matrix) if args.matrix else None

    sequences = load_sequences(args.input1)

    if args.subsample > 0:
        sequences = [seq[:args.subsample] for seq in sequences]
        logging.info(f"Subsampled to first {args.subsample} bases.")

    start_time = time.time()
    memory_usage = 0
    if psutil:
        process = psutil.Process(os.getpid())
        start_mem = process.memory_info().rss / 1024 ** 2

    if args.mode == "msa":
        if len(sequences) < 2:
            raise ValueError("Для MSA нужен файл с несколькими последовательностями.")
        aligned = multiple_sequence_alignment(sequences, args.match, args.mismatch, args.gap, scoring_matrix)
        end_time = time.time()
        if psutil:
            memory_usage = process.memory_info().rss / 1024 ** 2 - start_mem

        result = format_msa(aligned, clustal=args.clustal)
        result += f"\nTime: {end_time - start_time:.2f} sec\nMemory: {memory_usage:.2f} MB"

        with open(args.output, "w") as f:
            f.write(result)
        print("MSA сохранено в", args.output)
        return

    # Pairwise
    if not args.input2:
        raise ValueError("Для pairwise нужен --input2.")
    seq1 = sequences[0]
    seq2 = load_sequences(args.input2)[0]

    if args.subsample > 0:
        seq1 = seq1[:args.subsample]
        seq2 = seq2[:args.subsample]
        logging.info(f"Subsampled pairwise to first {args.subsample} bases.")

    if args.mode == "global":
        align1, align2, score = needleman_wunsch(seq1, seq2, args.match, args.mismatch, args.gap, scoring_matrix)
    else:
        align1, align2, score = smith_waterman(seq1, seq2, args.match, args.mismatch, args.gap, scoring_matrix)

    end_time = time.time()
    if psutil:
        memory_usage = process.memory_info().rss / 1024 ** 2 - start_mem

    result = f"Alignment Score: {score}\n"
    result += format_alignment(align1, align2)
    result += f"\nTime: {end_time - start_time:.2f} sec\nMemory: {memory_usage:.2f} MB"

    with open(args.output, "w") as f:
        f.write(result)

    print("выравнивание сохранено в", args.output)


if __name__ == "__main__":
    main()