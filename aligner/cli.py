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


TRANSLATIONS = {
    'en': {
        'description': "Tool for aligning genomic sequences. Supports pairwise (global/local) and multiple alignments (MSA). Use --mode to select, --matrix for scoring (e.g., BLOSUM62). For MSA - one FASTA with multiple seq.",
        'input1': "Path to FASTA file (for MSA — one file with multiple seq).",
        'input2': "Path to second FASTA (only for pairwise).",
        'mode': "Mode: global (NW), local (SW), msa (progressive).",
        'output': "Output file (txt with alignment and stats).",
        'match': "Match score (default 1).",
        'mismatch': "Mismatch score (default -1).",
        'gap': "Gap penalty (default -2).",
        'gap_open': "Gap open penalty (for affine, default None - disable affine).",
        'gap_extend': "Gap extend penalty (for affine, default None).",
        'matrix': "Scoring matrix (e.g., BLOSUM62 for proteins).",
        'clustal': "Output MSA in Clustal format (for --mode msa).",
        'verbose': "Enable detailed logging (debug level).",
        'subsample': "Subsample first N bases (for large genomes, default 0 - full).",
        'threads': "Number of threads for MSA (default cpu_count).",
        'lang': "Language for texts (en/ru).",
        'error_msa': "For MSA need file with multiple sequences.",
        'error_pairwise': "For pairwise need --input2.",
        'saved': "Alignment saved in",
        'msa_saved': "MSA saved in",
        'time': "Time",
        'memory': "Memory",
        'sec': "sec",
        'mb': "MB",
        'subsampled': "Subsampled to first {0} bases."
    },
    'ru': {
        'description': "Инструмент для выравнивания геномных последовательностей. Поддерживает pairwise (global/local) и multiple alignments (MSA). Используйте --mode для выбора, --matrix для scoring (e.g., BLOSUM62). Для MSA - один FASTA с multiple seq.",
        'input1': "Путь к FASTA-файлу (для MSA — один файл с multiple seq).",
        'input2': "Путь к второму FASTA (только для pairwise).",
        'mode': "Режим: global (NW), local (SW), msa (progressive).",
        'output': "Файл для вывода (txt with alignment and stats).",
        'match': "Score за совпадение (default 1).",
        'mismatch': "Score за несовпадение (default -1).",
        'gap': "Штраф за gap (default -2).",
        'gap_open': "Штраф за открытие gap (для affine, default None - отключить affine).",
        'gap_extend': "Штраф за расширение gap (для affine, default None).",
        'matrix': "Scoring matrix (e.g., BLOSUM62 for proteins).",
        'clustal': "Вывод MSA в Clustal format (для --mode msa).",
        'verbose': "Включить детальный logging (debug level).",
        'subsample': "Subsample first N bases (for large genomes, default 0 - full).",
        'threads': "Количество потоков для MSA (default cpu_count).",
        'lang': "Язык для текстов (en/ru).",
        'error_msa': "Для MSA нужен файл с несколькими последовательностями.",
        'error_pairwise': "Для pairwise нужен --input2.",
        'saved': "выравнивание сохранено в",
        'msa_saved': "MSA сохранено в",
        'time': "Время",
        'memory': "Память",
        'sec': "сек",
        'mb': "МБ",
        'subsampled': "Subsampled to first {0} bases."
    }
}


def main():
    # Основной CLI с subcommands и lang
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--lang', default='en', choices=['en', 'ru'], help=TRANSLATIONS['en']['lang'])
    args, remaining = parser.parse_known_args()
    lang = args.lang
    tr = TRANSLATIONS[lang]

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--output', default="alignment.txt", help=tr['output'])
    parent_parser.add_argument('--match', type=int, default=1, help=tr['match'])
    parent_parser.add_argument('--mismatch', type=int, default=-1, help=tr['mismatch'])
    parent_parser.add_argument('--gap', type=int, default=-2, help=tr['gap'])
    parent_parser.add_argument('--gap_open', type=int, default=None, help=tr['gap_open'])
    parent_parser.add_argument('--gap_extend', type=int, default=None, help=tr['gap_extend'])
    parent_parser.add_argument('--matrix', default=None, help=tr['matrix'])
    parent_parser.add_argument('--verbose', action="store_true", help=tr['verbose'])
    parent_parser.add_argument('--subsample', type=int, default=0, help=tr['subsample'])
    parent_parser.add_argument('--lang', default='en', choices=['en', 'ru'], help=tr['lang'])

    parser = argparse.ArgumentParser(description=tr['description'])
    subparsers = parser.add_subparsers(dest='mode', required=False)  # Changed to False to allow no mode

    global_parser = subparsers.add_parser('global', parents=[parent_parser])
    global_parser.add_argument("--input1", required=True, help=tr['input1'])
    global_parser.add_argument("--input2", required=True, help=tr['input2'])

    local_parser = subparsers.add_parser('local', parents=[parent_parser])
    local_parser.add_argument("--input1", required=True, help=tr['input1'])
    local_parser.add_argument("--input2", required=True, help=tr['input2'])

    msa_parser = subparsers.add_parser('msa', parents=[parent_parser])
    msa_parser.add_argument("--input1", required=True, help=tr['input1'])
    msa_parser.add_argument("--clustal", action="store_true", help=tr['clustal'])
    msa_parser.add_argument("--threads", type=int, default=os.cpu_count(), help=tr['threads'])

    # Parse remaining args with the full parser
    args = parser.parse_args(remaining)
    # Set lang from early parse (overrides default if specified)
    args.lang = lang

    # If no mode provided, show help and exit
    if not args.mode:
        parser.print_help()
        sys.exit(0)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    scoring_matrix = load_scoring_matrix(args.matrix) if args.matrix else None

    sequences = load_sequences(args.input1)

    if args.subsample > 0:
        sequences = [seq[:args.subsample] for seq in sequences]
        logging.info(tr['subsampled'].format(args.subsample))

    start_time = time.time()
    memory_usage = 0
    if psutil:
        process = psutil.Process(os.getpid())
        start_mem = process.memory_info().rss / 1024 ** 2

    if args.mode == "msa":
        if len(sequences) < 2:
            raise ValueError(tr['error_msa'])
        aligned = multiple_sequence_alignment(sequences, args.match, args.mismatch, args.gap, args.gap_open, args.gap_extend, scoring_matrix, args.threads)
        end_time = time.time()
        if psutil:
            memory_usage = process.memory_info().rss / 1024 ** 2 - start_mem

        result = format_msa(aligned, clustal=args.clustal)
        result += f"\n{tr['time']}: {end_time - start_time:.2f} {tr['sec']}\n{tr['memory']}: {memory_usage:.2f} {tr['mb']}"

        with open(args.output, "w") as f:
            f.write(result)
        print(tr['msa_saved'], args.output)
        return

    # Pairwise
    if not hasattr(args, 'input2'):
        raise ValueError(tr['error_pairwise'])
    seq1 = sequences[0]
    seq2 = load_sequences(args.input2)[0]

    if args.subsample > 0:
        seq1 = seq1[:args.subsample]
        seq2 = seq2[:args.subsample]
        logging.info(tr['subsampled'].format(args.subsample))

    if args.mode == "global":
        align1, align2, score = needleman_wunsch(seq1, seq2, args.match, args.mismatch, args.gap, args.gap_open, args.gap_extend, scoring_matrix)
    else:
        align1, align2, score = smith_waterman(seq1, seq2, args.match, args.mismatch, args.gap, scoring_matrix)

    end_time = time.time()
    if psutil:
        memory_usage = process.memory_info().rss / 1024 ** 2 - start_mem

    result = f"Alignment Score: {score}\n"
    result += format_alignment(align1, align2)
    result += f"\n{tr['time']}: {end_time - start_time:.2f} {tr['sec']}\n{tr['memory']}: {memory_usage:.2f} {tr['mb']}"

    with open(args.output, "w") as f:
        f.write(result)

    print(tr['saved'], args.output)


if __name__ == "__main__":
    main()