from Bio import SeqIO
from typing import List, Tuple
import os
import requests
import gzip
import logging
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_sequences(file_path: str) -> List[str]:
    """
     последовательности из FASTA (handle gz too).
    file_path: Путь к файлу .
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"файл не найден: {file_path}")

    sequences = []
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as handle:
                sequences = [str(record.seq) for record in SeqIO.parse(handle, "fasta")]
        else:
            sequences = [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]
    except Exception as e:
        logging.error(f"ошибка парсинга FASTA: {e}. проверьте, является ли файл валидным FASTA (не HTML).")
        raise

    if not sequences:
        logging.warning(
            f"нет последовательностей в файле {file_path}. возможно, скачали HTML вместо FASTA. Попробуйте скачать заново с headers.")
        raise ValueError(f"нет последовательностей в {file_path}. скачайте заново.")

    return sequences


def download_sequences(url: str, save_path: str, retries: int = 3) -> None:
    """
    Скачивает FASTA по URL и сохраняет (handle gz, with headers, retries).

    :param url: URL to FASTA (plain or gz).
    :param save_path: Path to save.
    :param retries: Number of tries on fail.
    """
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, stream=True)
            if response.status_code == 200:
                content = response.content
                if 'gzip' in response.headers.get('Content-Encoding', '') or url.endswith('.gz'):
                    with gzip.open(save_path, 'wb') as f:
                        f.write(content)
                    logging.info(f"Скачан gz файл: {save_path}")
                else:
                    with open(save_path, "wb") as f:
                        f.write(content)
                    logging.info(f"Скачан файл: {save_path}")
                # Check if FASTA
                with open(save_path, "r") as f:
                    head = f.read(10)
                    if not head.startswith('>'):
                        logging.warning(f"Файл {save_path} не начинается с '>', возможно HTML. Удаляю и retry.")
                        os.remove(save_path)
                        raise ValueError("Not FASTA")
                return
            else:
                logging.warning(f"Попытка {attempt + 1}: Ошибка {response.status_code}")
        except Exception as e:
            logging.warning(f"Попытка {attempt + 1}: Ошибка {e}")
        time.sleep(2)  # Delay
    raise ValueError(f"Не удалось скачать после {retries} попыток. Проверьте URL в браузере, возможно CAPTCHA.")


def format_alignment(align1: str, align2: str) -> str:
    match_line = ''.join('|' if a == b else ' ' for a, b in zip(align1, align2))
    return f"{align1}\n{match_line}\n{align2}"


def format_msa(alignments: List[str], clustal: bool = False) -> str:
    if clustal:
        result = "CLUSTAL format\n"
        for i in range(0, len(alignments[0]), 60):
            for j, align in enumerate(alignments):
                result += f"Seq{j + 1} {align[i:i + 60]}\n"
            result += "\n"
    else:
        result = ""
        for i, align in enumerate(alignments):
            result += f"Seq{i + 1}: {align}\n"
    return result