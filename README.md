# Genomic Sequence Aligner ENG

**Genomic Sequence Aligner** is a tool for aligning genomic sequences. The program is designed for analyzing and comparing DNA and protein sequences, supporting pairwise alignment and multiple sequence alignment (MSA).



**Genomic Sequence Aligner** is a tool for aligning genomic sequences with support for various methods and configuration parameters. Key features include:

- **Pairwise alignment**:
  - *Global alignment* (Needleman-Wunsch): aligns sequences over their entire length.
  - *Local alignment* (Smith-Waterman): finds the most similar regions of sequences.
- **Multiple sequence alignment (MSA)**: uses a progressive method to align multiple sequences.
- **Support for scoring matrices**: for example, BLOSUM62 for protein sequences.
- **Subsampling**: the ability to limit analysis to the first N bases for faster processing of large genomes.
- **Flexible output**: support for Clustal format for MSA and detailed logging for debugging.

## Program Structure

The program consists of modules, each performing a specific task:

- `algorithms.py`: Implementation of the Needleman-Wunsch and Smith-Waterman algorithms for pairwise alignment.
- `msa.py`: Implementation of the method for multiple sequence alignment (MSA).
- `scoring.py`: Loading of scoring matrices.
- `io_utils.py`: Functions for loading and saving sequences in FASTA format.
- `cli.py`: Command-line interface for running the program with parameters.

## Installation and Dependencies

### Requirements

- Python 3.6+
- Installed dependencies:
  - `numpy`
  - `biopython`
  - `requests`
  - `gzip`
  - `logging`
  - `time`
  - `os`
  - `sys`
  - `argparse`
  - `psutil`

## Usage

The program is run via the command line using the `cli.py` module. Basic syntax:

```bash
python -m aligner.cli [options]
```

### Required Parameters

- `--input1`: Path to the first FASTA file (for MSA, a file with multiple sequences).
- `--input2`: Path to the second FASTA file (only for pairwise alignment).
- `--mode`: Alignment mode (`global`, `local`, `msa`).
- `--output`: Path to the file to save the result.

### Optional Parameters

- `--match`: Score for a match (default 1).
- `--mismatch`: Score for a mismatch (default -1).
- `--gap`: Penalty for a gap (default -2).
- `--matrix`: Scoring matrix (e.g., BLOSUM62 for proteins).
- `--clustal`: Output MSA in Clustal format (only for `--mode msa`).
- `--verbose`: Enable detailed logging (debug level).
- `--subsample`: Use only the first N bases (default 0 — full analysis).

## Differences Between Alignment Modes

- **Global alignment (`global`)**:
  - Uses the Needleman-Wunsch algorithm.
  - Aligns sequences over their entire length.
  - Suitable for sequences of approximately the same length with a high degree of similarity.

- **Local alignment (`local`)**:
  - Uses the Smith-Waterman algorithm.
  - Finds and aligns the most similar regions of sequences.
  - Ideal for finding subsequences or local homologies.

- **Multiple sequence alignment (`msa`)**:
  - Uses a progressive method.
  - Aligns multiple sequences simultaneously.
  - Suitable for comparing multiple genomes or proteins, for example, for phylogenetic analysis.

## Test Files

The `data` directory contains a set of test FASTA files that can be used to verify the program's operation:

- `hemoglobin.fasta`: Human hemoglobin sequence (protein).
- `hemoglobin_alpha.fasta`: Mouse hemoglobin sequence (protein).
- `test_multi.fasta`: File with three sequences for MSA (human, mouse, chimpanzee).
- `Salmonella.fasta`: Complete Salmonella genome (DNA).
- `Escherichia_coli_str_K-12_substr_MG1655.fasta`: Complete E. coli genome (DNA).

## Example Commands

Below are examples of commands for various usage scenarios. File paths are relative to the `data` directory. Replace the paths with those relevant to your system.

1. **Global pairwise alignment**  
   Aligning human and mouse hemoglobin:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode global --output alignment_global.txt
   ```

2. **Local pairwise alignment**  
   Finding similar regions in human and mouse hemoglobin:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode local --output alignment_local.txt
   ```

3. **Multiple sequence alignment (MSA)**  
   Aligning three sequences from `test_multi.fasta`:

   ```bash
   python -m aligner.cli --input1 "data/test_multi.fasta" --mode msa --output msa.txt
   ```

4. **Pairwise alignment with BLOSUM62 matrix**  
   Global alignment of proteins using BLOSUM62:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode global --matrix BLOSUM62 --output alignment_blosum.txt
   ```

5. **Pairwise alignment with subsampling**  
   Global alignment of the first 1000 bases of Salmonella and E. coli genomes:

   ```bash
   python -m aligner.cli --input1 "data/Salmonella.fasta" --input2 "data/Escherichia_coli_str_K-12_substr_MG1655.fasta" --mode global --subsample 1000 --output alignment_subsampled.txt
   ```

6. **Multiple sequence alignment in Clustal format**  
   MSA with output in Clustal format:

   ```bash
   python -m aligner.cli --input1 "data/test_multi.fasta" --mode msa --clustal




## Genomic Sequence Aligner РУ

**Genomic Sequence Aligner** — инструмент для выравнивания геномных последовательностей. Программа предназначена для анализа и сравнения последовательностей ДНК и белков, поддерживая парное выравнивание (*pairwise alignment*) и множественное выравнивание (*multiple sequence alignment, MSA*).

# Общее описание

**Genomic Sequence Aligner** — это инструмент для выравнивания геномных последовательностей с поддержкой различных методов и параметров настройки. Основные возможности:

- **Парное выравнивание**:
  - *Глобальное выравнивание* (Needleman-Wunsch): выравнивает последовательности по всей длине.
  - *Локальное выравнивание* (Smith-Waterman): находит наиболее похожие участки последовательностей.
- **Множественное выравнивание (MSA)**: использует прогрессивный метод для выравнивания нескольких последовательностей.
- **Поддержка scoring matrices**: например, BLOSUM62 для белковых последовательностей.
- **Subsampling**: возможность ограничить анализ первыми N базами для ускорения работы с большими геномами.
- **Гибкость вывода**: поддержка формата Clustal для MSA и детализированный лог для отладки.

## Структура программы

Программа состоит из модулей, каждый из которых выполняет свою задачу:

- `algorithms.py`: Реализация алгоритмов Needleman-Wunsch и Smith-Waterman для парного выравнивания.
- `msa.py`: Реализация метода для множественного выравнивания (MSA).
- `scoring.py`: Загрузка *scoring matrices*.
- `io_utils.py`: Функции для загрузки и сохранения последовательностей в формате FASTA.
- `cli.py`: Командная строка для запуска программы с параметрами.

## Установка и зависимости

### Требования

- Python 3.6+
- Установленные зависимости:
  - `numpy`
  - `biopython`
  - `requests`
  - `gzip`
  - `logging`
  - `time`
  - `os`
  - `sys`
  - `argparse`
  - `psutil`

## Использование

Программа запускается через командную строку с использованием модуля `cli.py`. Основной синтаксис:

```bash
python -m aligner.cli [options]
```

### Обязательные параметры

- `--input1`: Путь к первому FASTA-файлу (для MSA — файл с несколькими последовательностями).
- `--input2`: Путь ко второму FASTA-файлу (только для парного выравнивания).
- `--mode`: Режим выравнивания (`global`, `local`, `msa`).
- `--output`: Путь к файлу для сохранения результата.

### Дополнительные параметры

- `--match`: Score за совпадение (по умолчанию 1).
- `--mismatch`: Score за несовпадение (по умолчанию -1).
- `--gap`: Штраф за gap (по умолчанию -2).
- `--matrix`: Scoring matrix (например, BLOSUM62 для белков).
- `--clustal`: Вывод MSA в формате Clustal (только для `--mode msa`).
- `--verbose`: Включить детализированный лог (уровень debug).
- `--subsample`: Использовать только первые N баз (по умолчанию 0 — полный анализ).

## Различия между режимами выравнивания

- **Глобальное выравнивание (`global`)**:
  - Использует алгоритм Needleman-Wunsch.
  - Выравнивает последовательности по всей длине.
  - Подходит для последовательностей примерно одинаковой длины с высокой степенью сходства.

- **Локальное выравнивание (`local`)**:
  - Использует алгоритм Smith-Waterman.
  - Находит и выравнивает наиболее похожие участки последовательностей.
  - Идеально для поиска подпоследовательностей или локальных гомологий.

- **Множественное выравнивание (`msa`)**:
  - Использует прогрессивный метод.
  - Выравнивает несколько последовательностей одновременно.
  - Подходит для сравнения нескольких геномов или белков, например, для филогенетического анализа.

## Тестовые файлы

В директории `data` представлен набор тестовых FASTA-файлов, которые можно использовать для проверки работы программы:

- `hemoglobin.fasta`: Последовательность гемоглобина человека (белок).
- `hemoglobin_alpha.fasta`: Последовательность гемоглобина мыши (белок).
- `test_multi.fasta`: Файл с тремя последовательностями для MSA (человек, мышь, шимпанзе).
- `Salmonella.fasta`: Полный геном Salmonella (ДНК).
- `Escherichia_coli_str_K-12_substr_MG1655.fasta`: Полный геном E. coli (ДНК).

## Примеры команд

Ниже приведены примеры команд для различных сценариев использования. Пути к файлам указаны относительно директории `data`. Замените пути на актуальные для вашей системы.

1. **Глобальное парное выравнивание**  
   Выравнивание гемоглобина человека и мыши:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode global --output alignment_global.txt
   ```

2. **Локальное парное выравнивание**  
   Поиск похожих участков гемоглобина человека и мыши:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode local --output alignment_local.txt
   ```

3. **Множественное выравнивание (MSA)**  
   Выравнивание трех последовательностей из `test_multi.fasta`:

   ```bash
   python -m aligner.cli --input1 "data/test_multi.fasta" --mode msa --output msa.txt
   ```

4. **Парное выравнивание с матрицей BLOSUM62**  
   Глобальное выравнивание белков с использованием BLOSUM62:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode global --matrix BLOSUM62 --output alignment_blosum.txt
   ```

5. **Парное выравнивание с subsampling**  
   Глобальное выравнивание первых 1000 баз геномов Salmonella и E. coli:

   ```bash
   python -m aligner.cli --input1 "data/Salmonella.fasta" --input2 "data/Escherichia_coli_str_K-12_substr_MG1655.fasta" --mode global --subsample 1000 --output alignment_subsampled.txt
   ```

6. **Множественное выравнивание в формате Clustal**  
   MSA с выводом в формате Clustal:

   ```bash
   python -m aligner.cli --input1 "data/test_multi.fasta" --mode msa --clustal --output msa_clustal.txt
   ```

7. **Парное выравнивание с детализированным логом**  
   Глобальное выравнивание с включением отладочной информации:

   ```bash
   python -m aligner.cli --input1 "data/hemoglobin.fasta" --input2 "data/hemoglobin_alpha.fasta" --mode global --verbose --output alignment_verbose.txt
   ```

## Как работает CLI

Командная строка (`cli.py`) — основной интерфейс программы. Вот как она обрабатывает запросы:

1. **Парсинг аргументов**: Использует `argparse` для обработки параметров, таких как `--input1`, `--mode` и т.д.
2. **Загрузка последовательностей**: Через `io_utils.py` загружает данные из FASTA-файлов.
3. **Обработка subsampling**: Если указан `--subsample`, обрезает последовательности до заданной длины.
4. **Выбор алгоритма**:
   - Для `--mode global`: Вызывает `needleman_wunsch`.
   - Для `--mode local`: Вызывает `smith_waterman`.
   - Для `--mode msa`: Вызывает `multiple_sequence_alignment`.
5. **Scoring**: Применяет *scoring matrix* (если указана через `--matrix`) или использует параметры `--match`, `--mismatch`, `--gap`.
6. **Вывод результата**: Форматирует выравнивание (`format_alignment` для парного или `format_msa` для MSA) и сохраняет в файл с указанием времени выполнения и потребления памяти (если установлен `psutil`).

## Дополнительные возможности

- **Scoring matrices**: Поддержка матриц, таких как BLOSUM62, позволяет учитывать биологические свойства аминокислот.
- **Subsampling**: Ускоряет анализ больших геномов, ограничивая выравнивание первыми N базами.
- **Формат Clustal**: Удобный для визуализации и совместимый с другими биоинформатическими инструментами.
- **Детальный лог**: Помогает в отладке и анализе работы программы.

## В ближайшем будущем добавлю UI
