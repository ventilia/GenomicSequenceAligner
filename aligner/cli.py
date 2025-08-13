import os
import sys
import time
import logging
from typing import List, Optional, Dict

try:
    import psutil
except ImportError:
    psutil = None
import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress
from rich.panel import Panel
from rich.text import Text
import inquirer
import yaml
from aligner.algorithms import needleman_wunsch, smith_waterman
from aligner.io_utils import load_sequences, format_alignment, format_msa
from aligner.msa import multiple_sequence_alignment
from aligner.scoring import load_scoring_matrix

console = Console()

# ascii-арт с цветами
ASCII_ART = Text("""
                                      ++++## 
                                     =*%*+#%%
                       ==++****+**+++++#**#%%
                     ==+#%%%@%*#####%%%%+*%%%
                     ==*#%@@%###*****#%%+#%% 
                     =+*#@@% #****  %#*+*%%  
                     =+*%%@%%% ##***  ++#%%  
                     =+*#%%%%%%%  #***+#%%   
                     ++*%%@%% %%%%  +*#%%    
                     ++*%% %%%  %%%*#%%%     
                   ++#+*#%  %%%#*##%%%       
       =++*##%%+==++**+##%%##%%%%%%%         
     =+*#%%#++*%%%@@%%+*%%%%%%%%%            
   ==+#%%@## #**+++**++#%%                   
   =+###  #***    ##*+*%%%                   
   =*#%@%#  #*+**   ++#%%                    
   =*## %%%%   ****++#%%                     
   +*#%%  %%%%   #*+*%%%                     
   +*%%@%#  %%%% +*#%%                       
   +*%% %%%% +#%##%%%                        
##*+*%%***%@@%%%%%%                          
+**+#%%@@%%%%%%%                             
 #*+#%%                                      
*  +#%% 
  ______                                          __                     
 /      \                                        |  \                    
|  $$$$$$\ ______  _______   ______  ______ ____  \$$ _______            
| $$ __\$$/      \|       \ /      \|      \    \|  \/       \           
| $$|    |  $$$$$$| $$$$$$$|  $$$$$$| $$$$$$\$$$$| $|  $$$$$$$           
| $$ \$$$| $$    $| $$  | $| $$  | $| $$ | $$ | $| $| $$                 
| $$__| $| $$$$$$$| $$  | $| $$__/ $| $$ | $$ | $| $| $$_____            
 \$$    $$\$$     | $$  | $$\$$    $| $$ | $$ | $| $$\$$     \           
  \$$$$$$  \$$$$$$$\$$   \$$ \$$$$$$ \$$  \$$  \$$\$$ \$$$$$$$           
  ______  __ __                                                    
 /      \|  |  \                                                     
|  $$$$$$| $$\$$ ______  _______   ______   ______                       
| $$__| $| $|  \/      \|       \ /      \ /      \                      
| $$    $| $| $|  $$$$$$| $$$$$$$|  $$$$$$|  $$$$$$\                     
| $$$$$$$| $| $| $$  | $| $$  | $| $$    $| $$   \$$                     
| $$  | $| $| $| $$__| $| $$  | $| $$$$$$$| $$                           
| $$  | $| $| $$\$$    $| $$  | $$\$$     | $$                           
 \$$   \$$\$$\$$_\$$$$$$$\$$   \$$ \$$$$$$$\$$                           
               |  \__| $$                                                
                \$$    $$                                                
                 \$$$$$$                                                 
""", style="bold green")

TRANSLATIONS = {
    'en': {
        'welcome': "Welcome to Aligner CLI!",
        'choose_lang': "Choose language (en/ru):",
        'choose_mode': "Choose mode:",
        'mode_global': "Global (Needleman-Wunsch): Full alignment of entire sequences. Best for aligning sequences of similar length and high similarity (e.g., orthologous genes).",
        'mode_local': "Local (Smith-Waterman): Finds the best local alignment (subsequence match). Ideal for finding conserved domains or motifs within larger, possibly unrelated sequences.",
        'mode_msa': "MSA (Multiple Sequence Alignment): Aligns 2 or more sequences simultaneously. Essential for phylogenetic tree construction and identifying conserved regions across a family of sequences.",
        'select_dir': "Select directory containing your FASTA files (default: current directory):",
        'select_file1': "Select the first FASTA file (input1):",
        'select_file2': "Select the second FASTA file (input2) for pairwise alignment:",
        'select_matrix': "Scoring Matrix (optional):\n  - For proteins: BLOSUM62, PAM250, etc.\n  - For DNA/RNA: Leave blank (None) to use simple match/mismatch scores.\n  Enter matrix name or 'None':",
        'match_score': "Match Score (integer, default 1):\n  - Bonus for identical characters.\n  - DNA/RNA: Often 1 or 2.\n  - Proteins: Often 2-5. Ignored if a matrix is used.",
        'mismatch_score': "Mismatch Score (integer, default -1):\n  - Penalty for different characters.\n  - DNA/RNA: Often -1 or -2.\n  - Proteins: Often -1 to -4. Ignored if a matrix is used.",
        'gap_penalty': "Gap Penalty (integer, default -2):\n  - Penalty for introducing a gap ('-').\n  - Linear gap cost: Total gap penalty = Gap Open + (Gap Length * Gap Extend).\n  - This is the 'Gap Extend' if 'Gap Open' is also set.",
        'gap_open': "Gap Open Penalty (integer, optional for affine gaps):\n  - Additional penalty for starting a new gap.\n  - Enables more biologically realistic affine gap penalties.\n  - Leave blank (None) for simple linear gaps.\n  - If set, 'Gap Penalty' becomes 'Gap Extend'.",
        'gap_extend': "Gap Extend Penalty (integer, optional for affine gaps):\n  - Penalty for each additional space in a gap.\n  - Must be set if 'Gap Open' is used.\n  - Leave blank (None) for simple linear gaps.",
        'subsample': "Subsample Length (integer, 0 for full sequence):\n  - Process only the first N characters of each sequence.\n  - Useful for quickly testing large files.\n  - Enter 0 to use the entire sequence:",
        'threads': "Number of Threads for MSA (integer, default is number of CPU cores):",
        'clustal': "Output MSA in Clustal Format? (y/n)\n  - Clustal is a standard format for viewing multiple sequence alignments.",
        'verbose': "Enable Verbose Logging? (y/n)\n  - Prints detailed steps and debug information to the console.",
        'preview_seq': "Preview First 100 Characters of Sequences? (y/n)\n  - Shows a snippet of the input sequences before alignment starts.",
        'tutorial': "Run Tutorial with Example Alignments? (y/n)",
        'batch_mode': "Run Batch Alignment on All FASTA Files in Directory? (y/n)\n  - Performs pairwise alignments between all unique pairs of files in the directory.",
        'config': "Load configuration from YAML file (path):",
        'processing': "Processing alignment...",
        'success': "Alignment completed successfully!",
        'error': "Error:",
        'error_msa': "MSA requires a FASTA file with at least 2 sequences.",
        'error_pairwise': "Pairwise alignment requires two FASTA files.",
        'error_no_files': "No FASTA files found in the selected directory. Please check the path or provide files.",
        'error_invalid_file': "Invalid file selected. Please choose a file with a .fasta, .fa, or .gz extension.",
        'error_negative': "Penalty values (gap, mismatch, gap_open, gap_extend) must be negative integers.",
        'time': "Execution Time",
        'memory': "Memory Usage",
        'sec': "seconds",
        'mb': "MB",
        'identity': "Sequence Identity",
        'gaps': "Number of Gaps",
        'config_load': "Load config from YAML file? (enter path or 'n' for none):",
        'config_save': "Save current parameters to YAML file? (y/n)",
        'config_saved': "Configuration saved to {path}",
        'output_file': "Output file path:",
        'subsampled': "Subsampled to first {0} bases",
        'affine_gap_note': "Note: Affine gap penalties (Gap Open & Gap Extend) are only supported for Global alignment.",
        'alignment_header': "Alignment Results",
        'msa_header': "Multiple Sequence Alignment",
        'pairwise_header': "Pairwise Alignment",
        'score_label': "Alignment Score",
        'identity_label': "Identity Percentage",
        'gaps_label': "Total Gap Characters",
        'time_label': "Processing Time",
        'memory_label': "Peak Memory Usage",
        'preview_title': "Sequence Preview",
        'stats_title': "Alignment Statistics"
    },
    'ru': {
        'welcome': "Добро пожаловать в Aligner CLI!",
        'choose_lang': "Выберите язык (en/ru):",
        'choose_mode': "Выберите режим выравнивания:",
        'mode_global': "Global (Needleman-Wunsch): Полное выравнивание всей последовательности. Лучше всего подходит для выравнивания последовательностей одинаковой длины с высоким сходством (например, ортологичные гены).",
        'mode_local': "Local (Smith-Waterman): Находит лучшее локальное выравнивание (совпадение подпоследовательности). Идеально для поиска консервативных доменов или мотивов внутри больших, возможно несвязанных последовательностей.",
        'mode_msa': "MSA (Множественное Выравнивание Последовательностей): Одновременное выравнивание 2 или более последовательностей. Необходимо для построения филогенетических деревьев и выявления консервативных регионов в семействе последовательностей.",
        'select_dir': "Выберите директорию с вашими FASTA-файлами (по умолчанию: текущая директория):",
        'select_file1': "Выберите первый FASTA-файл (input1):",
        'select_file2': "Выберите второй FASTA-файл (input2) для парного выравнивания:",
        'select_matrix': "Scoring Matrix (необязательно):\n  - Для белков: BLOSUM62, PAM250 и т.д.\n  - Для ДНК/РНК: Оставьте пустым (None) для использования простых счетов совпадения/несовпадения.\n  Введите название матрицы или 'None':",
        'match_score': "Счет за Совпадение (целое число, по умолчанию 1):\n  - Бонус за одинаковые символы.\n  - ДНК/РНК: Часто 1 или 2.\n  - Белки: Часто 2-5. Игнорируется, если используется матрица.",
        'mismatch_score': "Счет за Несовпадение (целое число, по умолчанию -1):\n  - Штраф за разные символы.\n  - ДНК/РНК: Часто -1 или -2.\n  - Белки: Часто -1 до -4. Игнорируется, если используется матрица.",
        'gap_penalty': "Штраф за Gap (целое число, по умолчанию -2):\n  - Штраф за вставку пробела ('-').\n  - Линейная стоимость gap: Общий штраф = Открытие Gap + (Длина Gap * Расширение Gap).\n  - Это 'Расширение Gap', если также установлено 'Открытие Gap'.",
        'gap_open': "Штраф за Открытие Gap (целое число, необязательно для аффинных gaps):\n  - Дополнительный штраф за начало нового gap.\n  - Позволяет использовать более биологически реалистичные аффинные штрафы.\n  - Оставьте пустым (None) для простых линейных gaps.\n  - Если установлено, 'Штраф за Gap' становится 'Расширение Gap'.",
        'gap_extend': "Штраф за Расширение Gap (целое число, необязательно для аффинных gaps):\n  - Штраф за каждый дополнительный символ в gap.\n  - Должно быть установлено, если используется 'Открытие Gap'.\n  - Оставьте пустым (None) для простых линейных gaps.",
        'subsample': "Длина Подвыборки (целое число, 0 для полной последовательности):\n  - Обрабатывать только первые N символов каждой последовательности.\n  - Полезно для быстрого тестирования больших файлов.\n  - Введите 0, чтобы использовать всю последовательность:",
        'threads': "Количество Потоков для MSA (целое число, по умолчанию количество ядер CPU):",
        'clustal': "Вывести MSA в формате Clustal? (y/n)\n  - Clustal - стандартный формат для просмотра множественных выравниваний.",
        'verbose': "Включить Подробный Лог? (y/n)\n  - Выводит подробные шаги и отладочную информацию в консоль.",
        'preview_seq': "Предпросмотр Первых 100 Символов Последовательностей? (y/n)\n  - Показывает фрагмент входных последовательностей перед началом выравнивания.",
        'tutorial': "Запустить Обучающий Режим с Примерами? (y/n)",
        'batch_mode': "Запустить Пакетное Выравнивание для Всех FASTA-файлов в Директории? (y/n)\n  - Выполняет парные выравнивания между всеми уникальными парами файлов в директории.",
        'config': "Загрузить конфигурацию из YAML файла (путь):",
        'processing': "Обработка выравнивания...",
        'success': "Выравнивание успешно завершено!",
        'error': "Ошибка:",
        'error_msa': "Для MSA требуется FASTA-файл с минимум 2 последовательностями.",
        'error_pairwise': "Для парного выравнивания нужны два FASTA-файла.",
        'error_no_files': "FASTA-файлы не найдены в выбранной директории. Проверьте путь или предоставьте файлы.",
        'error_invalid_file': "Выбран неверный файл. Выберите файл с расширением .fasta, .fa или .gz.",
        'error_negative': "Штрафы (gap, mismatch, gap_open, gap_extend) должны быть отрицательными целыми числами.",
        'time': "Время Выполнения",
        'memory': "Использование Памяти",
        'sec': "секунд",
        'mb': "МБ",
        'identity': "Идентичность Последовательностей",
        'gaps': "Количество Пробелов",
        'config_load': "Загрузить конфигурацию из YAML? (введите путь или 'n' для пропуска):",
        'config_save': "Сохранить текущие параметры в YAML? (y/n)",
        'config_saved': "Конфигурация сохранена в {path}",
        'output_file': "Путь к выходному файлу:",
        'subsampled': "Обрезано до первых {0} баз",
        'affine_gap_note': "Примечание: Аффинные штрафы (Открытие Gap & Расширение Gap) поддерживаются только для Global выравнивания.",
        'alignment_header': "Результаты Выравнивания",
        'msa_header': "Множественное Выравнивание Последовательностей",
        'pairwise_header': "Парное Выравнивание",
        'score_label': "Счет Выравнивания",
        'identity_label': "Процент Идентичности",
        'gaps_label': "Общее Количество Пробелов",
        'time_label': "Время Обработки",
        'memory_label': "Пиковое Использование Памяти",
        'preview_title': "Предпросмотр Последовательностей",
        'stats_title': "Статистика Выравнивания"
    }
}


def get_fasta_files(directory: str) -> List[str]:
    # возвращает список fasta/gz файлов в директории
    return [f for f in os.listdir(directory) if f.endswith(('.fasta', '.fa', '.gz'))]


def validate_file(file_path: str, tr: Dict) -> bool:
    # проверяет, существует ли файл и является ли он fasta/gz
    if not os.path.exists(file_path):
        console.print(f"{tr['error']} File {file_path} does not exist.", style="bold red")
        return False
    if not file_path.endswith(('.fasta', '.fa', '.gz')):
        console.print(f"{tr['error']} {tr['error_invalid_file']}", style="bold red")
        return False
    return True


def validate_params(params: Dict, tr: Dict) -> bool:
    # проверяет, что штрафы отрицательные
    for param in ['gap', 'mismatch', 'gap_open', 'gap_extend']:
        # Исправлено: корректная проверка None и типа
        if param in params and params[param] is not None:
            if not isinstance(params[param], int):
                console.print(f"{tr['error']} Invalid value for {param}. Must be an integer.", style="bold red")
                return False
            if params[param] > 0:
                console.print(f"{tr['error']} {tr['error_negative']}", style="bold red")
                return False
    # Проверка аффинных gaps
    if params['mode'] != 'global':
        if params.get('gap_open') is not None or params.get('gap_extend') is not None:
            console.print(f"{tr['error']} {tr['affine_gap_note']}", style="bold red")
            return False
    return True


def interactive_wizard() -> tuple[Dict, Dict]:
    # интерактивный wizard для параметров
    params = {}
    lang_choice = \
        inquirer.prompt([inquirer.List('lang', message=TRANSLATIONS['en']['choose_lang'], choices=['en', 'ru'])])[
            'lang']
    tr = TRANSLATIONS[lang_choice]  # Используем выбранный язык для всего последующего взаимодействия
    console.print(Panel(tr['welcome'], style="bold blue"))

    # режим с объяснениями
    mode_choices = [
        (tr['mode_global'], 'global'),
        (tr['mode_local'], 'local'),
        (tr['mode_msa'], 'msa')
    ]
    params['mode'] = inquirer.prompt([inquirer.List('mode', message=tr['choose_mode'], choices=mode_choices)])['mode']

    # batch mode?
    params['batch'] = inquirer.prompt([inquirer.Confirm('batch', message=tr['batch_mode'], default=False)])['batch']

    # директория и файлы
    while True:
        directory = inquirer.prompt(
            [inquirer.Path('dir', message=tr['select_dir'], path_type=inquirer.Path.DIRECTORY, default=os.getcwd())])[
            'dir']
        fasta_files = get_fasta_files(directory)
        if not fasta_files and not params['batch']:
            console.print(f"{tr['error']} {tr['error_no_files']}", style="bold red")
            continue
        break

    if not params['batch']:
        file1_choice = inquirer.prompt([inquirer.List('input1', message=tr['select_file1'], choices=fasta_files)])[
            'input1']
        params['input1'] = os.path.join(directory, file1_choice)
        if not validate_file(params['input1'], tr):
            sys.exit(1)
        if params['mode'] != 'msa':
            file2_choice = inquirer.prompt([inquirer.List('input2', message=tr['select_file2'], choices=fasta_files)])[
                'input2']
            params['input2'] = os.path.join(directory, file2_choice)
            if not validate_file(params['input2'], tr):
                sys.exit(1)
    else:
        # В batch режиме передаем директорию
        params['directory'] = directory
        # input_files будут загружены позже в run_batch_alignment

    # параметры выравнивания
    matrix_input = inquirer.prompt([inquirer.Text('matrix', message=tr['select_matrix'], default='None')])['matrix']
    params['matrix'] = matrix_input if matrix_input and matrix_input.lower() != 'none' else None

    params['match'] = int(inquirer.prompt([inquirer.Text('match', message=tr['match_score'], default='1')])['match'])
    params['mismatch'] = int(
        inquirer.prompt([inquirer.Text('mismatch', message=tr['mismatch_score'], default='-1')])['mismatch'])
    params['gap'] = int(inquirer.prompt([inquirer.Text('gap', message=tr['gap_penalty'], default='-2')])['gap'])

    gap_open_input = inquirer.prompt([inquirer.Text('gap_open', message=tr['gap_open'], default='None')])['gap_open']
    # Исправлено: корректное преобразование строки "None" в Python None
    params['gap_open'] = int(gap_open_input) if gap_open_input and gap_open_input.lower() != 'none' else None

    gap_extend_input = inquirer.prompt([inquirer.Text('gap_extend', message=tr['gap_extend'], default='None')])[
        'gap_extend']
    # Исправлено: корректное преобразование строки "None" в Python None
    params['gap_extend'] = int(gap_extend_input) if gap_extend_input and gap_extend_input.lower() != 'none' else None

    params['subsample'] = int(
        inquirer.prompt([inquirer.Text('subsample', message=tr['subsample'], default='0')])['subsample'])

    if params['mode'] == 'msa':
        params['threads'] = int(
            inquirer.prompt([inquirer.Text('threads', message=tr['threads'], default=str(os.cpu_count()))])['threads'])
        params['clustal'] = inquirer.prompt([inquirer.Confirm('clustal', message=tr['clustal'], default=False)])[
            'clustal']

    params['verbose'] = inquirer.prompt([inquirer.Confirm('verbose', message=tr['verbose'], default=False)])['verbose']
    params['preview'] = inquirer.prompt([inquirer.Confirm('preview', message=tr['preview_seq'], default=False)])[
        'preview']
    params['tutorial'] = inquirer.prompt([inquirer.Confirm('tutorial', message=tr['tutorial'], default=False)])[
        'tutorial']
    params['output'] = inquirer.prompt([inquirer.Text('output', message=tr['output_file'], default="alignment.txt")])[
        'output']

    # config load/save
    config_load = inquirer.prompt([inquirer.Text('config_load', message=tr['config_load'], default='n')])['config_load']
    if config_load != 'n':
        if os.path.exists(config_load):
            params = load_config(config_load)
            # Убедимся, что язык из конфига совместим, или используем выбранный
            # В данном случае язык уже выбран, поэтому просто обновляем tr
            tr = TRANSLATIONS[lang_choice]
        else:
            console.print(f"{tr['error']} Config file {config_load} not found.", style="bold red")
            sys.exit(1)

    if inquirer.prompt([inquirer.Confirm('config_save', message=tr['config_save'], default=False)])['config_save']:
        save_config(params, 'config.yaml')

    if not validate_params(params, tr):
        sys.exit(1)

    return params, tr


def run_tutorial(tr: Dict):
    # tutorial с примерами выравниваний
    console.print(Panel(
        "Example 1: Global alignment (Needleman-Wunsch)\n"
        "Sequences: 'AGC' and 'ACGC'\n"
        "Score: 1\n"
        "A-GC\n"
        " | |\n"
        "ACGC\n"
        "Example 2: Local alignment (Smith-Waterman)\n"
        "Sequences: 'ATGC' and 'TGCA'\n"
        "Score: 3\n"
        "TGC\n"
        "|||\n"
        "TGC\n"
        "Example 3: MSA (Multiple Sequence Alignment)\n"
        "Seq1: A-GC\n"
        "Seq2: ACGC\n"
        "Seq3: AGGC",
        title=tr['tutorial'], style="bold blue"
    ))


def load_config(path: str) -> Dict:
    # загружаем конфигурацию из yaml
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def save_config(params: Dict, path: str):
    # сохраняем параметры в yaml
    with open(path, 'w') as f:
        yaml.dump(params, f)
    console.print(f"{TRANSLATIONS[params.get('lang', 'en')]['config_saved'].format(path=path)}", style="green")


def print_alignment_table(align1: str, align2: str, tr: Dict):
    # таблица для pairwise выравнивания с цветами
    table = Table(title=tr['pairwise_header'])
    table.add_column("Seq1", style="cyan")
    table.add_column("Matches", style="magenta")
    table.add_column("Seq2", style="green")
    match_line = ''.join(['|' if a == b and a != '-' else ' ' for a, b in zip(align1, align2)])
    for i in range(0, len(align1), 60):
        table.add_row(align1[i:i + 60], match_line[i:i + 60], align2[i:i + 60])
    console.print(table)


def compute_stats(align1: str, align2: str) -> Dict:
    # статистика: % идентичности, gaps
    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    length = len(align1)
    gaps = align1.count('-') + align2.count('-')
    identity = (matches / length * 100) if length > 0 else 0
    return {'identity': identity, 'gaps': gaps}


def run_batch_alignment(directory: str, params: Dict, tr: Dict) -> str:
    # batch-режим: pairwise все-против-всех
    fasta_files = get_fasta_files(directory)
    if len(fasta_files) < 2:
        console.print(f"{tr['error']} {tr['error_pairwise']}", style="bold red")
        sys.exit(1)
    result = ""
    scoring_matrix = load_scoring_matrix(params['matrix']) if params['matrix'] else None
    with Progress() as progress:
        task = progress.add_task(tr['processing'], total=len(fasta_files) * (len(fasta_files) - 1) // 2)
        for i, file1 in enumerate(fasta_files):
            try:
                seqs1 = load_sequences(os.path.join(directory, file1))
                if not seqs1:
                    console.print(f"{tr['error']} No sequences found in {file1}. Skipping.", style="bold yellow")
                    continue
                seq1 = seqs1[0]  # Берем первую последовательность из файла
            except Exception as e:
                console.print(f"{tr['error']} Failed to load {file1}: {e}. Skipping.", style="bold red")
                continue

            if params['subsample'] > 0:
                seq1 = seq1[:params['subsample']]

            for j, file2 in enumerate(fasta_files[i + 1:], start=i + 1):
                try:
                    seqs2 = load_sequences(os.path.join(directory, file2))
                    if not seqs2:
                        console.print(f"{tr['error']} No sequences found in {file2}. Skipping.", style="bold yellow")
                        continue
                    seq2 = seqs2[0]  # Берем первую последовательность из файла
                except Exception as e:
                    console.print(f"{tr['error']} Failed to load {file2}: {e}. Skipping.", style="bold red")
                    continue

                if params['subsample'] > 0:
                    seq2 = seq2[:params['subsample']]

                console.print(f"\nProcessing: {file1} vs {file2}", style="bold blue")
                try:
                    if params['mode'] == 'global':
                        align1, align2, score = needleman_wunsch(
                            seq1, seq2, params['match'], params['mismatch'], params['gap'],
                            params.get('gap_open'), params.get('gap_extend'), scoring_matrix
                        )
                    else:  # local
                        align1, align2, score = smith_waterman(
                            seq1, seq2, params['match'], params['mismatch'], params['gap'], scoring_matrix
                        )
                    result += f"\nAlignment: {file1} vs {file2}\nScore: {score}\n"
                    print_alignment_table(align1, align2, tr)
                    stats = compute_stats(align1, align2)
                    result += f"{tr['identity']}: {stats['identity']:.2f}%\n{tr['gaps']}: {stats['gaps']}\n"
                except Exception as e:
                    console.print(f"{tr['error']} Alignment failed for {file1} vs {file2}: {e}", style="bold red")
                    result += f"\nAlignment: {file1} vs {file2}\nStatus: Failed - {e}\n"
                progress.update(task, advance=1)
    return result


@click.group(invoke_without_command=True)
@click.option('--config', type=str, help=TRANSLATIONS['en']['config'])
@click.pass_context
def cli(ctx, config):
    # дефолт: wizard с ascii-арт
    console.print(ASCII_ART)
    if config:
        if os.path.exists(config):
            params = load_config(config)
            # Исправлено: добавлены дефолтные значения для отсутствующих ключей
            params.setdefault('verbose', False)
            params.setdefault('preview', False)
            params.setdefault('tutorial', False)
            params.setdefault('batch', False)
            params.setdefault('clustal', False)
            params.setdefault('subsample', 0)
            params.setdefault('matrix', None)
            # Язык из конфига или по умолчанию
            lang = params.get('lang', 'en')
            tr = TRANSLATIONS[lang]
            if not validate_params(params, tr):
                sys.exit(1)
            if params.get('tutorial', False):
                run_tutorial(tr)
            run_alignment(params, tr)
        else:
            console.print(f"Error: Config file {config} not found.", style="bold red")
            sys.exit(1)
    elif ctx.invoked_subcommand is None:
        params, tr = interactive_wizard()
        if params['tutorial']:
            run_tutorial(tr)
        run_alignment(params, tr)
    else:
        pass


@cli.command(name='global')
@click.option('--input1', type=str, help=TRANSLATIONS['en']['select_file1'])
@click.option('--input2', type=str, help=TRANSLATIONS['en']['select_file2'])
@click.option('--directory', type=str, help="Directory for batch alignment")
@click.option('--output', default="alignment.txt", help="Output file")
@click.option('--match', type=int, default=1, help=TRANSLATIONS['en']['match_score'])
@click.option('--mismatch', type=int, default=-1, help=TRANSLATIONS['en']['mismatch_score'])
@click.option('--gap', type=int, default=-2, help=TRANSLATIONS['en']['gap_penalty'])
@click.option('--gap_open', type=int, default=None, help=TRANSLATIONS['en']['gap_open'])
@click.option('--gap_extend', type=int, default=None, help=TRANSLATIONS['en']['gap_extend'])
@click.option('--matrix', default=None, help=TRANSLATIONS['en']['select_matrix'])
@click.option('--subsample', type=int, default=0, help=TRANSLATIONS['en']['subsample'])
@click.option('--preview', is_flag=True, help=TRANSLATIONS['en']['preview_seq'])
@click.option('--verbose', is_flag=True, help=TRANSLATIONS['en']['verbose'])
@click.option('--batch', is_flag=True, help=TRANSLATIONS['en']['batch_mode'])
@click.option('--lang', default='en', type=click.Choice(['en', 'ru']), help=TRANSLATIONS['en']['choose_lang'])
def global_align(input1, input2, directory, output, match, mismatch, gap, gap_open, gap_extend, matrix, subsample,
                 preview, verbose, batch, lang):
    # subcommand для global выравнивания
    tr = TRANSLATIONS[lang]
    params = {
        'mode': 'global', 'input1': input1, 'input2': input2, 'directory': directory, 'output': output,
        'match': match, 'mismatch': mismatch, 'gap': gap, 'gap_open': gap_open, 'gap_extend': gap_extend,
        'matrix': matrix, 'subsample': subsample, 'preview': preview, 'verbose': verbose, 'batch': batch, 'lang': lang
    }
    if batch and not directory:
        console.print(f"{tr['error']} Directory required for batch mode.", style="bold red")
        sys.exit(1)
    if not batch and (not input1 or not input2):
        console.print(f"{tr['error']} {tr['error_pairwise']}", style="bold red")
        sys.exit(1)
    if not validate_params(params, tr):
        sys.exit(1)
    run_alignment(params, tr)


@cli.command()
@click.option('--input1', type=str, help=TRANSLATIONS['en']['select_file1'])
@click.option('--input2', type=str, help=TRANSLATIONS['en']['select_file2'])
@click.option('--directory', type=str, help="Directory for batch alignment")
@click.option('--output', default="alignment.txt", help="Output file")
@click.option('--match', type=int, default=1, help=TRANSLATIONS['en']['match_score'])
@click.option('--mismatch', type=int, default=-1, help=TRANSLATIONS['en']['mismatch_score'])
@click.option('--gap', type=int, default=-2, help=TRANSLATIONS['en']['gap_penalty'])
@click.option('--matrix', default=None, help=TRANSLATIONS['en']['select_matrix'])
@click.option('--subsample', type=int, default=0, help=TRANSLATIONS['en']['subsample'])
@click.option('--preview', is_flag=True, help=TRANSLATIONS['en']['preview_seq'])
@click.option('--verbose', is_flag=True, help=TRANSLATIONS['en']['verbose'])
@click.option('--batch', is_flag=True, help=TRANSLATIONS['en']['batch_mode'])
@click.option('--lang', default='en', type=click.Choice(['en', 'ru']), help=TRANSLATIONS['en']['choose_lang'])
def local(input1, input2, directory, output, match, mismatch, gap, matrix, subsample, preview, verbose, batch, lang):
    # subcommand для local выравнивания
    tr = TRANSLATIONS[lang]
    params = {
        'mode': 'local', 'input1': input1, 'input2': input2, 'directory': directory, 'output': output,
        'match': match, 'mismatch': mismatch, 'gap': gap, 'matrix': matrix, 'subsample': subsample,
        'preview': preview, 'verbose': verbose, 'batch': batch, 'lang': lang
    }
    if batch and not directory:
        console.print(f"{tr['error']} Directory required for batch mode.", style="bold red")
        sys.exit(1)
    if not batch and (not input1 or not input2):
        console.print(f"{tr['error']} {tr['error_pairwise']}", style="bold red")
        sys.exit(1)
    if not validate_params(params, tr):
        sys.exit(1)
    run_alignment(params, tr)


@cli.command()
@click.option('--input1', type=str, help=TRANSLATIONS['en']['select_file1'])
@click.option('--output', default="alignment.txt", help="Output file")
@click.option('--match', type=int, default=1, help=TRANSLATIONS['en']['match_score'])
@click.option('--mismatch', type=int, default=-1, help=TRANSLATIONS['en']['mismatch_score'])
@click.option('--gap', type=int, default=-2, help=TRANSLATIONS['en']['gap_penalty'])
@click.option('--gap_open', type=int, default=None, help=TRANSLATIONS['en']['gap_open'])
@click.option('--gap_extend', type=int, default=None, help=TRANSLATIONS['en']['gap_extend'])
@click.option('--matrix', default=None, help=TRANSLATIONS['en']['select_matrix'])
@click.option('--subsample', type=int, default=0, help=TRANSLATIONS['en']['subsample'])
@click.option('--threads', type=int, default=os.cpu_count(), help=TRANSLATIONS['en']['threads'])
@click.option('--clustal', is_flag=True, help=TRANSLATIONS['en']['clustal'])
@click.option('--preview', is_flag=True, help=TRANSLATIONS['en']['preview_seq'])
@click.option('--verbose', is_flag=True, help=TRANSLATIONS['en']['verbose'])
@click.option('--lang', default='en', type=click.Choice(['en', 'ru']), help=TRANSLATIONS['en']['choose_lang'])
def msa(input1, output, match, mismatch, gap, gap_open, gap_extend, matrix, subsample, threads, clustal, preview,
        verbose, lang):
    # subcommand для msa
    tr = TRANSLATIONS[lang]
    params = {
        'mode': 'msa', 'input1': input1, 'output': output, 'match': match, 'mismatch': mismatch, 'gap': gap,
        'gap_open': gap_open, 'gap_extend': gap_extend, 'matrix': matrix, 'subsample': subsample,
        'threads': threads, 'clustal': clustal, 'preview': preview, 'verbose': verbose, 'lang': lang
    }
    if not input1:
        console.print(f"{tr['error']} {tr['error_msa']}", style="bold red")
        sys.exit(1)
    if not validate_params(params, tr):
        sys.exit(1)
    run_alignment(params, tr)


def format_text_output(params: Dict, tr: Dict, score: Optional[int], stats: Optional[Dict],
                       align1: Optional[str] = None, align2: Optional[str] = None,
                       msa_result: Optional[List[str]] = None) -> str:
    """Форматирует текстовый вывод результатов выравнивания."""
    output_lines = []

    if params['mode'] == 'msa':
        output_lines.append("=" * 50)
        output_lines.append(f"{tr['msa_header']}")
        output_lines.append("=" * 50)
        if msa_result:
            # Для MSA просто выводим последовательности
            if params.get('clustal', False):
                output_lines.append("CLUSTAL format")
                for i in range(0, len(msa_result[0]), 60):
                    for j, seq in enumerate(msa_result):
                        output_lines.append(f"Seq{j + 1} {seq[i:i + 60]}")
                    output_lines.append("")  # Пустая строка между блоками
            else:
                for i, seq in enumerate(msa_result):
                    output_lines.append(f"Seq{i + 1}: {seq}")
    else:
        output_lines.append("=" * 50)
        output_lines.append(f"{tr['pairwise_header']}")
        output_lines.append("=" * 50)
        if align1 and align2:
            # Выравнивание в 3 строки
            match_line = ''.join(['|' if a == b and a != '-' else ' ' for a, b in zip(align1, align2)])
            # Разбиваем на блоки по 60 символов
            for i in range(0, len(align1), 60):
                output_lines.append(align1[i:i + 60])
                output_lines.append(match_line[i:i + 60])
                output_lines.append(align2[i:i + 60])
                output_lines.append("")  # Пустая строка между блоками

    output_lines.append("-" * 30)
    output_lines.append(f"{tr['stats_title']}")
    output_lines.append("-" * 30)

    if score is not None:
        output_lines.append(f"{tr['score_label']}: {score}")
    if stats:
        output_lines.append(f"{tr['identity_label']}: {stats['identity']:.2f}%")
        output_lines.append(f"{tr['gaps_label']}: {stats['gaps']}")

    return "\n".join(output_lines)


def run_alignment(params: Dict, tr: Dict):
    # выполнение выравнивания
    if params['verbose']:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    start_time = time.time()
    memory_usage = 0
    if psutil:
        process = psutil.Process(os.getpid())
        start_mem = process.memory_info().rss / 1024 ** 2

    score = None
    stats = None
    align1, align2 = None, None
    msa_result = None

    if params.get('batch', False):
        if 'directory' not in params:
            console.print(f"{tr['error']} Directory is required for batch mode.", style="bold red")
            sys.exit(1)
        result = run_batch_alignment(params['directory'], params, tr)
    else:
        scoring_matrix = load_scoring_matrix(params['matrix']) if params['matrix'] else None
        try:
            sequences = load_sequences(params['input1'])
        except Exception as e:
            console.print(f"{tr['error']} Failed to load sequences from {params['input1']}: {e}", style="bold red")
            sys.exit(1)

        if params['subsample'] > 0:
            sequences = [seq[:params['subsample']] for seq in sequences]
            console.print(tr['subsampled'].format(params['subsample']), style="yellow")

        if params.get('preview', False):
            preview_content = "\n".join([f"Seq{i + 1}: {seq[:100]}" for i, seq in enumerate(sequences)])
            console.print(Panel(preview_content, title=tr['preview_title'], style="blue"))

        with Progress() as progress:
            task = progress.add_task(tr['processing'], total=100)
            if params['mode'] == 'msa':
                if len(sequences) < 2:
                    console.print(f"{tr['error']} {tr['error_msa']}", style="bold red")
                    sys.exit(1)
                try:
                    msa_result = multiple_sequence_alignment(
                        sequences, params['match'], params['mismatch'], params['gap'],
                        params.get('gap_open'), params.get('gap_extend'), scoring_matrix,
                        params.get('threads', os.cpu_count())
                    )
                    result = format_msa(msa_result, clustal=params.get('clustal', False))
                    # Для текстового вывода MSA
                    text_result = format_text_output(params, tr, None, None, msa_result=msa_result)
                except Exception as e:
                    console.print(f"{tr['error']} MSA failed: {e}", style="bold red")
                    sys.exit(1)
            else:
                if 'input2' not in params:
                    console.print(f"{tr['error']} {tr['error_pairwise']}", style="bold red")
                    sys.exit(1)
                try:
                    seq1 = sequences[0]
                    seqs2 = load_sequences(params['input2'])
                    if not seqs2:
                        console.print(f"{tr['error']} No sequences found in {params['input2']}.", style="bold red")
                        sys.exit(1)
                    seq2 = seqs2[0]
                except Exception as e:
                    console.print(f"{tr['error']} Failed to load sequence from {params['input2']}: {e}",
                                  style="bold red")
                    sys.exit(1)

                if params['subsample'] > 0:
                    seq1 = seq1[:params['subsample']]
                    seq2 = seq2[:params['subsample']]
                    console.print(tr['subsampled'].format(params['subsample']), style="yellow")

                try:
                    if params['mode'] == 'global':
                        align1, align2, score = needleman_wunsch(
                            seq1, seq2, params['match'], params['mismatch'], params['gap'],
                            params.get('gap_open'), params.get('gap_extend'), scoring_matrix
                        )
                    else:
                        align1, align2, score = smith_waterman(
                            seq1, seq2, params['match'], params['mismatch'], params['gap'], scoring_matrix
                        )
                    stats = compute_stats(align1, align2)
                    print_alignment_table(align1, align2, tr)
                    # Для текстового вывода pairwise
                    text_result = format_text_output(params, tr, score, stats, align1, align2)
                    result = f"Score: {score}\n"
                    result += f"{tr['identity']}: {stats['identity']:.2f}%\n{tr['gaps']}: {stats['gaps']}\n"
                except Exception as e:
                    console.print(f"{tr['error']} Alignment failed: {e}", style="bold red")
                    sys.exit(1)
            progress.update(task, advance=100)

    end_time = time.time()
    if psutil:
        memory_usage = process.memory_info().rss / 1024 ** 2 - start_mem

    # Добавляем статистику в конец результата
    stats_line = f"\n{tr['time_label']}: {end_time - start_time:.2f} {tr['sec']}\n{tr['memory_label']}: {memory_usage:.2f} {tr['mb']}"
    result += stats_line
    if 'text_result' in locals():
        text_result += stats_line
    else:
        text_result = stats_line

    try:
        with open(params['output'], "w") as f:
            f.write(result)  # Основной результат в формате Rich/Clustal
        # Сохраняем также текстовый вывод в отдельный файл
        text_output_path = params['output'].replace(".txt", "_text.txt")
        with open(text_output_path, "w") as f:
            f.write(text_result)
        console.print(Panel(result, title=tr['success'], style="bold green"))
    except Exception as e:
        console.print(f"{tr['error']} Failed to write output to {params['output']}: {e}", style="bold red")
        sys.exit(1)


if __name__ == "__main__":
    cli()
