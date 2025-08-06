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

# зависимости: pip install click rich inquirer pyyaml biopython numpy numba psutil
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
  ______                                                                 
 /      \                                                                
|  $$$$$$\ ______   ______  __    __  ______  _______   _______  ______  
| $$___\$$/      \ /      \|  \  |  \/      \|       \ /       \/      \ 
 \$$    \|  $$$$$$|  $$$$$$| $$  | $|  $$$$$$| $$$$$$$|  $$$$$$|  $$$$$$\
 _\$$$$$$| $$    $| $$  | $| $$  | $| $$    $| $$  | $| $$     | $$    $$
|  \__| $| $$$$$$$| $$__| $| $$__/ $| $$$$$$$| $$  | $| $$_____| $$$$$$$$
 \$$    $$\$$     \\$$    $$\$$    $$\$$     | $$  | $$\$$     \\$$     \
  \$$$$$$  \$$$$$$$ \$$$$$$$ \$$$$$$  \$$$$$$$\$$   \$$ \$$$$$$$ \$$$$$$$
                        | $$                                             
  ______  __ __         | $$                                             
 /      \|  |  \         \$$                                             
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
        'mode_global': "Global (Needleman-Wunsch): Full alignment of entire sequences, suitable for similar genes.",
        'mode_local': "Local (Smith-Waterman): Best subsequences only, ideal for motif or domain search.",
        'mode_msa': "MSA: Multiple sequence alignment for 2+ sequences, for phylogenetic or structural analysis.",
        'select_dir': "Select directory (default: current):",
        'select_file1': "Select FASTA file for input1:",
        'select_file2': "Select FASTA file for input2 (pairwise only):",
        'select_matrix': "Scoring matrix (e.g., BLOSUM62 for proteins, None for default DNA scoring):",
        'match_score': "Match score (default 1): Value for matching characters. Recommend 1 for DNA, 5 for proteins.",
        'mismatch_score': "Mismatch score (default -1): Penalty for mismatches. Recommend -1 for DNA, -4 for proteins.",
        'gap_penalty': "Gap penalty (default -2): Penalty for gaps. Use negative values.",
        'gap_open': "Gap open penalty (affine, default None): Penalty for starting a gap. Set to enable affine gaps.",
        'gap_extend': "Gap extend penalty (affine, default None): Penalty for extending a gap.",
        'subsample': "Subsample first N bases (0 for full): For large files to speed up testing.",
        'threads': "Number of threads for MSA (default cpu_count):",
        'clustal': "Output MSA in Clustal format? (y/n)",
        'verbose': "Enable verbose logging for detailed steps? (y/n)",
        'preview_seq': "Preview first 100 bases of sequences? (y/n)",
        'tutorial': "Run tutorial with example alignments? (y/n)",
        'batch_mode': "Run batch alignment for all FASTA files in directory? (y/n)",
        'config': "Load configuration from YAML file (path):",
        'processing': "Processing alignment...",
        'success': "Alignment completed successfully!",
        'error': "Error:",
        'error_msa': "MSA requires a FASTA file with at least 2 sequences.",
        'error_pairwise': "Pairwise alignment requires two FASTA files.",
        'error_no_files': "No FASTA files found in directory. Please check or provide files.",
        'error_invalid_file': "Invalid file selected. Please choose a .fasta or .gz file.",
        'error_negative': "Penalty values (gap, mismatch, gap_open, gap_extend) must be negative.",
        'time': "Time",
        'memory': "Memory",
        'sec': "sec",
        'mb': "MB",
        'identity': "Identity %",
        'gaps': "Gaps count",
        'config_load': "Load config from YAML file? (enter path or 'n' for none):",
        'config_save': "Save current parameters to YAML file? (y/n)",
        'config_saved': "Configuration saved to {path}"
    },
    'ru': {
        'welcome': "Добро пожаловать в Aligner CLI!",
        'choose_lang': "Выберите язык (en/ru):",
        'choose_mode': "Выберите режим:",
        'mode_global': "Global (Needleman-Wunsch): Полное выравнивание всей последовательности, подходит для похожих генов.",
        'mode_local': "Local (Smith-Waterman): Только лучшие субпоследовательности, для поиска мотивов или доменов.",
        'mode_msa': "MSA: Множественное выравнивание для 2+ последовательностей, для филогенетики или структурного анализа.",
        'select_dir': "Выберите директорию (default: текущая):",
        'select_file1': "Выберите FASTA-файл для input1:",
        'select_file2': "Выберите FASTA-файл для input2 (только для pairwise):",
        'select_matrix': "Scoring matrix (e.g., BLOSUM62 для белков, None для ДНК по умолчанию):",
        'match_score': "Score за совпадение (default 1): Значение за совпадающие символы. Рекомендуется 1 для ДНК, 5 для белков.",
        'mismatch_score': "Score за несовпадение (default -1): Штраф за несовпадения. Рекомендуется -1 для ДНК, -4 для белков.",
        'gap_penalty': "Штраф за gap (default -2): Штраф за пробелы. Используйте отрицательные значения.",
        'gap_open': "Штраф за открытие gap (affine, default None): Штраф за начало пробела. Установите для включения affine gaps.",
        'gap_extend': "Штраф за расширение gap (affine, default None): Штраф за продолжение пробела.",
        'subsample': "Subsample первых N баз (0 для полного): Для больших файлов для ускорения тестирования.",
        'threads': "Количество потоков для MSA (default cpu_count):",
        'clustal': "Вывести MSA в формате Clustal? (y/n)",
        'verbose': "Включить детальный logging для подробных шагов? (y/n)",
        'preview_seq': "Предпросмотр первых 100 баз последовательностей? (y/n)",
        'tutorial': "Запустить tutorial с примерами выравниваний? (y/n)",
        'batch_mode': "Запустить batch-выравнивание для всех FASTA в директории? (y/n)",
        'config': "Загрузить конфигурацию из YAML файла (путь):",
        'processing': "Обработка выравнивания...",
        'success': "Выравнивание успешно завершено!",
        'error': "Ошибка:",
        'error_msa': "Для MSA нужен FASTA-файл с минимум 2 последовательностями.",
        'error_pairwise': "Для pairwise-выравнивания нужны два FASTA-файла.",
        'error_no_files': "FASTA-файлы не найдены в директории. Проверьте или предоставьте файлы.",
        'error_invalid_file': "Выбран неверный файл. Выберите файл .fasta или .gz.",
        'error_negative': "Штрафы (gap, mismatch, gap_open, gap_extend) должны быть отрицательными.",
        'time': "Время",
        'memory': "Память",
        'sec': "сек",
        'mb': "МБ",
        'identity': "Идентичность %",
        'gaps': "Количество gaps",
        'config_load': "Загрузить конфигурацию из YAML? (введите путь или 'n' для пропуска):",
        'config_save': "Сохранить текущие параметры в YAML? (y/n)",
        'config_saved': "Конфигурация сохранена в {path}"
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
        if param in params and params[param] is not None and params[param] > 0:
            console.print(f"{tr['error']} {tr['error_negative']}", style="bold red")
            return False
    return True


def interactive_wizard() -> tuple[Dict, Dict]:
    # интерактивный wizard для параметров
    params = {}
    lang = inquirer.prompt([inquirer.List('lang', message=TRANSLATIONS['en']['choose_lang'], choices=['en', 'ru'])])[
        'lang']
    tr = TRANSLATIONS[lang]
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
        params['directory'] = directory
        params['input_files'] = [os.path.join(directory, f) for f in fasta_files]

    # параметры выравнивания
    params['matrix'] = inquirer.prompt([inquirer.Text('matrix', message=tr['select_matrix'], default='None')])[
                           'matrix'] or None
    params['match'] = int(inquirer.prompt([inquirer.Text('match', message=tr['match_score'], default='1')])['match'])
    params['mismatch'] = int(
        inquirer.prompt([inquirer.Text('mismatch', message=tr['mismatch_score'], default='-1')])['mismatch'])
    params['gap'] = int(inquirer.prompt([inquirer.Text('gap', message=tr['gap_penalty'], default='-2')])['gap'])
    params['gap_open'] = inquirer.prompt([inquirer.Text('gap_open', message=tr['gap_open'], default='None')])[
                             'gap_open'] or None
    if params['gap_open'] != 'None':
        params['gap_open'] = int(params['gap_open'])
    params['gap_extend'] = inquirer.prompt([inquirer.Text('gap_extend', message=tr['gap_extend'], default='None')])[
                               'gap_extend'] or None
    if params['gap_extend'] != 'None':
        params['gap_extend'] = int(params['gap_extend'])
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
    params['output'] = inquirer.prompt([inquirer.Text('output', message="Output file:", default="alignment.txt")])[
        'output']

    # config load/save
    config_load = inquirer.prompt([inquirer.Text('config_load', message=tr['config_load'], default='n')])['config_load']
    if config_load != 'n':
        if os.path.exists(config_load):
            params = load_config(config_load)
            params['lang'] = lang  # сохраняем язык
            tr = TRANSLATIONS[lang]
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
        "ACGC\n\n"
        "Example 2: Local alignment (Smith-Waterman)\n"
        "Sequences: 'ATGC' and 'TGCA'\n"
        "Score: 3\n"
        "TGC\n"
        "|||\n"
        "TGC\n\n"
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
    console.print(f"{TRANSLATIONS[params['lang']]['config_saved'].format(path=path)}", style="green")


def print_alignment_table(align1: str, align2: str, tr: Dict):
    # таблица для pairwise выравнивания с цветами
    table = Table(title="Alignment")
    table.add_column("Seq1", style="cyan")
    table.add_column("Matches", style="magenta")
    table.add_column("Seq2", style="green")
    match_line = ''.join(
        ['[green]|[/green]' if a == b and a != '-' else '[red] [/red]' for a, b in zip(align1, align2)])
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
            seq1 = load_sequences(os.path.join(directory, file1))[0]
            if params['subsample'] > 0:
                seq1 = seq1[:params['subsample']]
            for j, file2 in enumerate(fasta_files[i + 1:], start=i + 1):
                seq2 = load_sequences(os.path.join(directory, file2))[0]
                if params['subsample'] > 0:
                    seq2 = seq2[:params['subsample']]
                console.print(f"\nProcessing: {file1} vs {file2}", style="bold blue")
                if params['mode'] == 'global':
                    align1, align2, score = needleman_wunsch(
                        seq1, seq2, params['match'], params['mismatch'], params['gap'],
                        params.get('gap_open'), params.get('gap_extend'), scoring_matrix
                    )
                else:
                    align1, align2, score = smith_waterman(
                        seq1, seq2, params['match'], params['mismatch'], params['gap'], scoring_matrix
                    )
                result += f"\nAlignment: {file1} vs {file2}\nScore: {score}\n"
                print_alignment_table(align1, align2, tr)
                stats = compute_stats(align1, align2)
                result += f"{tr['identity']}: {stats['identity']:.2f}%\n{tr['gaps']}: {stats['gaps']}\n"
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
            tr = TRANSLATIONS[params.get('lang', 'en')]
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
    # subcommand для global выравнивания (переименовано из 'global' во избежание конфликта с ключевым словом)
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

    if params.get('batch', False):
        result = run_batch_alignment(params['directory'], params, tr)
    else:
        scoring_matrix = load_scoring_matrix(params['matrix']) if params['matrix'] else None
        sequences = load_sequences(params['input1'])
        if params['subsample'] > 0:
            sequences = [seq[:params['subsample']] for seq in sequences]
            console.print(tr['subsampled'].format(params['subsample']), style="yellow")

        if params.get('preview', False):
            console.print(
                Panel("\n".join([f"Seq{i + 1}: {seq[:100]}" for i, seq in enumerate(sequences)]), title="Preview",
                      style="blue"))

        with Progress() as progress:
            task = progress.add_task(tr['processing'], total=100)

            if params['mode'] == 'msa':
                if len(sequences) < 2:
                    console.print(f"{tr['error']} {tr['error_msa']}", style="bold red")
                    sys.exit(1)
                aligned = multiple_sequence_alignment(
                    sequences, params['match'], params['mismatch'], params['gap'],
                    params.get('gap_open'), params.get('gap_extend'), scoring_matrix,
                    params.get('threads', os.cpu_count())
                )
                result = format_msa(aligned, clustal=params.get('clustal', False))
            else:
                if 'input2' not in params:
                    console.print(f"{tr['error']} {tr['error_pairwise']}", style="bold red")
                    sys.exit(1)
                seq1 = sequences[0]
                seq2 = load_sequences(params['input2'])[0]
                if params['subsample'] > 0:
                    seq1 = seq1[:params['subsample']]
                    seq2 = seq2[:params['subsample']]
                    console.print(tr['subsampled'].format(params['subsample']), style="yellow")
                if params['mode'] == 'global':
                    align1, align2, score = needleman_wunsch(
                        seq1, seq2, params['match'], params['mismatch'], params['gap'],
                        params.get('gap_open'), params.get('gap_extend'), scoring_matrix
                    )
                else:
                    align1, align2, score = smith_waterman(
                        seq1, seq2, params['match'], params['mismatch'], params['gap'], scoring_matrix
                    )
                result = f"Score: {score}\n"
                print_alignment_table(align1, align2, tr)
                stats = compute_stats(align1, align2)
                result += f"{tr['identity']}: {stats['identity']:.2f}%\n{tr['gaps']}: {stats['gaps']}\n"

            progress.update(task, advance=100)

    end_time = time.time()
    if psutil:
        memory_usage = process.memory_info().rss / 1024 ** 2 - start_mem

    result += f"\n{tr['time']}: {end_time - start_time:.2f} {tr['sec']}\n{tr['memory']}: {memory_usage:.2f} {tr['mb']}"
    with open(params['output'], "w") as f:
        f.write(result)
    console.print(Panel(result, title=tr['success'], style="bold green"))


if __name__ == "__main__":
    cli()