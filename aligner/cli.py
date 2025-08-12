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

# зависимости pip install click rich inquirer pyyaml biopython numpy numba psutil pexpect
console = Console()

# ascii-арт genomic aligner
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

  ______  __ __         |                                        
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
        'mode_global': "Global (NW): Full alignment of sequences.",
        'mode_local': "Local (SW): Align best parts.",
        'mode_msa': "MSA: Align multiple sequences.",
        'select_dir': "Select directory:",
        'select_file1': "Select input1 FASTA:",
        'select_file2': "Select input2 FASTA:",
        'select_matrix': "Scoring matrix (None default):",
        'match_score': "Match score (1 default):",
        'desc_match_score': "Rewards matching symbols, higher for important matches.",
        'mismatch_score': "Mismatch score (-1 default):",
        'desc_mismatch_score': "Penalizes mismatches, lower for strict alignment.",
        'gap_penalty': "Gap penalty (-2 default):",
        'desc_gap_penalty': "Penalizes gaps in sequences.",
        'gap_open': "Gap open (None default):",
        'desc_gap_open': "Cost to start a gap in affine mode.",
        'gap_extend': "Gap extend (None default):",
        'desc_gap_extend': "Cost to extend a gap in affine mode.",
        'subsample': "Subsample N bases (0 default):",
        'desc_subsample': "Limit sequences to N bases for testing large files.",
        'threads': "Threads for MSA:",
        'desc_threads': "Number of threads for parallel computation in MSA.",
        'clustal': "Clustal format? (y/n)",
        'desc_clustal': "Output MSA in Clustal format for compatibility.",
        'verbose': "Verbose logging? (y/n)",
        'desc_verbose': "Show detailed logs during execution.",
        'preview_seq': "Preview 100 bp? (y/n)",
        'desc_preview_seq': "Show first 100 bases of sequences.",
        'tutorial': "Run tutorial? (y/n)",
        'desc_tutorial': "Show examples and param explanations.",
        'batch_mode': "Batch mode? (y/n)",
        'desc_batch_mode': "Align all vs all in directory.",
        'config': "Load config (path):",
        'processing': "Processing...",
        'success': "Success!",
        'error': "Error:",
        'error_msa': "MSA needs 2+ sequences.",
        'error_pairwise': "Pairwise needs two files.",
        'error_no_files': "No FASTA files.",
        'error_invalid_file': "Invalid file.",
        'error_negative': "Penalties negative.",
        'error_config': "Invalid config path.",
        'confirm_run': "Run? (y/n)",
        'change_dir': "Change dir? (y/n)",
        'time': "Time",
        'memory': "Memory",
        'sec': "sec",
        'mb': "MB",
        'identity': "Identity %",
        'gaps': "Gaps count",
        'config_load': "Load config? (path/n):",
        'config_save': "Save params? (y/n)",
        'config_saved': "Saved to {path}"
    },
    'ru': {
        'welcome': "Добро пожаловать!",
        'choose_lang': "Язык (en/ru):",
        'choose_mode': "Режим:",
        'mode_global': "Global (NW): Полное выравнивание.",
        'mode_local': "Local (SW): Лучшие части.",
        'mode_msa': "MSA: Множественное.",
        'select_dir': "Директория:",
        'select_file1': "Input1 FASTA:",
        'select_file2': "Input2 FASTA:",
        'select_matrix': "Matrix (None default):",
        'match_score': "Match (1 default):",
        'desc_match_score': "Награда за совпадения, выше для важных.",
        'mismatch_score': "Mismatch (-1 default):",
        'desc_mismatch_score': "Штраф за несовпадения, ниже для строгого.",
        'gap_penalty': "Gap (-2 default):",
        'desc_gap_penalty': "Штраф за пробелы в последовательностях.",
        'gap_open': "Gap open (None default):",
        'desc_gap_open': "Стоимость открытия gap в affine.",
        'gap_extend': "Gap extend (None default):",
        'desc_gap_extend': "Стоимость расширения gap в affine.",
        'subsample': "Subsample N (0 default):",
        'desc_subsample': "Ограничить N баз для теста больших файлов.",
        'threads': "Потоки для MSA:",
        'desc_threads': "Количество потоков для параллели в MSA.",
        'clustal': "Clustal? (y/n)",
        'desc_clustal': "Вывод MSA в Clustal для совместимости.",
        'verbose': "Verbose? (y/n)",
        'desc_verbose': "Детальные логи во время работы.",
        'preview_seq': "Preview 100 bp? (y/n)",
        'desc_preview_seq': "Показать первые 100 баз последовательностей.",
        'tutorial': "Tutorial? (y/n)",
        'desc_tutorial': "Показать примеры и объяснения параметров.",
        'batch_mode': "Batch? (y/n)",
        'desc_batch_mode': "Выравнивание все vs все в директории.",
        'config': "Load config (path):",
        'processing': "Обработка...",
        'success': "Успех!",
        'error': "Ошибка:",
        'error_msa': "MSA 2+ seq.",
        'error_pairwise': "Pairwise два файла.",
        'error_no_files': "Нет FASTA.",
        'error_invalid_file': "Неверный файл.",
        'error_negative': "Штрафы отрицательны.",
        'error_config': "Неверный path config.",
        'confirm_run': "Запустить? (y/n)",
        'change_dir': "Сменить dir? (y/n)",
        'time': "Время",
        'memory': "Память",
        'sec': "сек",
        'mb': "МБ",
        'identity': "Идентичность %",
        'gaps': "Gaps",
        'config_load': "Load config? (path/n):",
        'config_save': "Save? (y/n)",
        'config_saved': "Сохранено {path}"
    }
}


def get_fasta_files(directory: str) -> List[str]:
    # список fasta/gz в dir
    return [f for f in os.listdir(directory) if f.endswith(('.fasta', '.fa', '.gz'))]


def validate_file(file_path: str, tr: Dict) -> bool:
    # проверка файла
    if not os.path.exists(file_path):
        console.print(f"{tr['error']} {file_path} not exist.", style="bold red")
        return False
    if not file_path.endswith(('.fasta', '.fa', '.gz')):
        console.print(f"{tr['error']} {tr['error_invalid_file']}", style="bold red")
        return False
    return True


def validate_params(params: Dict, tr: Dict) -> bool:
    # штрафы отрицательны
    for param in ['gap', 'mismatch', 'gap_open', 'gap_extend']:
        if param in params and params[param] is not None and params[param] > 0:
            console.print(f"{tr['error']} {tr['error_negative']}", style="bold red")
            return False
    return True


def interactive_wizard() -> tuple[Dict, Dict]:
    # wizard для params
    params = {}
    lang = inquirer.prompt([inquirer.List('lang', message=TRANSLATIONS['en']['choose_lang'], choices=['en', 'ru'])])[
        'lang']
    tr = TRANSLATIONS[lang]
    console.print(Panel(tr['welcome'], style="bold blue"))

    mode_choices = [
        (tr['mode_global'], 'global'),
        (tr['mode_local'], 'local'),
        (tr['mode_msa'], 'msa')
    ]
    params['mode'] = inquirer.prompt([inquirer.List('mode', message=tr['choose_mode'], choices=mode_choices)])['mode']

    params['batch'] = inquirer.prompt([inquirer.Confirm('batch', message=tr['batch_mode'], default=False)])['batch']

    directory = os.getcwd()
    while True:
        directory = inquirer.prompt(
            [inquirer.Path('dir', message=tr['select_dir'], path_type=inquirer.Path.DIRECTORY, default=directory)])[
            'dir']
        fasta_files = get_fasta_files(directory)
        if not fasta_files and not params['batch']:
            console.print(f"{tr['error']} {tr['error_no_files']}", style="bold red")
            change = inquirer.prompt([inquirer.Confirm('change', message=tr['change_dir'], default=True)])['change']
            if not change:
                sys.exit(1)
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

    params['matrix'] = inquirer.prompt([inquirer.Text('matrix', message=tr['select_matrix'], default='None')])['matrix']
    if params['matrix'] == 'None':
        params['matrix'] = None
    try:
        params['match'] = int(
            inquirer.prompt([inquirer.Text('match', message=tr['match_score'], default='1')])['match'])
        params['mismatch'] = int(
            inquirer.prompt([inquirer.Text('mismatch', message=tr['mismatch_score'], default='-1')])['mismatch'])
        params['gap'] = int(inquirer.prompt([inquirer.Text('gap', message=tr['gap_penalty'], default='-2')])['gap'])
        gap_open_input = inquirer.prompt([inquirer.Text('gap_open', message=tr['gap_open'], default='None')])[
            'gap_open']
        params['gap_open'] = int(gap_open_input) if gap_open_input != 'None' else None
        gap_extend_input = inquirer.prompt([inquirer.Text('gap_extend', message=tr['gap_extend'], default='None')])[
            'gap_extend']
        params['gap_extend'] = int(gap_extend_input) if gap_extend_input != 'None' else None
        params['subsample'] = int(
            inquirer.prompt([inquirer.Text('subsample', message=tr['subsample'], default='0')])['subsample'])
    except ValueError:
        console.print("Ошибка: введите число.", style="bold red")
        sys.exit(1)

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

    while True:
        config_load = inquirer.prompt([inquirer.Text('config_load', message=tr['config_load'], default='n')])[
            'config_load']
        if config_load == 'n':
            break
        if os.path.isfile(config_load):
            try:
                params = load_config(config_load)
                params['lang'] = lang
                tr = TRANSLATIONS[lang]
                break
            except (yaml.YAMLError, PermissionError) as e:
                console.print(f"{tr['error']} Invalid config: {str(e)}.", style="bold red")
        else:
            console.print(f"{tr['error']} {tr['error_config']}", style="bold red")
        retry = inquirer.prompt([inquirer.Confirm('retry', message="Retry? (y/n)", default=True)])['retry']
        if not retry:
            break

    if inquirer.prompt([inquirer.Confirm('config_save', message=tr['config_save'], default=False)])['config_save']:
        save_config(params, 'config.yaml', tr)

    if not validate_params(params, tr):
        sys.exit(1)

    confirm = inquirer.prompt([inquirer.Confirm('confirm', message=tr['confirm_run'], default=True)])['confirm']
    if not confirm:
        sys.exit(0)

    return params, tr


def run_tutorial(tr: Dict):
    # tutorial с примерами и описаниями params
    param_desc = (
        f"{tr['match_score']}: {tr['desc_match_score']}\n"
        f"{tr['mismatch_score']}: {tr['desc_mismatch_score']}\n"
        f"{tr['gap_penalty']}: {tr['desc_gap_penalty']}\n"
        f"{tr['gap_open']}: {tr['desc_gap_open']}\n"
        f"{tr['gap_extend']}: {tr['desc_gap_extend']}\n"
        f"{tr['subsample']}: {tr['desc_subsample']}\n"
        f"{tr['threads']}: {tr['desc_threads']}\n"
        f"{tr['clustal']}: {tr['desc_clustal']}\n"
        f"{tr['verbose']}: {tr['desc_verbose']}\n"
        f"{tr['preview_seq']}: {tr['desc_preview_seq']}\n"
        f"{tr['batch_mode']}: {tr['desc_batch_mode']}"
    )
    console.print(Panel(param_desc, title="Param descriptions", style="bold yellow"))
    console.print(Panel(
        "Example 1: Global\n"
        "Seq: 'AGC' 'ACGC'\n"
        "Score: 1\n"
        "A-GC\n"
        " | |\n"
        "ACGC\n\n"
        "Example 2: Local\n"
        "Seq: 'ATGC' 'TGCA'\n"
        "Score: 3\n"
        "TGC\n"
        "|||\n"
        "TGC\n\n"
        "Example 3: MSA\n"
        "Seq1: A-GC\n"
        "Seq2: ACGC\n"
        "Seq3: AGGC",
        title=tr['tutorial'], style="bold blue"
    ))


def load_config(path: str) -> Dict:
    # load yaml config
    with open(path, 'r') as f:
        return yaml.safe_load(f)


def save_config(params: Dict, path: str, tr: Dict):
    # save params to yaml
    with open(path, 'w') as f:
        yaml.dump(params, f)
    console.print(f"{tr['config_saved'].format(path=path)}", style="green")


def print_alignment_table(align1: str, align2: str, tr: Dict):
    # таблица выравнивания с цветами
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
    # статистика идентичность gaps
    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    length = len(align1)
    gaps = align1.count('-') + align2.count('-')
    identity = (matches / length * 100) if length > 0 else 0
    return {'identity': identity, 'gaps': gaps}


def run_batch_alignment(directory: str, params: Dict, tr: Dict) -> str:
    # batch pairwise all vs all
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
                console.print(f"Processing: {file1} vs {file2}", style="bold blue")
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
    # дефолт wizard с арт
    console.print(ASCII_ART)
    if config:
        if os.path.isfile(config):
            try:
                params = load_config(config)
                tr = TRANSLATIONS[params.get('lang', 'en')]
                if not validate_params(params, tr):
                    sys.exit(1)
                if params.get('tutorial', False):
                    run_tutorial(tr)
                run_alignment(params, tr)
            except (yaml.YAMLError, PermissionError) as e:
                console.print(f"{tr['error']} Invalid config: {str(e)}.", style="bold red")
                sys.exit(1)
        else:
            console.print("Error: Config not file.", style="bold red")
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
@click.option('--directory', type=str, help="Dir for batch")
@click.option('--output', default="alignment.txt", help="Output")
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
    # global subcommand
    tr = TRANSLATIONS[lang]
    params = {
        'mode': 'global', 'input1': input1, 'input2': input2, 'directory': directory, 'output': output,
        'match': match, 'mismatch': mismatch, 'gap': gap, 'gap_open': gap_open, 'gap_extend': gap_extend,
        'matrix': matrix, 'subsample': subsample, 'preview': preview, 'verbose': verbose, 'batch': batch, 'lang': lang
    }
    if batch and not directory:
        console.print(f"{tr['error']} Dir for batch.", style="bold red")
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
@click.option('--directory', type=str, help="Dir for batch")
@click.option('--output', default="alignment.txt", help="Output")
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
    # local subcommand
    tr = TRANSLATIONS[lang]
    params = {
        'mode': 'local', 'input1': input1, 'input2': input2, 'directory': directory, 'output': output,
        'match': match, 'mismatch': mismatch, 'gap': gap, 'matrix': matrix, 'subsample': subsample,
        'preview': preview, 'verbose': verbose, 'batch': batch, 'lang': lang
    }
    if batch and not directory:
        console.print(f"{tr['error']} Dir for batch.", style="bold red")
        sys.exit(1)
    if not batch and (not input1 or not input2):
        console.print(f"{tr['error']} {tr['error_pairwise']}", style="bold red")
        sys.exit(1)
    if not validate_params(params, tr):
        sys.exit(1)
    run_alignment(params, tr)


@cli.command()
@click.option('--input1', type=str, help=TRANSLATIONS['en']['select_file1'])
@click.option('--output', default="alignment.txt", help="Output")
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
    # msa subcommand
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


@cli.command(name='help')
@click.option('--lang', default='en', type=click.Choice(['en', 'ru']), help=TRANSLATIONS['en']['choose_lang'])
def help_cmd(lang):
    # help subcommand for param desc
    tr = TRANSLATIONS[lang]
    desc = (
        f"{tr['match_score']}: {tr['desc_match_score']}\n"
        f"{tr['mismatch_score']}: {tr['desc_mismatch_score']}\n"
        f"{tr['gap_penalty']}: {tr['desc_gap_penalty']}\n"
        f"{tr['gap_open']}: {tr['desc_gap_open']}\n"
        f"{tr['gap_extend']}: {tr['desc_gap_extend']}\n"
        f"{tr['subsample']}: {tr['desc_subsample']}\n"
        f"{tr['threads']}: {tr['desc_threads']}\n"
        f"{tr['clustal']}: {tr['desc_clustal']}\n"
        f"{tr['verbose']}: {tr['desc_verbose']}\n"
        f"{tr['preview_seq']}: {tr['desc_preview_seq']}\n"
        f"{tr['tutorial']}: {tr['desc_tutorial']}\n"
        f"{tr['batch_mode']}: {tr['desc_batch_mode']}"
    )
    console.print(Panel(desc, title="Parameters help", style="bold yellow"))


def run_alignment(params: Dict, tr: Dict):
    # запуск выравнивания
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
            console.print(f"Subsampled to {params['subsample']} bases.", style="yellow")

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
                    console.print(f"Subsampled to {params['subsample']} bases.", style="yellow")
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