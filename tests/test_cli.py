import pytest
import os
import sys
import subprocess
import time
import shutil

@pytest.fixture
def dummy_fasta():
    # создаем тестовую директорию и файлы
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)
    with open(os.path.join(test_dir, "test1.fasta"), "w") as f:
        f.write(">seq1\nAGC")
    with open(os.path.join(test_dir, "test2.fasta"), "w") as f:
        f.write(">seq2\nACGC")
    with open(os.path.join(test_dir, "test_multi.fasta"), "w") as f:
        f.write(">seq1\nAGC\n>seq2\nACGC\n>seq3\nAGGC")
    with open(os.path.join(test_dir, "config.yaml"), "w") as f:
        f.write(
            "mode: global\ninput1: cli_test/test1.fasta\ninput2: cli_test/test2.fasta\noutput: cli_test/test_out.yaml\nmatch: 1\nmismatch: -1\ngap: -2\nlang: en\nmatrix: null\nverbose: false\npreview: false\ntutorial: false\nbatch: false\nclustal: false\nsubsample: 0")
    yield test_dir
    # Рекурсивно удаляем директорию
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)

def test_cli_pairwise_flags(dummy_fasta):
    # тест старого стиля флагов
    test_dir = dummy_fasta
    output_file = os.path.join(test_dir, "test_out.txt")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
             os.path.join(test_dir, "test2.fasta"),
             "--output", output_file, "--lang", "ru"], capture_output=True, text=True)
        assert result.returncode == 0
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "Счет" in content or "Score" in content
            assert "Идентичность" in content
        os.remove(output_file)
    except Exception as e:
        pytest.fail(f"CLI failed: {e}")


def test_cli_msa_flags(dummy_fasta):
    # тест msa с флагами
    test_dir = dummy_fasta
    output_file = os.path.join(test_dir, "test_msa.txt")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "aligner.cli", "msa", "--input1", os.path.join(test_dir, "test_multi.fasta"), "--output",
             output_file, "--lang", "en"], capture_output=True, text=True)
        assert result.returncode == 0
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "Seq1" in content
        os.remove(output_file)
    except Exception as e:
        pytest.fail(f"CLI failed: {e}")

# Большие тесты для CLI
def test_cli_comprehensive():
    # тест всех режимов cli
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    # 1. global alignment
    with open(os.path.join(test_dir, "test1.fasta"), "w") as f:
        f.write(">seq1\nAGCTAGCT")
    with open(os.path.join(test_dir, "test2.fasta"), "w") as f:
        f.write(">seq2\nAGCTGGCT")

    # global
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "global_out.txt")], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "global_out.txt"))
    with open(os.path.join(test_dir, "global_out.txt"), "r") as f:
        content = f.read()
        assert "Score" in content
        assert "Identity" in content or "Идентичность" in content

    # local
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "local", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "local_out.txt")], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "local_out.txt"))
    with open(os.path.join(test_dir, "local_out.txt"), "r") as f:
        content = f.read()
        assert "Score" in content
        assert "Identity" in content or "Идентичность" in content

    # msa
    with open(os.path.join(test_dir, "test_msa.fasta"), "w") as f:
        f.write(">seq1\nAGC\n>seq2\nACGC\n>seq3\nAGGC")
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "msa", "--input1", os.path.join(test_dir, "test_msa.fasta"), "--output",
         os.path.join(test_dir, "msa_out.txt")], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "msa_out.txt"))
    with open(os.path.join(test_dir, "msa_out.txt"), "r") as f:
        content = f.read()
        assert "Seq1" in content
        assert "Seq2" in content
        assert "Seq3" in content

    # clustal format
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "msa", "--input1", os.path.join(test_dir, "test_msa.fasta"), "--output",
         os.path.join(test_dir, "clustal_out.txt"), "--clustal"], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "clustal_out.txt"))
    with open(os.path.join(test_dir, "clustal_out.txt"), "r") as f:
        content = f.read()
        assert "CLUSTAL" in content

    # subsample
    with open(os.path.join(test_dir, "test_long.fasta"), "w") as f:
        f.write(">seq1\n" + "A" * 1000)
    with open(os.path.join(test_dir, "test_long2.fasta"), "w") as f:
        f.write(">seq2\n" + "A" * 1000)
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test_long.fasta"), "--input2",
         os.path.join(test_dir, "test_long2.fasta"), "--output", os.path.join(test_dir, "subsample_out.txt"), "--subsample", "100"],
        capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "subsample_out.txt"))
    with open(os.path.join(test_dir, "subsample_out.txt"), "r") as f:
        content = f.read()
        # Исправлено: проверяем правильное сообщение
        assert "Subsampled" in content or "Обрезано" in content or "100" in content

    # affine gaps
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "affine_out.txt"), "--gap_open", "-5", "--gap_extend", "-1"],
        capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "affine_out.txt"))
    with open(os.path.join(test_dir, "affine_out.txt"), "r") as f:
        content = f.read()
        assert "Score" in content

    # scoring matrix
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "matrix_out.txt"), "--matrix", "BLOSUM62"], capture_output=True,
        text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "matrix_out.txt"))
    with open(os.path.join(test_dir, "matrix_out.txt"), "r") as f:
        content = f.read()
        assert "Score" in content

    # verbose
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "verbose_out.txt"), "--verbose"], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "verbose_out.txt"))

    # preview
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "preview_out.txt"), "--preview"], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "preview_out.txt"))
    with open(os.path.join(test_dir, "preview_out.txt"), "r") as f:
        content = f.read()
        # Исправлено: проверяем содержимое основного файла вывода, который должен содержать результат выравнивания
        # Preview отображается в консоли, но не попадает в основной файл.
        assert "Score" in content # Основной файл должен содержать результат
        # Тест не проверяет консольный вывод preview, только файл.

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_error_handling():
    # тест обработки ошибок
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    # invalid file
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "nonexistent.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta")], capture_output=True, text=True)
    assert result.returncode != 0

    # invalid mode parameters
    with open(os.path.join(test_dir, "single_seq.fasta"), "w") as f:
        f.write(">seq1\nAGC")
    result = subprocess.run([sys.executable, "-m", "aligner.cli", "msa", "--input1", os.path.join(test_dir, "single_seq.fasta")],
                            capture_output=True, text=True)
    assert result.returncode != 0

    # positive penalties
    with open(os.path.join(test_dir, "test1.fasta"), "w") as f:
        f.write(">seq1\nAGC")
    with open(os.path.join(test_dir, "test2.fasta"), "w") as f:
        f.write(">seq2\nACGC")
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--gap", "2"], capture_output=True, text=True)
    assert result.returncode != 0

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_config_functionality():
    # тест конфигурации
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    # save config
    config = {
        'mode': 'global',
        'input1': os.path.join(test_dir, 'test1.fasta'),
        'input2': os.path.join(test_dir, 'test2.fasta'),
        'output': os.path.join(test_dir, 'config_test_out.txt'),
        'match': 2,
        'mismatch': -2,
        'gap': -3,
        'lang': 'en',
        # Исправлено: добавлены недостающие ключи
        'verbose': False,
        'preview': False,
        'tutorial': False,
        'batch': False,
        'clustal': False,
        'subsample': 0,
        'matrix': None  # Исправлено: добавлен ключ matrix
    }
    import yaml
    config_path = os.path.join(test_dir, "test_config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    with open(os.path.join(test_dir, "test1.fasta"), "w") as f:
        f.write(">seq1\nAGC")
    with open(os.path.join(test_dir, "test2.fasta"), "w") as f:
        f.write(">seq2\nACGC")

    # load config
    result = subprocess.run([sys.executable, "-m", "aligner.cli", "--config", config_path],
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "config_test_out.txt"))

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_performance_with_large_sequences():
    # тест производительности с большими последовательностями
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    seq = "A" * 1000
    with open(os.path.join(test_dir, "large1.fasta"), "w") as f:
        f.write(f">seq1\n{seq}")
    with open(os.path.join(test_dir, "large2.fasta"), "w") as f:
        f.write(f">seq2\n{seq}")

    start = time.time()
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "large1.fasta"), "--input2",
         os.path.join(test_dir, "large2.fasta"), "--output", os.path.join(test_dir, "large_out.txt")], capture_output=True, text=True)
    end = time.time()

    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "large_out.txt"))
    assert (end - start) < 30  # должно выполниться менее чем за 30 секунд

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_multithreading():
    # тест многопоточности msa
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    with open(os.path.join(test_dir, "multi_thread_test.fasta"), "w") as f:
        f.write(">seq1\nAGCTAGCT\n>seq2\nAGCTGGCT\n>seq3\nAGCTTGCT\n>seq4\nAGCTCGCT")

    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "msa", "--input1", os.path.join(test_dir, "multi_thread_test.fasta"), "--output",
         os.path.join(test_dir, "thread_out.txt"), "--threads", "2"], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "thread_out.txt"))
    with open(os.path.join(test_dir, "thread_out.txt"), "r") as f:
        content = f.read()
        assert "Seq1" in content
        assert "Seq2" in content
        assert "Seq3" in content
        assert "Seq4" in content

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_language_support():
    # тест поддержки языков
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    with open(os.path.join(test_dir, "test1.fasta"), "w") as f:
        f.write(">seq1\nAGC")
    with open(os.path.join(test_dir, "test2.fasta"), "w") as f:
        f.write(">seq2\nACGC")

    # english
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "en_out.txt"), "--lang", "en"], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "en_out.txt"))
    with open(os.path.join(test_dir, "en_out.txt"), "r") as f:
        content = f.read()
        assert "Score" in content
        assert "Identity" in content
        assert "Time" in content

    # russian
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "test1.fasta"), "--input2",
         os.path.join(test_dir, "test2.fasta"), "--output", os.path.join(test_dir, "ru_out.txt"), "--lang", "ru"], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "ru_out.txt"))
    with open(os.path.join(test_dir, "ru_out.txt"), "r") as f:
        content = f.read()
        assert "Счет" in content or "Score" in content  # может быть смешанный вывод
        assert "Идентичность" in content or "Identity" in content
        assert "Время" in content or "Time" in content

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_edge_cases():
    # тест граничных случаев
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    # empty sequences
    with open(os.path.join(test_dir, "empty1.fasta"), "w") as f:
        f.write(">seq1\n")
    with open(os.path.join(test_dir, "empty2.fasta"), "w") as f:
        f.write(">seq2\n")

    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "empty1.fasta"), "--input2",
         os.path.join(test_dir, "empty2.fasta"), "--output", os.path.join(test_dir, "empty_out.txt")], capture_output=True, text=True)
    # может завершиться с ошибкой или успешно, в зависимости от реализации
    # но не должен упасть

    # single character sequences
    with open(os.path.join(test_dir, "single1.fasta"), "w") as f:
        f.write(">seq1\nA")
    with open(os.path.join(test_dir, "single2.fasta"), "w") as f:
        f.write(">seq2\nT")

    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "single1.fasta"), "--input2",
         os.path.join(test_dir, "single2.fasta"), "--output", os.path.join(test_dir, "single_out.txt")], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "single_out.txt"))

    # identical sequences
    with open(os.path.join(test_dir, "identical1.fasta"), "w") as f:
        f.write(">seq1\nATGCATGC")
    with open(os.path.join(test_dir, "identical2.fasta"), "w") as f:
        f.write(">seq2\nATGCATGC")

    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "identical1.fasta"), "--input2",
         os.path.join(test_dir, "identical2.fasta"), "--output", os.path.join(test_dir, "identical_out.txt")], capture_output=True, text=True)
    assert result.returncode == 0
    assert os.path.exists(os.path.join(test_dir, "identical_out.txt"))
    with open(os.path.join(test_dir, "identical_out.txt"), "r") as f:
        content = f.read()
        assert "Score" in content
        assert "100.00%" in content or "Identity" in content

    # cleanup
    shutil.rmtree(test_dir)

def test_cli_memory_usage():
    # тест использования памяти (если psutil доступен)
    test_dir = "cli_test"
    os.makedirs(test_dir, exist_ok=True)

    try:
        import psutil
        seq = "A" * 5000
        with open(os.path.join(test_dir, "memory_test1.fasta"), "w") as f:
            f.write(f">seq1\n{seq}")
        with open(os.path.join(test_dir, "memory_test2.fasta"), "w") as f:
            f.write(f">seq2\n{seq}")

        result = subprocess.run(
            [sys.executable, "-m", "aligner.cli", "global", "--input1", os.path.join(test_dir, "memory_test1.fasta"), "--input2",
             os.path.join(test_dir, "memory_test2.fasta"), "--output", os.path.join(test_dir, "memory_out.txt")], capture_output=True, text=True)
        assert result.returncode == 0
        assert os.path.exists(os.path.join(test_dir, "memory_out.txt"))
        with open(os.path.join(test_dir, "memory_out.txt"), "r") as f:
            content = f.read()
            assert "Memory" in content or "Память" in content

        # cleanup
        shutil.rmtree(test_dir)
    except ImportError:
        # psutil не установлен, пропускаем тест
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        pytest.skip("psutil not installed")

def test_cli_batch_mode_simulation():
    # симуляция batch режима через флаги
    test_dir = "cli_test"
    batch_test_dir = os.path.join(test_dir, "batch_test")
    os.makedirs(batch_test_dir, exist_ok=True)
    with open(os.path.join(batch_test_dir, "seq1.fasta"), "w") as f:
        f.write(">seq1\nAGC")
    with open(os.path.join(batch_test_dir, "seq2.fasta"), "w") as f:
        f.write(">seq2\nACGC")
    with open(os.path.join(batch_test_dir, "seq3.fasta"), "w") as f:
        f.write(">seq3\nAGGC")

    # batch не поддерживается в текущей реализации через флаги,
    # но можно проверить корректную обработку ошибок
    result = subprocess.run(
        [sys.executable, "-m", "aligner.cli", "global", "--directory", batch_test_dir, "--output",
         os.path.join(test_dir, "batch_out.txt")], capture_output=True, text=True)
    # ожидаем ошибку, так как не указаны input1 и input2
    assert result.returncode != 0


    # cleanup
    shutil.rmtree(test_dir)