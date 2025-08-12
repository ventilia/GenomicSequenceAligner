import pytest
import os
import sys
import pexpect
from subprocess import run, CalledProcessError

@pytest.fixture
def dummy_fasta():
    # создаем тестовую директорию и файлы
    os.makedirs("test_dir", exist_ok=True)
    with open("test_dir/test1.fasta", "w") as f:
        f.write(">seq1\nAGC")
    with open("test_dir/test2.fasta", "w") as f:
        f.write(">seq2\nACGC")
    with open("test_dir/test_multi.fasta", "w") as f:
        f.write(">seq1\nAGC\n>seq2\nACGC\n>seq3\nAGGC")
    with open("test_dir/config.yaml", "w") as f:
        f.write("mode: global\ninput1: test_dir/test1.fasta\ninput2: test_dir/test2.fasta\noutput: test_out.yaml\nmatch: 1\nmismatch: -1\ngap: -2\nlang: en")
    yield
    for f in ["test_dir/test1.fasta", "test_dir/test2.fasta", "test_dir/test_multi.fasta", "test_dir/config.yaml"]:
        if os.path.exists(f):
            os.remove(f)
    if os.path.exists("test_dir"):
        os.rmdir("test_dir")
    if os.path.exists("alignment.txt"):
        os.remove("alignment.txt")
    if os.path.exists("config.yaml"):
        os.remove("config.yaml")

def test_cli_wizard_global(dummy_fasta):
    # тест интерактивного режима для global
    child = pexpect.spawn(sys.executable, ["-m", "aligner.cli"])
    child.expect("Choose language")
    child.sendline("en")
    child.expect("Choose mode")
    child.sendline("global")
    child.expect("Run batch mode")
    child.sendline("n")
    child.expect("Select directory")
    child.sendline("test_dir")
    child.expect("Select FASTA file for input1")
    child.sendline("test1.fasta")
    child.expect("Select FASTA file for input2")
    child.sendline("test2.fasta")
    for _ in range(8):  # параметры
        child.sendline("")
    child.expect("Verbose logging")
    child.sendline("n")
    child.expect("Preview first 100 bp")
    child.sendline("y")
    child.expect("Run tutorial")
    child.sendline("n")
    child.expect("Output file")
    child.sendline("")
    child.expect("Load config")
    child.sendline("n")
    child.expect("Save current params")
    child.sendline("n")
    child.expect("Run with these params")
    child.sendline("y")
    child.expect("Success")
    child.close()
    assert os.path.exists("alignment.txt")
    with open("alignment.txt", "r") as f:
        content = f.read()
        assert "Score" in content
        assert "Identity" in content
    os.remove("alignment.txt")

def test_cli_wizard_msa(dummy_fasta):
    # тест для msa
    child = pexpect.spawn(sys.executable, ["-m", "aligner.cli"])
    child.expect("Choose language")
    child.sendline("ru")
    child.expect("Выберите режим")
    child.sendline("msa")
    child.expect("Run batch mode")
    child.sendline("n")
    child.expect("Select directory")
    child.sendline("test_dir")
    child.expect("Select FASTA file for input1")
    child.sendline("test_multi.fasta")
    for _ in range(8):
        child.sendline("")
    child.expect("Threads for MSA")
    child.sendline("")
    child.expect("Clustal format for MSA")
    child.sendline("n")
    child.expect("Verbose logging")
    child.sendline("n")
    child.expect("Preview first 100 bp")
    child.sendline("n")
    child.expect("Run tutorial")
    child.sendline("n")
    child.expect("Output file")
    child.sendline("")
    child.expect("Load config")
    child.sendline("n")
    child.expect("Save current params")
    child.sendline("n")
    child.expect("Run with these params")
    child.sendline("y")
    child.expect("Успех")
    child.close()
    assert os.path.exists("alignment.txt")
    with open("alignment.txt", "r") as f:
        content = f.read()
        assert "Seq1" in content
    os.remove("alignment.txt")

def test_cli_batch_mode(dummy_fasta):
    # тест batch
    child = pexpect.spawn(sys.executable, ["-m", "aligner.cli"])
    child.expect("Choose language")
    child.sendline("en")
    child.expect("Choose mode")
    child.sendline("global")
    child.expect("Run batch mode")
    child.sendline("y")
    child.expect("Select directory")
    child.sendline("test_dir")
    for _ in range(8):
        child.sendline("")
    child.expect("Verbose logging")
    child.sendline("n")
    child.expect("Preview first 100 bp")
    child.sendline("n")
    child.expect("Run tutorial")
    child.sendline("n")
    child.expect("Output file")
    child.sendline("")
    child.expect("Load config")
    child.sendline("n")
    child.expect("Save current params")
    child.sendline("n")
    child.expect("Run with these params")
    child.sendline("y")
    child.expect("Success")
    child.close()
    assert os.path.exists("alignment.txt")
    with open("alignment.txt", "r") as f:
        content = f.read()
        assert "test1.fasta vs test2.fasta" in content
    os.remove("alignment.txt")

def test_cli_config_load(dummy_fasta):
    # тест загрузки config
    child = pexpect.spawn(sys.executable, ["-m", "aligner.cli"])
    child.expect("Choose language")
    child.sendline("en")
    child.expect("Choose mode")
    child.sendline("global")
    child.expect("Run batch mode")
    child.sendline("n")
    child.expect("Select directory")
    child.sendline("test_dir")
    child.expect("Select FASTA file for input1")
    child.sendline("test1.fasta")
    child.expect("Select FASTA file for input2")
    child.sendline("test2.fasta")
    for _ in range(8):
        child.sendline("")
    child.expect("Verbose logging")
    child.sendline("n")
    child.expect("Preview first 100 bp")
    child.sendline("n")
    child.expect("Run tutorial")
    child.sendline("n")
    child.expect("Output file")
    child.sendline("")
    child.expect("Load config")
    child.sendline("test_dir/config.yaml")
    child.expect("Save current params")
    child.sendline("n")
    child.expect("Run with these params")
    child.sendline("y")
    child.expect("Success")
    child.close()
    assert os.path.exists("test_out.yaml")
    with open("test_out.yaml", "r") as f:
        content = f.read()
        assert "Score" in content
    os.remove("test_out.yaml")

def test_cli_pairwise_flags(dummy_fasta):
    # тест флагов для pairwise
    try:
        run([sys.executable, "-m", "aligner.cli", "global", "--input1", "test_dir/test1.fasta", "--input2", "test_dir/test2.fasta", "--output", "test_out.txt", "--lang", "ru"], check=True)
        assert os.path.exists("test_out.txt")
        with open("test_out.txt", "r") as f:
            content = f.read()
            assert "Score" in content or "Счет" in content
            assert "Identity" in content or "Идентичность" in content
        os.remove("test_out.txt")
    except CalledProcessError as e:
        pytest.fail(f"CLI failed: {e}")

def test_cli_msa_flags(dummy_fasta):
    # тест флагов для msa
    try:
        run([sys.executable, "-m", "aligner.cli", "msa", "--input1", "test_dir/test_multi.fasta", "--output", "test_msa.txt", "--lang", "en"], check=True)
        assert os.path.exists("test_msa.txt")
        with open("test_msa.txt", "r") as f:
            content = f.read()
            assert "Seq1" in content
        os.remove("test_msa.txt")
    except CalledProcessError as e:
        pytest.fail(f"CLI failed: {e}")

def test_cli_error_handling(dummy_fasta):
    # тест ошибок
    child = pexpect.spawn(sys.executable, ["-m", "aligner.cli"])
    child.expect("Choose language")
    child.sendline("en")
    child.expect("Choose mode")
    child.sendline("msa")
    child.expect("Run batch mode")
    child.sendline("n")
    child.expect("Select directory")
    child.sendline("nonexistent")
    child.expect("No FASTA files found")
    child.expect("Change directory")
    child.sendline("n")
    child.close()

def test_cli_preview(dummy_fasta):
    # тест preview
    child = pexpect.spawn(sys.executable, ["-m", "aligner.cli"])
    child.expect("Choose language")
    child.sendline("en")
    child.expect("Choose mode")
    child.sendline("global")
    child.expect("Run batch mode")
    child.sendline("n")
    child.expect("Select directory")
    child.sendline("test_dir")
    child.expect("Select FASTA file for input1")
    child.sendline("test1.fasta")
    child.expect("Select FASTA file for input2")
    child.sendline("test2.fasta")
    for _ in range(8):
        child.sendline("")
    child.expect("Verbose logging")
    child.sendline("n")
    child.expect("Preview first 100 bp")
    child.sendline("y")
    child.expect("Seq1: AGC")
    child.expect("Run tutorial")
    child.sendline("n")
    child.expect("Output file")
    child.sendline("")
    child.expect("Load config")
    child.sendline("n")
    child.expect("Save current params")
    child.sendline("n")
    child.expect("Run with these params")
    child.sendline("y")
    child.expect("Success")
    child.close()
    assert os.path.exists("alignment.txt")
    os.remove("alignment.txt")