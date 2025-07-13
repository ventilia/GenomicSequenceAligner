import pytest
from subprocess import run, CalledProcessError
import os
import sys



@pytest.fixture
def dummy_fasta():
    with open("test1.fasta", "w") as f:
        f.write(">seq1\nAGC")
    with open("test2.fasta", "w") as f:
        f.write(">seq2\nACGC")
    with open("test_multi.fasta", "w") as f:
        f.write(">seq1\nAGC\n>seq2\nACGC\n>seq3\nAGGC")
    yield
    os.remove("test1.fasta")
    os.remove("test2.fasta")
    os.remove("test_multi.fasta")

def test_cli_pairwise(dummy_fasta):
    try:
        run([sys.executable, "-m", "aligner.cli", "--input1", "test1.fasta", "--input2", "test2.fasta", "--mode", "global", "--output", "test_out.txt"], check=True)
        assert os.path.exists("test_out.txt")
        with open("test_out.txt", "r") as f:
            content = f.read()
            assert "Alignment Score" in content  # Check output
        os.remove("test_out.txt")
    except CalledProcessError as e:
        pytest.fail(f"CLI failed: {e}")

def test_cli_msa(dummy_fasta):
    try:
        run([sys.executable, "-m", "aligner.cli", "--input1", "test_multi.fasta", "--mode", "msa", "--output", "test_msa.txt"], check=True)
        assert os.path.exists("test_msa.txt")
        with open("test_msa.txt", "r") as f:
            content = f.read()
            assert "Seq1" in content
        os.remove("test_msa.txt")
    except CalledProcessError as e:
        pytest.fail(f"CLI failed: {e}")