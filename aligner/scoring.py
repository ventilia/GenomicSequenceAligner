from Bio.Align import substitution_matrices
from typing import Dict, Tuple, Optional

def load_scoring_matrix(name: str = "BLOSUM62") -> Dict[Tuple[str, str], int]:

    #Загружает scoring matrix из биопит


    try:
        matrix = substitution_matrices.load(name)
        scoring_dict = {}
        for i, row in enumerate(matrix):
            for j, score in enumerate(row):
                key = (matrix.alphabet[i], matrix.alphabet[j])
                scoring_dict[key] = int(score)
        return scoring_dict
    except ValueError:
        raise ValueError(f"маттрица '{name}' не найдена в biopython.")