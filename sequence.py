"""
sequence.py

Sequence base class: Abstract representation of biological sequences
ProteinSequence: Protein-specific calculations and validation
DNASequence: DNA operations (complement, translation)
Responsibility: Sequence representation, validation, basic operations

"""
import datetime 
import abc
from collections import Counter
# import database as ds

class Sequence(abc.ABC):
    sequence = ""
    id = ""
    description = ""
    creation_date = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    def __init__(self, sequence, id, description):
        #TODO: add validation for sequence and id AFTER writing database.py
        self.sequence = sequence
        self.id = id 
        self.description = description
    def __str__(self):
        return str(self.sequence)
    def __repr__(self):
        return f"Sequence {self.id}: {self.sequence}"
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self, index):
        return self.sequence[index]
    @abc.abstractmethod
    def validate(self):
        pass
    @abc.abstractmethod
    def get_composition(self):
        pass
    def to_fasta(self):
        seq = self.sequence.replace(" ", "")
        return f">{self.description}{seq}"
    def find_motif(self, pattern):
        found_at = []
        k = len(pattern)
        for i in range(0, len(self.sequence)-k):
            if self.sequence[i:i+k] == pattern:
                found_at.append(i)
        return found_at
    @classmethod
    def from_fasta(cls, fasta_string):
        header, _, seq = fasta_string.partition("\n") # splits into id with description - header and sequence 
        seq.replace(" ", "") 
        parts = header[1:].split(maxsplit=1) # removes ">" and split
        seq_id = parts[0] 
        d = parts[1] if len(parts) > 1 else "" # description 
        return cls(seq, seq_id, d)
    def get_subsequence(self, start, end):
        subsequence = self.sequence[start:end]
        id = f"{self.sequence.id}_{start}-{end}"
        description = f"Subsequence of sequence [id = {self.sequence.id}] from position {start} to {end}"
        return self.__class__(subsequence, id, description)
    def get_sequence(self):
        return self.sequence
    def set_sequence(self, seq):
        self.sequence = seq #TODO: add validation like in the constructor

class ProteinSequence(Sequence):
    AMINO_ACIDS = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}
    AMINO_ACID_WEIGHTS = {
        "A": 89.09,  "R": 174.20, "N": 132.12, "D": 133.10,
        "C": 121.16, "E": 147.13, "Q": 146.14, "G": 75.07,
        "H": 155.15, "I": 131.17, "L": 131.17, "K": 146.19,
        "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
        "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15
    }
    # Kyte - Doolittle Hydropathy scale
    HYDROPATHY_INDEX = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5,
        "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
        "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
        "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
        "T": -0.7,"W": -0.9, "Y": -1.3,"V": 4.2
    } 
    def validate(self):
        if not all(aa in self.AMINO_ACIDS for aa in self.sequence):
            raise ValueError("Sequence contains invalid aminoacids. ")
        return True
    def calculate_molecular_weight(self):
        total = 0
        for value in self.sequence:
            total += self.AMINO_ACID_WEIGHTS[value]
        return total - (len(self.sequence) - 1) * 18.015
    def calculate_hydropathy(self):
        total = 0
        for value in self.sequence:
            total += self.HYDROPATHY_INDEX[value]
        return total / len(self.sequence)
    def get_composition(self):
        counts = Counter(self.sequence)
        total = len(self.sequence)
        return {aa: round(count/total, 2) for aa, count in counts.items()}