"""
sequence.py

Sequence base class: Abstract representation of biological sequences
ProteinSequence: Protein-specific calculations and validation
DNASequence: DNA operations (complement, translation)
Responsibility: Sequence representation, validation, basic operations

"""
import datetime 
import abc
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
    @classmethod
    def from_fasta(cls, fasta_string):
        header, _, seq = fasta_string.partition("\n") # splits into id with description - header and sequence 
        seq.replace(" ", "") 
        parts = header[1:].split(maxsplit=1) # removes ">" and split
        seq_id = parts[0] 
        d = parts[1] if len(parts) > 1 else "" # description 
        return cls(seq, seq_id, d)
    def get_sequence(self):
        return self.sequence
    def set_sequence(self, seq):
        self.sequence = seq #TODO: add validation like in the constructor

    