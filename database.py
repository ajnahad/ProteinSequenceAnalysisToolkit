"""
database.py

SequenceDatabase: Collection management
Responsibility: Store, retrieve, search, filter sequences, file I/O

"""
from sequence import ProteinSequence, DNASequence
class SequenceDatabase:
    def __init__(self, name):
        self.name = name
        self.sequences = {}
    def add_sequence(self, sequence, id, description):
        if id in self.sequences.keys():
            raise KeyError("This ID is taken.")
        elif sequence in self.sequences.values():
            raise ValueError("This sequence is already in the database.")
        else:
            if "ATCG" in self.sequence:
                sequence[id] = DNASequence(sequence, id, description)
            else:
                sequence[id] = ProteinSequence(sequence, id, description)
    def remove_sequence(self, id):
        if id not in self.sequence.keys():
            raise KeyError("This ID does not exist in the database. ")
        else:
            del self.sequences[id]
    def get_sequence(self, id):
        if id not in self.sequence.keys():
            raise KeyError("This ID does not exist in the database. ")
        else:
            return self.sequences[id]
    def search_by_keyword(self, keyword):
        for seq in self.sequences.values():
            if keyword.lower() in seq.description:
                return seq.id
        print(f"No sequence descriptions containing {keyword} found. ")
    def filter_by_length(self, min_length, max_length):
        return {k: v for k, v in self.sequences.items() if min_length <= len(v) <= max_length}
    def get_statistics(self):
        lengths = [len(s.sequence) for s in self.sequences.values()]
        if not lengths:
            return {"count" : 0, "avg_length" : 0, "min_length" : 0, "max_length" : 0}
        return{
            "count" : len(lengths),
            "avg_length" : round(sum(lengths) / len(lengths), 2),
            "min_length" : min(lengths),
            "max_length" : max(lengths)
        }
    def __len__(self):
        return len(self.sequences)
    def __iter__(self):
        return iter(self.sequences.values())