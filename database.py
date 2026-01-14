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
    def add_sequence(self, sequence, s_id, description):
        if s_id in self.sequences.keys():
            raise KeyError("This ID is taken.")
        elif any(sequence == seq.sequence for seq in self.sequences.values()):
            raise ValueError("This sequence is already in the database.")
        else:
            if set(sequence.upper()) <= set("ATCG"):
                obj = DNASequence(sequence, s_id, description)
            else:
                obj = ProteinSequence(sequence, s_id, description)
        self.sequences[s_id] = obj    
    def remove_sequence(self, s_id):
        if s_id not in self.sequences.keys():
            raise KeyError("This ID does not exist in the database. ")
        else:
            del self.sequences[s_id]
    def get_sequence(self, s_id):
        if s_id not in self.sequences.keys():
            raise KeyError("This ID does not exist in the database. ")
        else:
            return self.sequences[s_id]
    def search_by_keyword(self, keyword):
        for seq in self.sequences.values():
            if keyword.lower() in seq.description:
                return seq.seq_id
        print(f"No sequence descriptions containing {keyword} found. ")
    def filter_by_length(self, min_length, max_length):
        return {k: v for k, v in self.sequences.items() if min_length <= len(v.sequence) <= max_length}
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
    @staticmethod
    def validate_fasta(filepath, seqtype):
        if seqtype == "nucleotides":
            valid_set = set("ATCGX")
        elif seqtype == "aminoacids":
            valid_set = set("ACDEFGHIKLMNPQRSTVWY*")
        else:
            raise ValueError("seqtype must be 'nucleotides' or 'aminoacids'")
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                if not set(line.upper()) <= valid_set:
                    raise ValueError(f"Invalid character in line: {line}")
    @staticmethod
    def export_to_fasta(sequences, filepath):
        with open(filepath, "w") as f:
            for seq in sequences:
                header = ">" + seq.seq_id + " " + seq.description
                f.write(header + "\n") 
                f.write(seq.sequence + "\n")
    @classmethod
    def merge_databases(cls, d1, d2, database_name):
        db = cls(database_name)
        for seq in d1:
            db.add_sequence(seq.sequence, seq.s_id, seq.description)
        for seq in d2:
            db.add_sequence(seq.sequence, seq.s_id, seq.description)
        return db
    @classmethod
    def from_fasta_file(cls, filepath, name, seqtype):
        cls.validate_fasta(filepath, seqtype)
        db = cls(name)
        seq_id, descr, seq = None, "", []
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_id is not None:
                        if seqtype == "nucleotides":
                           obj = DNASequence("".join(seq), seq_id, descr)
                        elif seqtype == "aminoacids":
                            obj = ProteinSequence("".join(seq), seq_id, descr)
                        db.add_sequence(obj.sequence, obj.s_id, obj.description)
                    header = line[1:]
                    parts = header.split(maxsplit=1)
                    seq_id = parts[0]
                    descr = parts[1] if len(parts) > 1 else ""
                    seq.clear()
                else:
                    seq.append(line)
            if seq_id is not None:
                if seqtype == "nucleotides":
                    obj = DNASequence("".join(seq), seq_id, descr)
                elif seqtype == "aminoacids":
                    obj = ProteinSequence("".join(seq), seq_id, descr)
                db.add_sequence(obj.sequence, obj.s_id, obj.description)
        return db                        
