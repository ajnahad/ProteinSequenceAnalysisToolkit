"""
comparison.py

SequenceComparator: Multi-sequence comparison
Responsibility: Find conserved residues, calculate distances, consensus sequence

"""
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
import Levenshtein

class SequenceComparator:
    def __init__(self, sequences):
        """
        sequences: list of aligned sequence strings
        """
        self.sequences = sequences
        self.comparison_results = {}
        self.msa = MultipleSeqAlignment(
            [SeqRecord(Seq(seq), id=f"seq{i}") for i, seq in enumerate(self.sequences)]
        )
    def find_conserved_residues(self, threshold=0.9):
        """
        Returns a consensus string of highly conserved residues.
        """
        summary = AlignInfo.SummaryInfo(self.msa)
        consensus = summary.gap_consensus(threshold=threshold)
        return consensus
    def calculate_pairwise_distances(self):
        """
        Calculates a pairwise distance matrix using BLOSUM62.
        Stores results in self.comparison_results.
        """
        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(self.msa)
        names = [f"seq{i}" for i in range(len(self.sequences))]
        for i, name1 in enumerate(names):
            self.comparison_results[name1] = {}
            for j, name2 in enumerate(names):
                if i == j:
                    self.comparison_results[name1][name2] = 0.0
                else:
                    self.comparison_results[name1][name2] = dm[i, j]
        return self.comparison_results
    def group_by_similarity(self, threshold):
        """
        Groups sequences whose similarity (1 - distance) is >= threshold.
        Returns list of lists of sequences.
        """
        if not self.comparison_results:
            self.calculate_pairwise_distances()
        groups = []
        unassigned = set(range(len(self.sequences)))
        while unassigned:
            i = unassigned.pop()
            group = {i}
            for j in list(unassigned):
                similarity = 1 - self.comparison_results[f"seq{i}"][f"seq{j}"]
                if similarity >= threshold:
                    group.add(j)
                    unassigned.remove(j)
            groups.append([self.sequences[idx] for idx in group])
        return groups
    def get_unique_sequences(self):
        return set(self.sequences)
    @staticmethod
    def hamming_distance(seq1, seq2):
        """
        Counts number of positions where sequences differ.
        Only works if sequences are the same length.
        """
        return sum(a != b for a, b in zip(seq1, seq2))
    @staticmethod
    def levenshtein_distance(seq1, seq2):
        """
        Returns edit distance between sequences.
        """
        return Levenshtein.distance(seq1, seq2)
    @classmethod
    def from_database(cls, database, sequence_ids):
        """
        Builds a SequenceComparator from a SequenceDatabase object.
        Extracts sequence strings from database objects.
        """
        sequences = []
        for seq_id in sequence_ids:
            if seq_id in database.sequences:
                seq_obj = database.sequences[seq_id]
                sequences.append(str(seq_obj.sequence))
        return cls(sequences)
