"""
alignment.py

AlignmentScorer: Pairwise alignment
Responsibility: Load scoring matrices, perform Needleman-Wunsch, score alignments

"""
from Bio.Align import substitution_matrices, PairwiseAligner
class AlignmentScorer:
    matrix = substitution_matrices.load("BLOSUM62")
    def __init__(self, gap_penalty=-10, gap_extension=-1):
        """
        Initializes the alignment scorer
        -parameters:
            gap_penalty [int]: penalty for opening a gap
            gap_extension [int]: penalty for extending an existing gap
        """
        self.gap_penalty = gap_penalty
        self.gap_extension = gap_extension
    def align(self, seq1, seq2, algorithm): # can be "global" or "local"
        """
        Aligns two sequences using global or local alignment
        -parameters:
            seq1 [str]: first sequence
            seq2 [str]: second sequence
            algorithm [str]: "global" or "local"
        -returns a tuple: (aligned_seq1, aligned_seq2, alignment_score)
        """
        aligner = PairwiseAligner()
        if algorithm not in {"global", "local"}:
            raise ValueError("algorithm must be 'global' or 'local'")
        aligner.mode = "global" if algorithm == "global" else "local" 
        aligner.substitution_matrix = self.matrix
        aligner.open_gap_score = self.gap_penalty
        aligner.extend_gap_score = self.gap_extension
        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0] # takes the best scoring alignment
        aligned_seq1 = str(best_alignment[0])
        aligned_seq2 = str(best_alignment[1])
        score = best_alignment.score
        return aligned_seq1, aligned_seq2, score
    def score_alignment(self, seq1, seq2, algorithm):
        """
        Manually computes the alignment score using gap penalties
        and the substitution matrix
        -parameters:
            seq1 [str]: first sequence
            seq2 [str]: second sequence
            algorithm [str]: "global" or "local"
        -returns a float - total alignment score
        """
        a_seq1, a_seq2, _ = self.align(seq1, seq2, algorithm)
        total_score = 0
        gap_in_seq1 = False
        gap_in_seq2 = False
        for pos in range(len(a_seq1)):
            char1 = a_seq1[pos]
            char2 = a_seq2[pos]
            if char1 != "-" and char2 != "-":
                total_score += self.matrix[char1, char2]
                gap_in_seq1 = False
                gap_in_seq2 = False
            elif char1 == "-" and char2 != "-":
                if not gap_in_seq1:
                    total_score += self.gap_penalty
                    gap_in_seq1 = True
                else:
                    total_score += self.gap_extension
            elif char1 != "-" and char2 == "-":
                if not gap_in_seq2:
                    total_score+=self.gap_penalty
                    gap_in_seq2 = True
                else:
                    total_score += self.gap_extension
        return total_score
    @staticmethod
    def calculate_identity(a_seq1, a_seq2):
        """
        Calculates percent identity between two aligned sequences
        -parameters:
            a_seq1 [str]: aligned sequence 1
            a_seq2 [str]: aligned sequence 2
        -returns a float - percent identity (0-100)
        """
        matches, valid_positions = 0, 0
        for pos in range(len(a_seq1)):
            if a_seq1[pos] == "-" or a_seq2[pos] == "-":
                continue
            else:
                valid_positions+=1
                if a_seq1[pos] == a_seq2[pos]:
                    matches += 1
        if valid_positions == 0:
            return 0.0
        return (matches/valid_positions) * 100
    @staticmethod
    def calculate_similarity(a_seq1, a_seq2, matrix): # requires a preloaded matrix object
        """
        Calculates percent similarity using a substitution matrix
        A position is considered similar if the matrix score > 0
        -parameters:
            a_seq1 [str]: aligned sequence 1
            a_seq2 [str]: aligned sequence 2
            matrix: substitution matrix object
        -returns a float: percent similarity (0-100)
        """
        similar, valid_positions = 0, 0
        for pos in range(len(a_seq1)):
            if a_seq1[pos] == "-" or a_seq2[pos] == "-":
                continue
            else:
                valid_positions+=1
                if matrix[a_seq1[pos], a_seq2[pos]] > 0:
                    similar += 1
        if valid_positions == 0:
            return 0.0
        return (similar/valid_positions) * 100
