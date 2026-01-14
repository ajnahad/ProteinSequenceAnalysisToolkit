"""
analysis.py

SequenceAnalyzer: Biochemical property calculations
Responsibility: Molecular weight, pI, hydropathy, charge at pH, instability index

"""
from collections import Counter 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
class SequenceAnalyzer:
    pKa_N_term = 9.69
    pKa_C_term = 2.34
    pKa_sidechains = {
        "K": 10.53,  # Lys
        "R": 12.48,  # Arg
        "H": 6.00,   # His
        "D": 3.86,   # Asp
        "E": 4.25,   # Glu
        "C": 8.33,   # Cys
        "Y": 10.07   # Tyr
    }
    CHOU_FASMAN_PROPENSITY = {
        "A": {"H": 1.42, "E": 0.83, "C": 0.66},
        "R": {"H": 0.98, "E": 0.93, "C": 0.95},
        "N": {"H": 0.67, "E": 0.89, "C": 1.56},
        "D": {"H": 1.01, "E": 0.54, "C": 1.46},
        "C": {"H": 0.70, "E": 1.19, "C": 1.19},
        "Q": {"H": 1.11, "E": 1.10, "C": 0.98},
        "E": {"H": 1.51, "E": 0.37, "C": 0.74},
        "G": {"H": 0.57, "E": 0.75, "C": 1.64},
        "H": {"H": 1.00, "E": 0.87, "C": 1.01},
        "I": {"H": 1.08, "E": 1.60, "C": 0.47},
        "L": {"H": 1.21, "E": 1.30, "C": 0.59},
        "K": {"H": 1.16, "E": 0.74, "C": 1.01},
        "M": {"H": 1.45, "E": 1.05, "C": 0.60},
        "F": {"H": 1.13, "E": 1.38, "C": 0.59},
        "P": {"H": 0.57, "E": 0.55, "C": 1.52},  # helix breaker
        "S": {"H": 0.77, "E": 0.75, "C": 1.43},
        "T": {"H": 0.83, "E": 1.19, "C": 1.19},
        "W": {"H": 1.08, "E": 1.37, "C": 0.62},
        "Y": {"H": 0.69, "E": 1.47, "C": 1.14},
        "V": {"H": 1.06, "E": 1.70, "C": 0.50}
    }
    DOMAIN_PATTERNS = {
        # Zinc Finger motif (Cys-Cys-His-His) simplified
        "Zinc_Finger": "CXXCXXXXXXXXXXXXHXXXH",
        # SH3 domain (Pro-rich motif)
        "SH3": "PXXPXXP",
        # ATP binding motif (Walker A motif, simplified)
        "Walker_A": "GXXXXGK",
        # EF-hand calcium binding motif
        "EF_hand": "DXDXXD",
        # Leucine zipper (every 7th residue is L)
        "Leucine_zipper": "LXXXXXXLXXXXXXL",
        # Helix-turn-helix motif
        "HTH": "XNXAXXKX",
    }
    # Kyte - Doolittle Hydropathy scale
    HYDROPATHY_INDEX = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5,
        "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
        "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
        "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
        "T": -0.7,"W": -0.9, "Y": -1.3,"V": 4.2
    } 
    AMINO_ACID_WEIGHTS = {
        "A": 89.09,  "R": 174.20, "N": 132.12, "D": 133.10,
        "C": 121.16, "E": 147.13, "Q": 146.14, "G": 75.07,
        "H": 155.15, "I": 131.17, "L": 131.17, "K": 146.19,
        "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
        "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15
    }
    def __init__(self, sequence):
        """
        Parameter - sequence [str]: protein sequence
        """
        self.sequence = sequence
        self.counts = Counter(sequence)
        self.ionizable_counts = Counter(
            aa for aa in sequence if aa in self.pKa_sidechains
        )
        self._protparam = ProteinAnalysis(sequence)
    def net_charge_at_ph(self, pH):
        """
        Computes net charge at a given pH using Henderson-Hasselbalch
        """
        charge = 0.0
        # Positive groups (N-terminus, K, R, H)
        charge += +1 / (1+10**(pH-self.pKa_N_term)) # N-terminus
        charge += self.ionizable_counts["K"] * (1 / (1 + 10**(pH-self.pKa_sidechains["K"])))
        charge += self.ionizable_counts["R"] * (1 / (1 + 10**(pH-self.pKa_sidechains["R"])))
        charge += self.ionizable_counts["H"] * (1 / (1 + 10**(pH-self.pKa_sidechains["H"])))
        # Negative groups (C-terminus, D, E, C, Y)
        charge += -1 / (1 + 10**(self.pKa_C_term-pH)) # C terminus
        charge += self.ionizable_counts["D"] * (-1 / (1 + 10**(self.pKa_sidechains["D"] - pH)))
        charge += self.ionizable_counts["E"] * (-1 / (1 + 10**(self.pKa_sidechains["E"] - pH)))
        charge += self.ionizable_counts["C"] * (-1 / (1 + 10**(self.pKa_sidechains["C"] - pH)))
        charge += self.ionizable_counts["Y"] * (-1 / (1 + 10**(self.pKa_sidechains["Y"] - pH)))
        return charge
    def calculate_isoelectric_point(self):
        """
        Estimates the isoelectric point by brute-force pH scanning
        """
        pH = 0.0
        step = 0.01
        while pH <= 14.0:
            if abs(self.net_charge_at_ph(pH)) < 0.01: # tolerance
                return round(pH, 2)
            pH += step
        return None
    def predict_secondary_structure(self):
        """
        Predicts secondary structure using Chou-Fasman propensities
        -returns a str - a sequence of "H" (helix), "S" (sheet), and "C" (coil)
        """
        secondary_sequence = ""
        for residue in self.sequence:
            helix_score = self.CHOU_FASMAN_PROPENSITY[residue]["H"]
            sheet_score = self.CHOU_FASMAN_PROPENSITY[residue]["E"]
            coil_score = self.CHOU_FASMAN_PROPENSITY[residue]["C"]
            max_score = max(helix_score, sheet_score, coil_score)
            if max_score == helix_score:
                secondary_sequence += "H"
            elif max_score == sheet_score:
                secondary_sequence += "S"
            else:
                secondary_sequence += "C"
        return secondary_sequence
    def identify_domain_patterns(self): 
        """
        Identifies simplified domain motifs using pattern matching
        -returns a dict - domain names mapped to lists of (start, end) positions
        """
        matches = {}
        for name, pattern in self.DOMAIN_PATTERNS.items():
            positions = []
            for i in range(len(self.sequence)-len(pattern)+1): # defining the window size
                match = True  
                for j in range(len(pattern)): # look at each individual aa within the window
                    if pattern[j] != "X" and pattern[j] != self.sequence[i+j]: # checking if pattern matches with X or the exact aminoacid that needs to be at the exact position at the protein
                        match = False
                        break
                if match:
                    positions.append((i, i+len(pattern)))
            if positions:
                matches[name] = positions
        return matches
    def calculate_instability_index(self): # Guruprasad scoring method 
        """
        Calculates protein instability index (Guruprasad method)
        """
        return self._protparam.instability_index()
    def calculate_gravy(self):
        """
        Computes GRAVY (average hydropathy)
        """
        if not self.sequence:
            return 0.0
        total = sum(self.HYDROPATHY_INDEX.get(aa, 0) for aa in self.sequence)
        return total / len(self.sequence)
    @property
    def basic_properties(self):
        """
        Returns basic sequence statistics.
        """
        mw = sum(self.AMINO_ACID_WEIGHTS.get(aa, 0) for aa in self.sequence)
        return {
            "length": len(self.sequence),
            "molecular_weight" : round(mw, 2),
            "GRAVY" : round(self.calculate_gravy(), 2)
        }
