"""
main.py

Main workflow
Toolkit usage demonstration

"""
from database import SequenceDatabase
from alignment import AlignmentScorer
from analysis import SequenceAnalyzer
from comparison import SequenceComparator
from reporting import AnalysisReport
def main():
    db = SequenceDatabase("Example Database")
    db.add_sequence("MKTFFVAGVILLLSAAPVFA", "prot1", "Signal peptide protein")
    db.add_sequence("MKTFFVAGVILLLSAAPIFA", "prot2", "Similar signal peptide")
    seq1 = db.get_sequence("prot1")
    seq2 = db.get_sequence("prot2")
    aligner = AlignmentScorer()
    aligned1, aligned2, score = aligner.align(seq1.sequence, seq2.sequence, "global")
    identity = aligner.calculate_identity(aligned1, aligned2)
    analyzer = SequenceAnalyzer(seq1.sequence)
    props = analyzer.basic_properties
    pI = analyzer.calculate_isoelectric_point()
    comparator = SequenceComparator([aligned1, aligned2])
    distances = comparator.calculate_pairwise_distances()
    report = AnalysisReport("Protein Analysis Report")
    report.add_section("Alignment Score", str(score))
    report.add_section("Percent Identity", f"{identity:.2f}%")
    report.add_section("Basic Properties", str(props))
    report.add_section("Isoelectric Point", str(pI))
    report.add_section("Distances", str(distances))
    report.save_to_file("analysis_report.txt")
if __name__ == "__main__":
    main()
