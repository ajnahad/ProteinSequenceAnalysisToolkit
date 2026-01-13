"""
reporting.py

AnalysisReport: Format and output results
Responsibility: Generate formatted text reports, save to files

"""
import datetime
class AnalysisReport:
    def __init__(self, title, format_style="text"):
        self.title = title
        self.format_style = format_style
        self.sections = {}
        self.timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    def add_section(self, section_name, content):
        """ 
        Adds text content to a section.
        """
        self.sections[section_name] = content
    def add_table(self, section_name, data, headers):
        """
        Adds a table to a section.
        data: list of lists (rows)
        headers: list of column headers
        """
        table_str = self.format_table(data, headers, self.format_style)
        self.sections[section_name] = table_str
    def generate(self):
        lines = [self.title, f"Generated: {self.timestamp}", ""]
        for section, content in self.sections.items():
            lines.append(section)
            lines.append("-" * len(section))  
            lines.append(content)
            lines.append("")  
        return "\n".join(lines)
    def save_to_file(self, filepath):
        with open(filepath, "w") as f:
            f.write(self.generate())
    @staticmethod
    def format_table(data, headers, style="text"):
        lines = ["\t".join(headers)]
        for row in data:
            lines.append("\t".join(str(item) for item in row))
        return "\n".join(lines)
    @staticmethod
    def create_summary_statistics(sequences):
        lengths = [len(seq) for seq in sequences]
        return {
            "count": len(sequences),
            "avg_length": round(sum(lengths)/len(lengths), 2) if lengths else 0,
            "min_length": min(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0
        }
    @classmethod
    def from_sequence_analysis(cls, sequences, title):
        report = cls(title)
        stats = cls.create_summary_statistics(sequences)
        report.add_section("Summary Statistics", str(stats))
        return report
    @classmethod
    def from_alignment(cls, alignment_result, title):
        report = cls(title)
        report.add_section("Alignment", str(alignment_result))
        return report
    