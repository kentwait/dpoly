from pipeline.model import Sequence

class TestSequence:
    def setup(self):
        self.seq = Sequence('test', 'ATGCATGCATGCAAA')

    def teardown(self):
        pass

    def test_fasta_format(self):
        # No line width specified
        fasta_string = self.seq.fasta_format()
        assert fasta_string == '>test\nATGCATGCATGCAAA'

        # Line width is shorter than the sequence.
        # Sequence should wrap according to the line width while the ID
        # line remains unaffected.
        fasta_string = self.seq.fasta_format(line_width=6)
        assert fasta_string == '>test\nATGCAT\nGCATGC\nAAA'

        # Line width is longer than the sequence.
        # No wrapper is expected.
        fasta_string = self.seq.fasta_format(line_width=20)
        assert fasta_string == '>test\nATGCATGCATGCAAA'
