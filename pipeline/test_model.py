from pipeline.model import Sequence, NuclSequence, CodonSequence
from pipeline.model import Alignment


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

        self.seq.description = 'test description'
        fasta_string = self.seq.fasta_format()
        assert fasta_string == '>test test description\nATGCATGCATGCAAA'

    def test_len(self):
        assert len(self.seq) == 15

    def test_getitem(self):
        assert self.seq[0] == 'A'
        assert self.seq[10] == 'G'
        assert self.seq[14] == 'A'
        assert self.seq[-1] == 'A'

        assert self.seq[0:1] == 'A'
        assert self.seq[0:3] == 'ATG'
        assert self.seq[12:] == 'AAA'
        assert self.seq[-3:] == 'AAA'
        assert self.seq[-6:-3] == 'TGC'

    def test_contains(self):
        assert 'ATGCAT' in self.seq
        assert 'ATGCATGCATGCAAA' in self.seq
        assert 'ATGCGG' not in self.seq

    def test_str(self):
        assert str(self.seq) == 'ATGCATGCATGCAAA'


class TestCodonSequence:
    def setup(self):
        self.seq = CodonSequence('test', 'ATGCATGCATGCAAA')

    def test_len(self):
        assert len(self.seq) == 5

    def test_getitem(self):
        assert self.seq[0] == 'ATG'
        assert self.seq[1] == 'CAT'
        assert self.seq[4] == 'AAA'
        assert self.seq[-1] == 'AAA'

        assert self.seq[0:1] == 'ATG'
        assert self.seq[0:2] == 'ATGCAT'
        assert self.seq[3:4] == 'TGC'
        assert self.seq[3:] == 'TGCAAA'
        assert self.seq[-2:] == 'TGCAAA'
        assert self.seq[-3:-1] == 'GCATGC'

    def test_contains(self):
        assert 'ATGCAT' in self.seq
        assert 'ATGCATGCATGCAAA' in self.seq
        assert 'ATGCGG' not in self.seq

    def test_str(self):
        assert str(self.seq) == 'ATG CAT GCA TGC AAA', print(str(self.seq))


class TestAlignmentEmpty:
    def setup(self):
        self.aln = Alignment('test')

    def test_add_sequence_obj(self):
        seq_obj = NuclSequence('test', 'ATGCATGCATGCAAA')
        self.aln.add_sequence_obj(seq_obj)
        assert len(self.aln._records) == 1
        assert len(self.aln._records_lookup_d) == 1
        assert len(self.aln._aln_matrix) == 1

    def test_add_sequence(self):
        self.aln.add_sequence('test', 'ATGCATGCATGCAAA', 'nucleotide')
        assert len(self.aln._records) == 1
        assert len(self.aln._records_lookup_d) == 1
        assert len(self.aln._aln_matrix) == 1
