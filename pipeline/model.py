from copy import deepcopy
from collections import namedtuple
import numpy as np


SequenceAnnotation = namedtuple('SequenceAnnotation',
                                'name, description, seq_type')

class Sequence(object):
    """Represents a biological sequence of characters.

    The Sequence object is indexable using [] and can be
    sliced like a string or a list.

    """
    def __init__(self, name, sequence, description=None, seq_type=None):
        """Creates a new instance of the Sequence object.

        Parameters
        ----------
        name : str
            Name of the sequence. This corresponds to the text found in the
            identifier line of a FASTA-formatted file.
        sequence : str
        description : str, optional
            Description or other miscellenous information about the sequence.
            In the FASTA format, this is found in the identifier line, but is
            separated from the identifier by a whitespace.
        seq_type: str, optional
            Indicates whether the sequence is a nucleotide, protein
            (amino acid), or codon sequence.

        """
        self.name = name
        self.description = description
        self.seq_type = seq_type
        self._sequence = sequence

    @property
    def i(self):
        """Always returns a read-only array of single characters,
        even for codon sequences.
        """
        return np.array(self._sequence)

    @property
    def sequence(self):
        """Returns a read-only string representation of the
        sequence.
        """
        return self._sequence

    def fasta_format(self, line_width=None):
        """Output the sequence as a FASTA-formatted string.

        Parameters
        ----------
        line_width : int
            Number of characters per line.

        Returns
        -------
        str

        """
        string = '>' + self.name
        if self.description:
            string += ' ' + self.description
        string += '\n'

        if line_width:
            last = 0
            for i in range(line_width, len(self), line_width):
                string += self._sequence[i-line_width:i] + '\n'
                last = i
            string += self._sequence[last:]
            return string
        string += self._sequence
        return string

    def __len__(self):
        return len(self._sequence)

    def __getitem__(self, i):
        return self._sequence[i]

    def __iter__(self):
        return iter(self._sequence)

    def __contains__(self, x):
        return x in self._sequence

    def __str__(self):
        return self._sequence


class NuclSequence(Sequence):
    """Represents a nucleotide sequence.
    """
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         seq_type='nucleotide')

class ProtSequence(Sequence):
    """Represents a protein or amino acid sequence.
    """
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         seq_type='protein')

class CodonSequence(Sequence):
    """Represents a codon-based coding sequence.
    """
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         seq_type='codon')

    def __len__(self):
        return int(len(self._sequence) / 3)

    def __getitem__(self, i):
        if isinstance(i, int):
            if i >= 0:
                i = i * 3
                return self._sequence[i:i+3]
            i = len(self._sequence) + 1 - 1
            return self._sequence[i-3:i]
        elif isinstance(i, slice):
            start = i.start * 3
            end = (i.stop * 3) if isinstance(i.stop, int) \
                  else len(self._sequence)
            return self._sequence[start:end]
        return IndexError()

    def __iter__(self):
        return (self._sequence[i:i+3] for i in range(0, len(self._sequence), 3))

    def __str__(self):
        return ' '.join(list(self))

class Alignment(object):
    """Represents an alignment of biological sequences

    The Alignment object is indexable using [] like a string or a list.
    If the key used is an object, it acts like a column-wise positional index.
    But if the key used is a string, it retrieves like a dictionary key and
    looks up the list of sequence names for a match.

    """
    def __init__(self, name, description=None, aln_type=None):
        """Creates a new alignment.

        Sequences cannot be added directly at instantiation.
        Use the `add_sequence_obj` or `add_sequence` methods to add
        sequences into the alignment.

        Parameters
        ----------
        name : str
        description : str, optional
        aln_type : str, optional

        """
        self.name = name
        self.description = description
        self.aln_type = aln_type
        self._records = []
        self._records_lookup_d = dict()
        self._aln_matrix = np.array([])
        self.filters = dict()

    @property
    def sequences(self):
        """Returns a read-only copy of the alignment as a list of
        strings.
        """
        return [s for s in self]

    @property
    def i(self):
        """Returns a read-only copy of the alignment matrix.
        """
        return self._aln_matrix

    def add_sequence_obj(self, sequence_obj):
        """Adds a Sequence object containing a single aligned sequence
        to the alignment.

        Assumes that the sequence to be added is a Sequence object.
        If adding a sequence string, use the `add_sequence` method
        instead.

        Parameters
        ----------
        sequence_obj : Sequence

        See also
        --------
        add_sequence

        """
        assert sequence_obj.name not in self._records_lookup_d.keys()
        a = SequenceAnnotation(sequence_obj.name,
                               sequence_obj.description,
                               sequence_obj.seq_type)
        if not self.aln_type:
            self.aln_type = a.seq_type
        assert self.aln_type == a.seq_type
        self._records_lookup_d[a.name] = len(self._records)
        self._records.append(a)
        seq_array = np.array(list(sequence_obj.sequence))
        if self._aln_matrix.shape[-1] == 0:
            self._aln_matrix = np.array([seq_array,])
        else:
            self._aln_matrix = np.vstack((self._aln_matrix, seq_array))

    def add_sequence(self, name, sequence, seq_type, description=None):
        """Adds a single alignmed sequence to the alignment.

        In this method, the sequence is a string and not yet wrapped
        as a Sequence object. For adding sequence objects, use the
        `add_sequence_obj` method instead.

        Parameters
        ----------
        name : str
        sequence : str
        seq_type : str
        description : str, optional

        See aslo
        --------
        add_sequence_obj

        """
        s = Sequence(name, sequence, description=description, seq_type=seq_type)
        self.add_sequence_obj(s)

    def filter_sites(self, *filter_names, exclude_char='X'):
        """Filters out sites in the alignmnet using a given list of filters.

        Assumes that all filter tracks use the same excluding character to mark
        sites. The default excluding character is "X".

        Parameters
        ----------
        filter_names: str or Filter object
            Names of registered filters or Filter objects to be used
            to filter the alignment.
        exclude_char : str
            Filter track character that marks sites to be excluded
            from the returned alignment.

        Returns
        -------
        Alignment
            New alignment object with excluded sites removed from the
            alignment. The number of entries in the alignment does not change.

        """
        keep_coords = set(range(self._aln_matrix.shape[0]))
        for filt in filter_names:
            coords = set()
            if isinstance(filt, str) and filt in self.filters.keys():
                coords = set(self.filters[filt].coords(
                    exclude_char, inverse=True))
            elif isinstance(filt, Filter):
                coords = filt.coords(exclude_char, inverse=True)
            else:
                raise ValueError()
            keep_coords = keep_coords.intersection(coords)
        new_aln = deepcopy(self)
        new_aln._aln_matrix = self._aln_matrix[:, coords]  # pylint: disable=protected-access
        return new_aln

    def filter_sequences(self, *sequence_names):
        """Filters the alignmnet by returning only the specified entries.

        This does not change the number of sites in the alignment.

        Parameters
        ----------
        sequence_names
            Filter track character that marks sites to be excluded
            from the returned alignment

        Returns
        -------
        Alignment
            New alignment object with excluded sites removed from the
            alignment. The number of entries in the alignment does not change.

        """
        positions = []
        new_records = []
        for name in sequence_names:
            if name in self._records_lookup_d.keys():
                i = self._records_lookup_d[name]
                positions.append(i)
                new_records.append(self._records[i])

        new_aln = deepcopy(self)
        new_aln._aln_matrix = self._aln_matrix[positions]  # pylint: disable=protected-access
        new_aln._records = new_records  # pylint: disable=protected-access
        new_aln._records_lookup_d = {k:v for k, v in self._records_lookup_d # pylint: disable=protected-access
                                     if k in sequence_names}
        return new_aln

    def use_all_filters(self, exclude_char='X'):
        return self.filter_sites(*list(self.filters.keys()),
                                 exclude_char=exclude_char)

    def __len__(self):
        return len(self._aln_matrix.shape[-1])

    def __getitem__(self, i):
        if isinstance(i, int):
            return self._aln_matrix[:, i]
        elif isinstance(i, slice):
            return self._aln_matrix[:, i]
        elif isinstance(i, str):
            if i in self._records_lookup_d.keys():
                pos = self._records_lookup_d[i]
                return self._aln_matrix[pos]
        return IndexError()

    def __iter__(self):
        return iter(self._aln_matrix.transpose())

    # def __repr__(self):
    #     pass


class NuclAlignment(Alignment):
    def __init__(self, name, description=None):
        super().__init__(name, description=description, aln_type='nucleotide')

    def add_sequence_obj(self, sequence_obj: NuclSequence):
        assert sequence_obj.seq_type == 'nucleotide'
        super().add_sequence_obj(sequence_obj)


class ProtAlignment(Alignment):
    def __init__(self, name, description=None):
        super().__init__(name, description=description, aln_type='protein')

    def add_sequence_obj(self, sequence_obj: NuclSequence):
        assert sequence_obj.seq_type == 'nucleotide'
        super().add_sequence_obj(sequence_obj)


class CodonAlignment(Alignment):
    def __init__(self, name, description=None):
        super().__init__(name, description=description, aln_type='codon')

    def add_sequence_obj(self, sequence_obj: NuclSequence):
        assert sequence_obj.seq_type == 'nucleotide'
        super().add_sequence_obj(sequence_obj)

    def __len__(self):
        return int(len(self._aln_matrix.shape[-1]) / 3)

    def __getitem__(self, i):
        if isinstance(i, int):
            x = int(i/3)
            return self._aln_matrix[:, x:x+3]
        elif isinstance(i, slice):
            start = i.start * 3
            end = (i.stop * 3) + 3
            return self._aln_matrix[:, start:end]
        elif isinstance(i, str):
            if i in self._records_lookup_d.keys():
                pos = self._records_lookup_d[i]
                return self._aln_matrix[pos]
        return None

    def __iter__(self):
        return iter(self._aln_matrix)


class Filter(object):
    def __init__(self, name, sequence, description=None, chars=None):
        self.pos_list = ()
        self.char_list = ()
        self.name = name
        self.description = description
        if chars:
            self.check_sequence(sequence, chars)
        self.encode(sequence)

    def encode(self, sequence):
        pos_list = []
        char_list = []
        last = 0
        for i, c in enumerate(sequence):
            if not char_list or char_list[-1] != c:
                pos_list.append(str(i))
                char_list.append(c)
                last = 1
        pos_list.append(last+1)
        self.pos_list = tuple(pos_list)
        self.char_list = tuple(char_list)

    @staticmethod
    def check_sequence(sequence, chars):
        for c in sequence:
            if c not in chars:
                raise Exception(
                    'character {} is not one of the allowed ' \
                    'characters in the sequence ({})'.format(
                        c, ', '.join(chars))
                    )

    def full_sequence(self):
        len_list = [int(self.pos_list[i+1]) - int(self.pos_list[i])
                    for i in range(len(self.pos_list)-1)]
        return ''.join(
            [char * length for length, char in zip(len_list, self.char_list)]
        )

    def encoded_sequence(self):
        return ''.join(zip(list(map(str, self.pos_list)), self.char_list))

    def coords(self, char, inverse=False):
        first = 0
        coords_list = []
        for last, c in zip(self.pos_list[1:], self.char_list):
            if (char == c and not inverse) or (char != c and inverse):
                coords_list += list(range(first, last))
            first = last
        return np.array(coords_list)

    def filter(self, sequence, exclude_char):
        assert len(sequence) == len(self)
        seq_array = np.array(sequence)
        coords = self.coords(exclude_char, inverse=True)
        return ''.join(seq_array[coords])

    def mask(self, sequence, exclude_char, mask_char='_'):
        assert len(sequence) == len(self)
        seq_array = np.array(sequence)
        coords = self.coords(exclude_char, inverse=False)
        seq_array[coords] = mask_char
        return ''.join(seq_array)

    def __repr__(self):
        return self.name + ':' + self.encoded_sequence()

    def __len__(self):
        return self.pos_list[-1]

class ConsAlignFilter(Filter):
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         chars=('-', 'X'))

    def consistent_sites(self, sequence):
        return self.filter(sequence, '-')

    def inconsistent_sites(self, sequence):
        return self.filter(sequence, 'X')

    def mask_consistent_sites(self, sequence, mask_char='_'):
        return self.mask(sequence, '-', mask_char=mask_char)

    def mask_inconsistent_sites(self, sequence, mask_char='_'):
        return self.mask(sequence, 'X', mask_char=mask_char)

    def consistent_site_coords(self):
        return self.coords('-')

    def inconsistent_site_coords(self):
        return self.coords('X')

class GapFilter(Filter):
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description,
                         chars=('-', 'X'))

    def remove_gaps(self, sequence):
        return self.filter(sequence, 'X')
