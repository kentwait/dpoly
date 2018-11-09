import re
from collections import namedtuple
import numpy as np
from copy import deepcopy

SequenceAnnotation = namedtuple('SequenceAnnotation',
                                'name, description, seq_type')

class Sequence(object):
    def __init__(self, name, sequence, description=None, seq_type=None):
        self.name = name
        self.description = description
        self.seq_type = seq_type
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, i):
        return self.sequence[i]

    def __iter__(self):
        return iter(self.sequence)

    def __contains__(self, s):
        return s in self.sequence

    def __repr__(self):
        return self.sequence

    def fasta_format(self, line_width=None):
        string = '>{id} {description}\n'.format(
            id=self.name,
            description=self.description,
        )
        if line_width:
            last = 0
            for i in range(0, len(self), line_width):
                string += self.sequence[i:i+line_width]
                last = i
            string += self.sequence[last:]
            return string
        string += self.sequence
        return string

class NuclSequence(Sequence):
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description, 
                         seq_type='nucleotide')

class ProtSequence(Sequence):
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description, 
                         seq_type='protein')

class CodonSequence(Sequence):
    def __init__(self, name, sequence, description=None):
        super().__init__(name, sequence, description=description, 
                         seq_type='codon')

class Alignment(object):
    def __init__(self, name, description=None, aln_type=None):
        self.name = name
        self.description = description
        self.aln_type = aln_type
        self._records = []
        self._records_lookup_d = dict()
        self.aln_matrix = np.array()
        self.filters = dict()

    def add_sequence(self, sequence_obj):
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
        self.aln_matrix = np.vstack((self.aln_matrix, seq_array))

    def filter_sites(self, *filter_names, exclude_char='X'):
        keep_coords = set(range(self.aln_matrix.shape[0]))
        for filt in filter_names:
            if filt in self.filters.keys():
                coords = set(self.filters[filt].coords(
                    exclude_char, inverse=True))
                keep_coords = keep_coords.intersection(coords)
        new_aln = deepcopy(self)
        new_aln.aln_matrix = self.aln_matrix[:, coords]
        return new_aln

    def filter_sequences(self, *sequence_names):
        positions = []
        new_records = []
        for name in sequence_names:
            if name in self._records_lookup_d.keys():
                i = self._records_lookup_d[name]
                positions.append(i)
                new_records.append(self._records[i])

        new_aln = deepcopy(self)
        new_aln.aln_matrix = self.aln_matrix[positions]
        new_aln._records = new_records
        new_aln._records_lookup_d = {k:v for k, v in self._records_lookup_d
                                     if k in sequence_names}
        return new_aln

    def use_all_filters(self, exclude_char='X'):
        return self.filter_sites(*list(self.filters.keys()),
                                 exclude_char=exclude_char)

    def __len__(self):
        return len(self.aln_matrix.shape[-1])

    def __getitem__(self, i):
        if isinstance(i, int):
            return self.aln_matrix[:, i]
        elif isinstance(i, str):
            if i in self._records_lookup_d.keys():
                pos = self._records_lookup_d[i]
                return self.aln_matrix[pos]
        return None

    def __iter__(self):
        return iter(self.aln_matrix.transpose())

    # def __repr__(self):
    #     pass


class NuclAlignment(Alignment):
    def __init__(self, name, description=None):
        super().__init__(name, description=description, aln_type='nucleotide')


class ProtAlignment(Alignment):
    def __init__(self, name, description=None):
        super().__init__(name, description=description, aln_type='protein')


class CodonAlignment(Alignment):
    def __init__(self, name, description=None):
        super().__init__(name, description=description, aln_type='codon')

    def __len__(self):
        return len(self.aln_matrix.shape[-1]) / 3

    def __getitem__(self, i):
        if isinstance(i, int):
            x = i/3
            return self.aln_matrix[:, x:x+3]
        elif isinstance(i, str):
            if i in self._records_lookup_d.keys():
                pos = self._records_lookup_d[i]
                return self.aln_matrix[pos]
        return None

    def __iter__(self):
        return iter(self.aln_matrix.transpose())


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
            if len(char_list) == 0 or char_list[-1] != c:
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
