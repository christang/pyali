import numpy as np
import pandas as pd


# a target, reference alignment of interest
test_refs = [
    'abcde-',
    '-bcdef',
]

# a set of alignments based on one of the reference sequences
# to merge onto the reference alignment
test_alis = [
    [
        'ab-cde',
        '-bbcd-',
    ], [
        'bcdef',
        '-cde-',
    ]
]


def sequence_tuples(seq):
    return [(a, i) for i, a in enumerate(seq)]

def alignment_frame(ali, base=0):
    df = pd.DataFrame(map(list, ali)).transpose()
    df.columns = range(base, base+len(ali))
    return df


class Alignment(object):

    GAP = '-'

    @classmethod
    def from_reference(cls, refs):
        obj = cls()
        obj.n_refs = len(refs)
        obj.alignment = pd.DataFrame([sequence_tuples(s) for s in refs]).transpose()
        return obj

    def merge(self, index, ali):
        self.alignment = self.merge_on(index, ali)

    def merge_on(self, index, ali):
        ali_index = self.get_indexed(index, ali[0])
        ali_frame = alignment_frame(ali[1:], 1 + self.alignment.columns.max())
        ali_frame[index] = ali_index
        return self.alignment.merge(ali_frame, on=[index], how='left').fillna(Alignment.GAP)

    def get_indexed(self, index, seq):
        ref_ungapped = [(a, i) for a, i in self.alignment[index] if a != Alignment.GAP]
        np_seq = np.array([None for a in seq])
        np_seq[[i for i, a in enumerate(seq) if a != Alignment.GAP]] = ref_ungapped
        return np_seq

    def build_str(self):
        i_refs = range(self.n_refs)
        i_seqs = [i for i in self.alignment.columns if i not in i_refs]
        s = []
        for col in i_refs:
            s.append('%d: %s' % (col, ''.join(self.alignment[col].apply(lambda e: e[0]))))
        for col in i_seqs:
            s.append('%d: %s' % (col, ''.join(self.alignment[col])))
        return '\n'.join(s)

    def __str__(self):
        return self.build_str()


# test code
a = Alignment.from_reference(test_refs)
a.merge(0, test_alis[0])
a.merge(1, test_alis[1])
print a

