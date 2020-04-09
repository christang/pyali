# Author: Antoniya A. Aleksandrova
# Language: Python 3
# Description: Takes pairwise alignments and outputs a multiple alignment
# Usage: python align_merger.py -in <list with fasta files> -out <output file> -width <width of lines in output fasta> -ref <name of reference sequence>

import argparse
import glob
import os
import sys

try:
    from pyali.mrgali import Alignment
except:
    print("[ERROR]: pyali has not been installed. To install, run `python setup.py install` from project directory.")
from raise_aln import raise_seq


def align_merger(file_list, outname, width, reference_seq):
    refs, alis, alis_names = [], [], []
    for f in file_list:
        if reference_seq != '':
            raise_seq(f, f + '.tmp', reference_seq)
            alignment = open(f + '.tmp', 'r')
        else:
            alignment = open(f, 'r')
        flag = 0
        sequence1 = ""
        sequence2 = ""
        alis_element = []
        for line in alignment:
            if line.startswith(">") and flag == 2:
                alis_element.append(sequence2)
                alis_names.append(name2)
                sequence2 = ""
                name2 = line.strip()
            if line.startswith(">") and flag == 1:
                flag = 2
                alis_element.append(sequence1)
                name2 = line.strip()
            if line.startswith(">") and flag == 0:
                flag = 1
                name1 = line.strip()
            if not line.startswith(">") and flag == 1:
                sequence1 = sequence1 + line.strip().upper()
            if not line.startswith(">") and flag == 2:
                sequence2 = sequence2 + line.strip().upper()
        alignment.close()
        if reference_seq != '':
            os.remove(f + '.tmp')
        alis_element.append(sequence2)
        alis.append(alis_element)
        alis_names.append(name2)

    refs = [''.join([s for s in seqs[0] if s != '-']) for seqs in alis]
    if refs.count(refs[0]) != len(refs):
        print("The reference sequences in all the provided alignments are not identical.")
        for i, r in enumerate(refs[1:]):
            for j, s in enumerate(refs[0]):
                if s!=refs[i+1][j]:
                    print(file_list[0] +": (" + str(j) + "," + s + "), " + file_list[i+1] + ": " + r[j])
        raise SystemExit("References need to be the same to proceed.")

    a = Alignment.from_reference(refs)
    for i in range(len(alis)):
        a.merge(i, alis[i])

    flds = str(a).split('\n')

    aligned_list = []
    out = open(outname, 'w')
    for i, ln in enumerate(flds):
        if i == 0:
            s = ln[ln.index(':') + 2:]
            out.write(name1 + '\n')
            aligned_list.append((name1, s))
            while len(s) > 0:
                out.write(s[:width] + '\n')
                s = s[width:]
        if i >= len(refs):
            s = ln[ln.index(':') + 2:]
            out.write(alis_names[i - len(refs)] + '\n')
            aligned_list.append((alis_names[i - len(refs)], s))
            while len(s) > 0:
                out.write(s[:width] + '\n')
                s = s[width:]
    out.close()
    return aligned_list


if __name__ == "__main__":

    # Remove previously generated merged alignments
    if os.path.isfile('merged_alignments.fasta'):
        os.remove('merged_alignments.fasta')

    # Parse user input and set defaults
    parser = argparse.ArgumentParser()

    parser.add_argument('-in', '--list_of_fasta_files', nargs='?')  # file with a list of fasta files (one per line)
    parser.add_argument('-out', '--out', nargs='?')  # output file
    parser.add_argument('-width', '--width', nargs='?')  # fasta line width
    parser.add_argument('-ref', '--reference', nargs='?')  # name of the reference sequence

    parser.set_defaults(list_of_fasta_files=glob.glob('*.fasta'))
    parser.set_defaults(out='merged_alignments.fasta')
    parser.set_defaults(width='72')
    parser.set_defaults(reference='')

    parsed = parser.parse_args()
    if type(parsed.__dict__['list_of_fasta_files']) == list:
        file_list = parsed.__dict__['list_of_fasta_files']
        print("[INFO]: Pairwise alignments will be taken from all fasta files in current directory.")
    else:
        with open(parsed.__dict__['list_of_fasta_files'], 'r') as f:
            file_list = f.read().strip().split('\n')
        print("[INFO]: Pairwise alignment files read from " + parsed.__dict__['list_of_fasta_files'] + ".")
    print("[INFO]: The following files were read: " + ', '.join(file_list))
    if len(file_list) == 0:
        raise SystemExit("No fasta files available for alignment.")
    outname = parsed.__dict__['out']
    print("[INFO]: Merged alignment will be written to " + outname + ".")
    width = int(parsed.__dict__['width'])
    ref_seq = parsed.__dict__['reference']

    align_merger(file_list, outname, width, ref_seq)

    print("[INFO]: Done.")
