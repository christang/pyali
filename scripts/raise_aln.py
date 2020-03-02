# Author: Antoniya A. Aleksandrova
# Language: Python 3.5/2.7
# Description: Move a specified alignment to the top of a fasta file 
# Usage: python raise_aln.py -in <fasta file> -out <output file> -seq <sequence name>  

import os
import sys
import argparse

def raise_seq(infile, outfile, seqn):
	aligns = []
	name = ""
	seq = ""
	if '>' not in seqn:
		seqn = '>' + seqn
	f = open(infile, 'r')
	for line in f:
		if line.startswith(">"):
			if seq != "":
				aligns.append((name, seq))
			name = line.strip()
			seq = ""
		elif not line.startswith("#"):
			seq = seq + line
	aligns.append((name, seq))
	f.close()
	
	index = [x for x, y in enumerate(aligns) if y[0] == seqn] # locate the top sequence
	if len(index)==0:
		raise SystemExit(seqn + " cannot be located in " + infile)
	else:
		index = index[0]
	out = open(outfile, 'w')
	out.write(aligns[index][0] + '\n' + aligns[index][1].strip())
	for a in aligns:
		if a[0]!=seqn:
			out.write('\n' + a[0] + '\n' + a[1].strip())
	out.close()



if __name__ == "__main__":
	if len(sys.argv)<3:
		raise SystemExit("Usage: python raise_fasta.py -in <fasta file> -seq <sequence name>\n\tOptional inputs: -out <output file>")

	parser = argparse.ArgumentParser()

	parser.add_argument('-in', '--fasta_file', nargs='?')
	parser.add_argument('-out', '--out', nargs='?') # output file
	parser.add_argument('-seq', '--seqname', nargs='?') # name of the sequence that should be moved to the top

	parser.set_defaults(out = '_raised.fasta')
	parsed = parser.parse_args()
	infile = parsed.__dict__['fasta_file']
	filename, file_extension = os.path.splitext(infile)
	outfile = parsed.__dict__['out']
	if outfile == '_raised.fasta':
		outfile = filename + '_raised' + file_extension
	seqnanme = '>' + parsed.__dict__['seqname']
	raise_seq(infile, outfile, seqname)