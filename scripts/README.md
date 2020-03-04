# Working with FASTA files 
The scripts here use **mrgali** (pyali) but can pre- and post-process files to make it easier to work directly with fasta files.

`python raise_aln.py -in <fasta file> -out <output file> -seq <sequence name>` 

will move the alignment specified by the -seq tag to the top of a fasta file.


`python align_merger.py -in <list with fasta files> -out <output file> -width <width of lines in output fasta> -ref <name of reference sequence>`

takes a list of alignment fasta files that share one common sequence (-ref) and merges them into a single multiple sequence alignment. See *tests/examples* for an example. 

Tests can be used by running `pytest -v`.

If you're using the scripts separately from the repo, make sure that you've installed *pyali* and that you've changed the dependency, as specified, in *align_merger.py*