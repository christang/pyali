import pytest
import sys
import os
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[2]) + '/pyali')

testsdir = str(Path(__file__).resolve().parents[1]) + '/tests/'

@pytest.mark.mrgali
@pytest.mark.parametrize('refs,alis,msa', [(['IGE','IGE','IGE'],
                                                  [
                                                      [
                                                          'I-G-E',
                                                          'IRGIS',
                                                          'TT---'
                                                      ], [
                                                          'IGE',
                                                          'TTT'
                                                      ], [
                                                          'I-G-E',
                                                          'F-F--'
                                                      ]
                                                  ],
                                                  ['0: I-G-E', '1: I-G-E', '2: I-G-E',
                                                   '3: IRGIS', '4: TT---', '5: T-T-T', '6: F-F--']
                                                  )])
def test_pyali_Alignmnet(refs, alis, msa):
    from mrgali import Alignment
    a = Alignment.from_reference(refs)
    a.merge(0, alis[0])
    a.merge(1, alis[1])
    a.merge(2, alis[2])
    assert str(a).split('\n') == msa


@pytest.mark.scripts
@pytest.mark.parametrize('infile,outfile,refn', [(testsdir + 'examples/a2.fa', testsdir + 'examples/a2_tmp.fa','template')])
def test_raise_seq(infile, outfile, refn):
    from raise_aln import raise_seq
    assert os.path.isfile(infile) == True, '%s does not exist.'%(infile)
    raise_seq(infile, outfile, refn)
    assert os.path.isfile(outfile), '%s does not exist.'%(outfile)
    with open(outfile) as f:
        content = f.read().splitlines()
    assert content[0] == '>' + refn, 'The sequence of %s was not raised to the top of the alignments.' % (refn)
    os.remove(outfile)

@pytest.mark.scripts
@pytest.mark.parametrize('file_list,outname,width,refn,result', [([testsdir + 'examples/a1.fa', testsdir + 'examples/a2.fa'],
                                                           testsdir + 'examples/a1_a2_merged.fa', 80, 'template',
                                                                         ['>template',
                                                                         '--DERE-T',
                                                                         '>t1',
                                                                         '--CC-CCC',
                                                                         '>t2',
                                                                         '--D-DDDD',
                                                                         '>t3',
                                                                         '--EEE-EE',
                                                                         '>t4',
                                                                         'FFF-FF-F']
                                                                  )])
def test_align_merger(file_list, outname, width, refn, result):
    from align_merger import align_merger
    align_merger(file_list, outname, width, refn)
    assert os.path.isfile(outname), 'No output was produced'
    with open(outname) as f:
        content = f.read().splitlines()
    assert content == result

