import pytest

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
    from pyali.mrgali import *
    a = Alignment.from_reference(refs)
    a.merge(0, alis[0])
    a.merge(1, alis[1])
    a.merge(2, alis[2])
    assert str(a).split('\n') == msa

