from nose.tools import eq_

import files
import gff


def test_reader():
    p = files.path('gff')
    r = gff.Reader(p)

    # test that the Reader iterator can be reused
    for i in xrange(2):
        features = list(r)
        eq_({
            'seqid': 'Chr1',
            'source': 'TAIR10',
            'type': 'chromosome',
            'start': 1,
            'end': 30427671,
            'score': '.',
            'strand': '.',
            'phase': '.',
            'raw_attributes': 'ID=Chr1;Name=Chr1',
            'ID': 'Chr1',
            'Name': 'Chr1',
        }, features[0].__dict__)
        eq_({
            'seqid': 'Chr1',
            'source': 'TAIR10',
            'type': 'gene',
            'start': 3631,
            'end': 5899,
            'score': '.',
            'strand': '+',
            'phase': '.',
            'raw_attributes': 'ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010',
            'ID': 'AT1G01010',
            'Note': 'protein_coding_gene',
            'Name': 'AT1G01010',
        }, features[1].__dict__)
