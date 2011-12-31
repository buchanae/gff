from nose.tools import eq_

import files
import gff


def test_parse_attributes():
    eq_({'Foo': 'Bar', 'BaZ': 'bat'}, gff.parse_attributes('Foo=Bar;BaZ=bat'))
    eq_({}, gff.parse_attributes(''))

def test_reader():
    p = files.path('gff')
    lines = list(gff.reader(p))
    eq_({
        'seqid': 'Chr1',
        'source': 'TAIR10',
        'type': 'chromosome',
        'start': 1,
        'end': 30427671,
        'score': '.',
        'strand': '.',
        'phase': '.',
        'attributes': 'ID=Chr1;Name=Chr1',
        'ID': 'Chr1',
        'Name': 'Chr1',
    }, lines[0])
    eq_({
        'seqid': 'Chr1',
        'source': 'TAIR10',
        'type': 'gene',
        'start': 3631,
        'end': 5899,
        'score': '.',
        'strand': '+',
        'phase': '.',
        'attributes': 'ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010',
        'ID': 'AT1G01010',
        'Note': 'protein_coding_gene',
        'Name': 'AT1G01010',
    }, lines[1])
