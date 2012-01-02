from nose.tools import eq_

import files
import gff


def test_reader():
    p = files.path('gff')
    r = gff.Reader(p)
    eq_([
        'Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1',
        'Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010'], list(r))


def test_feature_from_string():
    p = files.path('gff')
    r = gff.Reader(p)
    r = list(r)
    a = gff.Feature.from_string(r[0])

    eq_('Chr1', a.seqid)
    eq_('TAIR10', a.source)
    eq_('chromosome', a.type)
    eq_(1, a.start)
    eq_(30427671, a.end)
    eq_(30427671, a.length)
    eq_('.', a.score)
    eq_('.', a.strand)
    eq_('.', a.phase)
    eq_('ID=Chr1;Name=Chr1', a.raw_attributes)
    eq_('Chr1', a.attributes['ID'])
    eq_('Chr1', a.attributes['Name'])
