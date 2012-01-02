from nose.tools import eq_, raises

import files
import gff


def test_reader():
    p = files.path('gff')
    r = gff.reader(p)
    eq_([
        'Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1',
        'Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010'], list(r))


def test_feature_from_string():
    p = files.path('gff')
    r = gff.reader(p)
    a = gff.Feature.from_string(r.next())

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

def test_feature_to_string():
    p = files.path('gff')
    r = gff.reader(p)
    l = r.next()
    a = gff.Feature.from_string(l)
    eq_(l, str(a))

@raises(gff.InvalidGFFString)
def test_invalid_columns():
    p = files.path('invalid')
    r = gff.reader(p)
    a = gff.Feature.from_string(r.next())

@raises(gff.InvalidGFFString)
def test_invalid_start():
    p = files.path('invalid')
    r = gff.reader(p)
    r.next()
    a = gff.Feature.from_string(r.next())

@raises(gff.InvalidGFFString)
def test_invalid_end():
    p = files.path('invalid')
    r = gff.reader(p)
    r.next()
    r.next()
    a = gff.Feature.from_string(r.next())
