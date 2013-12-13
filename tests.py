from collections import OrderedDict
import itertools
from StringIO import StringIO

import nose
from nose.tools import eq_, ok_, raises, assert_raises

import gff
from gff.parser import ParseError
from gff.tree import build_tree


RAW_GFF = """
Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1
Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1	TAIR10	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
Chr1	TAIR10	protein	3760	5630	.	+	.	ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	five_prime_UTR	3631	3759	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	3760	3913	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	3996	4276	.	+	2	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1
""".strip()

def make_gff(seqid='seqid', source='source', type='type', start='1', end='2',
             score='1.1', strand='+', phase='0',
             attrs='ID=two;Parent=four,five'):

    return '\t'.join([seqid, source, type, start, end, score,
                      strand, phase, attrs])

VALID = make_gff()


def test_feature_from_string():
    a = gff.GFF.from_string(VALID)

    eq_('seqid', a.seqid)
    eq_('source', a.source)
    eq_('type', a.type)
    eq_(1, a.start)
    eq_(2, a.end)
    eq_(1.1, a.score)
    eq_('+', a.strand)
    eq_(0, a.phase)

    expected_attributes = OrderedDict([('ID', 'two'), ('Parent', 'four,five')])
    eq_(expected_attributes, a.attributes)


def test_feature_to_GFF_string():
    a = gff.GFF.from_string(VALID)
    eq_(VALID, str(a))


def test_ID_property():
    a = gff.GFF.from_string(VALID)
    eq_(a.ID, 'two')
    a.ID = 'foo'
    eq_(a.ID, 'foo')
    updated = make_gff(attrs='ID=foo;Parent=four,five')
    eq_(str(a), updated)


def test_parent_IDs_property():
    a = gff.GFF.from_string(VALID)
    eq_(a.parent_IDs, ['four', 'five'])
    a.parent_IDs = ['six', 'seven']
    updated = make_gff(attrs='ID=two;Parent=six,seven')
    eq_(str(a), updated)
    

def test_parent_ID_property():
    s = make_gff(attrs='Parent=foo')
    a = gff.GFF.from_string(s)
    eq_(a.parent_ID, 'foo')
    a.parent_ID = 'bar'
    updated = make_gff(attrs='Parent=bar')
    eq_(str(a), updated)

    with assert_raises(gff.GFF.MultipleParents):
        s = make_gff(attrs='Parent=foo,bar')
        a = gff.GFF.from_string(s)
        a.parent_ID


@raises(ParseError)
def test_invalid_columns():
    s = '\t'.join(['.'] * 8)
    a = gff.GFF.from_string(s)


@raises(ParseError)
def test_invalid_start():
    s = make_gff(start='foo')
    a = gff.GFF.from_string(s)


@raises(ParseError)
def test_invalid_end():
    s = make_gff(end='foo')
    a = gff.GFF.from_string(s)


def test_from_file():
    fh = StringIO(RAW_GFF)
    reader = gff.GFF.from_file(fh)
    a = reader.next()
    eq_('Chr1', a.attributes['ID'])


def test_tree():
    fh = StringIO(RAW_GFF)
    reader = gff.GFF.from_file(fh)
    records = list(reader)

    t = build_tree(records)

    all_recs = [node.record for node in t.walk()]

    expected = [
        # root
        None,

        # Chromosome
        records[0],

        records[1],
        records[2],
        records[4],
        records[5],
        records[6],
        records[7],
        records[8],
        records[9],

        # Note this is out of order because TAIR doesn't seem to use the Parent
        # attribute correctly for protein records.
        records[3],
        # Notice that the CDS nodes are duplicated because they have multiple
        # parents.
        records[6],
        records[8],
    ]

    for a, b in itertools.izip_longest(all_recs, expected):
        print repr(a).ljust(50), repr(b)

    eq_(all_recs, expected)


if __name__ == '__main__':
    nose.main()
