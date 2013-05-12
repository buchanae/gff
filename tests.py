from collections import OrderedDict
import itertools
from StringIO import StringIO

import nose
from nose.tools import eq_, ok_, raises

import gff


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

VALID_COLUMNS = [
    'seqid',
    'source',
    'type',
    '1',
    '2',
    '1.1',
    '+',
    '0',
    'one=two;three=four',
]

VALID = '\t'.join(VALID_COLUMNS)


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

    expected_attributes = OrderedDict([('one', 'two'), ('three', 'four')])
    eq_(expected_attributes, a.attributes)


def test_feature_to_GFF_string():
    a = gff.GFF.from_string(VALID)
    eq_(VALID, str(a))


@raises(gff.GFF.ParseError)
def test_invalid_columns():
    missing_col_str = '\t'.join(VALID_COLUMNS[:-1])
    a = gff.GFF.from_string(missing_col_str)


@raises(gff.GFF.ParseError)
def test_invalid_start():
    cols = list(VALID_COLUMNS)
    cols[3] = 'foo'
    invalid_start_str = '\t'.join(cols)
    a = gff.GFF.from_string(invalid_start_str)


@raises(gff.GFF.ParseError)
def test_invalid_end():
    cols = list(VALID_COLUMNS)
    cols[4] = 'foo'
    invalid_start_str = '\t'.join(cols)
    a = gff.GFF.from_string(invalid_start_str)


def test_Reader():
    fh = StringIO(RAW_GFF)
    reader = gff.Reader(fh)
    a = reader.next()
    eq_('Chr1', a.attributes['ID'])


def test_Tree_Node():
    t = gff.Tree()

    a = t.Node()
    b = t.Node()
    c = t.Node()
    d = t.Node()

    c.parent = a
    d.parent = b

    eq_(a.children, [c])
    eq_(b.children, [d])


def test_GFFTree():
    fh = StringIO(RAW_GFF)
    reader = gff.Reader(fh)
    records = list(reader)

    t = gff.GFFTree(records)

    all_recs_prefix = list(t.walk())

    expected = [
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
        records[6],
        records[8],
    ]

    for a, b in itertools.izip_longest(all_recs_prefix, expected):
        print repr(a).ljust(50), repr(b)

    eq_(all_recs_prefix, expected)


if __name__ == '__main__':
    nose.main()
