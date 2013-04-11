from collections import OrderedDict
from StringIO import StringIO

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
    'reference',
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

    eq_('reference', a.reference)
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
