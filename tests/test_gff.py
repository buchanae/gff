from collections import OrderedDict
from tempfile import NamedTemporaryFile

from nose.tools import eq_, raises

from gff import Attributes, Feature


valid = OrderedDict([
    ('seqid', 'SeqID'), ('source', 'SOURCE'), ('type', 'type'), ('start', '23'), 
    ('end', '41'), ('score', '.'), ('strand', '-'), ('phase', '.'), 
    ('attribs', 'ID=Foo;Parent=Bar,Baz;Note=FOO')
])

def mk(**kwargs):
    v = valid.copy()
    v.update(kwargs)
    return '\t'.join(v.values())

valid_str = mk()
valid_b_str = mk(seqid='SeqIDB')
missing_col_str = '\t'.join(valid.values()[:-1])
invalid_start_str = mk(start='foo')
invalid_end_str = mk(end='foo')

dummy_file = NamedTemporaryFile(delete=False)
dummy_file.write('##comment\n')
dummy_file.write('##comment\n')
dummy_file.write(valid_str + '\n')
dummy_file.write(missing_col_str + '\n')
dummy_file.write(invalid_start_str + '\n')
dummy_file.write(invalid_end_str + '\n')
dummy_file.write(valid_b_str + '\n')
dummy_file.close()


def test_reader():
    eq_([valid_str, valid_b_str], [str(x) for x in Feature.from_file(dummy_file.name)])

def test_feature_from_string():
    a = Feature.from_string(valid_str)

    eq_('SeqID', a.seqid)
    eq_('SOURCE', a.source)
    eq_('type', a.type)
    eq_(23, a.start)
    eq_(41, a.end)
    eq_(19, a.length)
    eq_('.', a.score)
    eq_('-', a.strand)
    eq_('.', a.phase)
    eq_('ID=Foo;Parent=Bar,Baz;Note=FOO', str(a.attributes))
    eq_('Foo', a.attributes['ID'])
    eq_('FOO', a.attributes['Note'])
    eq_(['Bar', 'Baz'], a.attributes['Parent'])

def test_feature_to_string():
    a = Feature(*valid.values())
    eq_(valid_str, str(a))

@raises(Feature.ParseError)
def test_invalid_columns():
    a = Feature.from_string(missing_col_str)

@raises(Feature.ParseError)
def test_invalid_start():
    a = Feature.from_string(invalid_start_str)

@raises(Feature.ParseError)
def test_invalid_end():
    a = Feature.from_string(invalid_end_str)

def test_attributes():
    a = Attributes.from_string('ID=Csa1M000010.1.exon2; Parent=Csa1M000010.1')
    b = Attributes([('ID', 'Csa1M000010.1.exon2'), ('Parent', 'Csa1M000010.1')])
    eq_(b, a)
