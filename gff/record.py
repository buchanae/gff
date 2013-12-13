from __future__ import absolute_import

from gff.parser import Parser
from gff.formatter import Formatter


def add_parser(Parser_cls, *args, **kwargs):
    """
    Class decorator that instantiates Parser_cls and adds to the GFF class.
    """
    def decorator(cls):
        cls.parser = parser = Parser_cls(cls, *args, **kwargs)
        cls.from_string = parser.parse
        cls.from_file = parser.parse_many
        return cls
    return decorator

def add_formatter(Formatter_cls, *args, **kwargs):
    """
    Class decorator that instantiates Formatter_cls and adds to the GFF class.
    """
    def decorator(cls):
        cls.formatter = Formatter_cls(*args, **kwargs)
        return cls
    return decorator


@add_parser(Parser)
@add_formatter(Formatter)
class GFF(object):

    """
    Represents a record from GFF record.

    Most commonly, you'll want to read GFF records from a file

    >>> gff_file_path = 'path/to/annotation.gff'
    >>> file_handle = open(gff_file_path)
    >>> records = GFF.from_file(file_handle)


    You can create records directly too, for example
    "GeneX", defined by "Acme", on chromsome 'Chr1'
    on the forward strand at positions 100-200:

    >>> GFF('Chr1', 'Acme', 'gene', 100, 200, '.', '+', '.', {'ID': 'GeneX'})


    GFF records can also be parsed from strings:

    >>> s = 'Chr1	Acme	gene	100	200	.	+	.	ID=GeneX'
    >>> GFF.from_string(s)


    Helpers are provided for accessing common attributes:

    >>> s = 'Chr1	Acme	gene	100	200	.	+	.	ID=GeneX;Parent=Chr1'
    >>> g = GFF.from_string(s)
    >>> g.ID
    'GeneX'
    >>> g.parent_IDs
    ['Chr1']
    >>> g.parent_ID
    'Chr1'
    >>> g.parent_ID = 'foo'
    >>> print g
    Chr1	Acme	gene	100	200	.	+	.	ID=GeneX;Parent=foo
    >>> g.parent_IDs = ['foo', 'bar']
    >>> print g
    Chr1	Acme	gene	100	200	.	+	.	ID=GeneX;Parent=foo,bar


    Note that GFF.parent_IDs (plural) will always return a list,
    and GFF.parent_ID (singular) will raise an exception if the
    record has multiple parent IDs.

    >>> s = 'Chr1	Acme	gene	100	200	.	+	.	ID=GeneX;Parent=foo,bar'
    >>> g = GFF.from_string(s)
    >>> g.parent_IDs
    ['foo', 'bar']
    >>> g.parent_ID
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "gff/record.py", line 183, in parent_ID
        raise self.MultipleParents()
    gff.record.MultipleParents


    Attributes are stored in a dictionary, which means their original order
    is lost. If you need to maintain ordering, try this:

    >>> from collections import OrderedDict
    >>> GFF.parser.attributes_cls = OrderedDict
    >>> s = 'Chr1	Acme	gene	3631	5899	.	+	.	one=1;two=2'
    >>> g = GFF.from_string(s)
    >>> g.attributes
    OrderedDict([('one', '1'), ('two', '2')])

    Note that using OrderedDict causes a substantial drop in performance
    (~3x slower to read records)
    """

    class MultipleParents(Exception): pass

    # I'm using __slots__ to save on memory usage, since typically thousands
    # of GFF records are created at a time.
    __slots__ = ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
                 'phase', 'attributes')

    def __init__(self, seqid, source, feature_type, start, end,
                 score, strand, phase, attributes=None):

        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes or {}

    @property
    def _key(self):
        return (self.seqid, self.source, self.type, self.start, self.end, 
                self.score, self.strand, self.phase, self.attributes)

    def __eq__(self, other):
        return self._key == other._key

    def __str__(self):
        """Return a GFF formatted string.

        If you need custom attributes ordering, see Formatter.format_attributes().
        """
        return self.formatter.format(self)

    # TODO full repr
    def __repr__(self):
        tpl = 'GFF({}, {}, {}, {})'
        return tpl.format(self.seqid, self.type, self.start, self.end)

    @property
    def ID(self):
        """A helper for accessing the "ID" attribute, which is very common."""
        return self.attributes.get('ID', '')

    @ID.setter
    def ID(self, value):
        """A helper for setting the "ID" attribute, which is very common."""
        self.attributes['ID'] = value

    @property
    def parent_IDs(self):
        """
        A helper for accessing the "Parent" attribute, which is very common.
        Returns a list of parent IDs.

        Note, a record may have multiple parents, so this property
        always returns a list. If you want to ensure only a single parent,
        see GFF.parent_ID (singular).
        """
        parent_IDs = self.attributes.get('Parent')
        if parent_IDs:
            return parent_IDs.split(',')

    @parent_IDs.setter
    def parent_IDs(self, value):
        """
        Given a list of parent IDs, join them with a comma
        and set the "Parent" attribute.
        """
        self.parent_ID = ','.join(value)

    @property
    def parent_ID(self):
        """
        A helper for access a _single_ parent ID.
        If this record has multiple parents, raise a MultipleParents exception.

        Note, GFF v3 allows records to have multiple parents, but having a
        single parent is much more common and useful, so this helper provides
        easy access to a single parent. If you use this on a record with
        parents, we raise an exception so that you know you're probably
        missing something you didn't expect.
        """
        parent_IDs = self.parent_IDs

        if parent_IDs:
            if len(parent_IDs) > 1:
                # Note: this library doesn't yet know how to handle multiple parents.
                #       We raise an exception to ensure the user knows that.
                raise self.MultipleParents()
            else:
                return parent_IDs[0]

    @parent_ID.setter
    def parent_ID(self, value):
        """Given a parent ID, set the "Parent" attribute."""
        self.attributes['Parent'] = value
