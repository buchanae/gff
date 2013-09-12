from collections import defaultdict, OrderedDict
import itertools


__version__ = '1.0.0'


def parse_attributes_string(raw):
    attributes = []
    for token in raw.split(';'):
        token = token.strip()
        if token != '':
            i = token.find('=')
            k = token[:i]
            v = token[i + 1:]
            attributes.append((k, v))
    return attributes


class GFF(object):

    '''TODO'''

    class ParseError(Exception): pass

    @classmethod
    def from_string(cls, raw):
        '''Parse a GFF3 line.'''

        cols = raw.split('\t')

        if len(cols) != 9:
            raise GFF.ParseError('Invalid number of columns in raw GFF string')

        # replace columns containing only '.' with None values
        cols = [col if col != '.' else None for col in cols]

        ref, source, ftype, start, end, score, strand, phase, raw_attributes = cols

        try:
            start = int(start)
        except ValueError:
            raise GFF.ParseError("Couldn't parse start column as an integer")
        except TypeError:
            start = None

        try:
            end = int(end)
        except ValueError:
            raise GFF.ParseError("Couldn't parse end column as an integer")
        except TypeError:
            end = None

        try:
            score = float(score)
        except ValueError:
            raise GFF.ParseError("Couldn't parse score column as a float")
        except TypeError:
            score = None

        try:
            phase = int(phase)
        except ValueError:
            raise GFF.ParseError("Couldn't parse phase column as a int")
        except TypeError:
            phase = None

        attributes = ()
        if raw_attributes:
            attributes = parse_attributes_string(raw_attributes)

        return cls(ref, source, ftype, start, end, score, strand, phase, attributes)


    def __init__(self, seqid, source, feature_type, start, end,
                 score, strand, phase, attributes=()):

        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = OrderedDict(attributes)

    @property
    def _key(self):
        return (self.seqid, self.source, self.type, self.start, self.end, self.score,
                self.strand, self.phase, self.attributes)

    def __eq__(self, other):
        return self._key == other._key

    def __str__(self):

        '''Return a GFF3 feature string.'''

        attributes = ';'.join([k + '=' + v for k, v in self.attributes.items()])

        cols = [self.seqid, self.source, self.type, self.start, self.end,
                self.score, self.strand, self.phase, attributes]

        # If any of the columns are None value, repalce them with '.'
        cols = [col if col is not None else '.' for col in cols]

        # Convert all the columns to strings
        cols = [str(col) for col in cols]

        return '\t'.join(cols)

    def __repr__(self):
        return 'GFF({}, {}, {}, {})'.format(self.seqid, self.type, self.start, self.end)

    @classmethod
    def from_stream(cls, stream):
        '''Read a GFF3 stream, returning a GFF for every valid line.'''
        for line in stream:
            # skip GFF comment lines
            if not line.startswith('#'):
                yield cls.from_string(line.strip())


class Tree(object):

    class Node(object):
        def __init__(self):
            self._parent = None
            self.children = []

        def __str__(self):
            return 'Node()'

        @property
        def parent(self):
            return self._parent

        @parent.setter
        def parent(self, value):
            # If this node already has a parent,
            # delete this node from that parent's children
            if self._parent:
                self._parent.children.remove(self)
            self._parent = value
            self._parent.children.append(self)


    def __init__(self):
        self.roots = []

    def walk(self):
        def fn(node):
            yield node
            for child in node.children:
                for node in fn(child):
                    yield node

        return itertools.chain.from_iterable(fn(root) for root in self.roots)


class GFFTree(Tree):

    class Node(GFF, Tree.Node):

        def __init__(self, record):
            GFF.__init__(self, record.seqid, record.source, record.type,
                         record.start, record.end, record.score, record.strand,
                         record.phase, record.attributes)
            Tree.Node.__init__(self)

        def __repr__(self):
            return 'GFFTree.Node({})'.format(GFF.__repr__(self))


    def __init__(self, records):
        super(GFFTree, self).__init__()

        # First we index records by their parent ID.
        by_parent_ID = defaultdict(list)

        # orphans are records without a parent.
        orphans = []

        for record in records:

            parent_IDs = self.record_parent_IDs(record)
            if parent_IDs:
                for parent_ID in parent_IDs:
                    # Records will multiple parents are duplicated
                    record = self.Node(record)
                    by_parent_ID[parent_ID].append(record)
            else:
                record = self.Node(record)
                orphans.append(record)

        # Now we make a second pass, linking the nodes and building the tree.
        #
        # We do make two passes because we can't guarantee that
        # a parent is defined before it's children.
        # (GFF files a frequently a mess)
        #
        # This two-phase process should be simple and robust,
        # although probably less efficient.
        #
        # TODO some GFF files are huge though, and take up tons of memory.
        #      how could we offer a more efficient (or possibly distributed) version?

        def link_children(node):
            ID = self.record_ID(node)
            for child in by_parent_ID[ID]:
                child.parent = node
                link_children(child)

        for orphan in orphans:
            link_children(orphan)

        self.roots = orphans

    def record_ID(self, record):
        return record.attributes.get('ID')

    def record_parent_IDs(self, record):
        p = record.attributes.get('Parent')

        # Frequently GFF files are a mess, and records like genes
        # don't have a "Parent" attribute, but they clearly belong
        # to the "seqid" column.
        if not p and self.record_ID(record) != record.seqid:
            p = record.seqid

        # Records can have multiple parents, separated by a comma
        if p:
            return p.split(',')
