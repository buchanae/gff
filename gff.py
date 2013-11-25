from collections import defaultdict


__version__ = '2.0.0'


class MultipeParents(Exception): pass

def parse_attributes_string(raw):
    attrs = {}
    for token in raw.split(';'):
        token = token.strip()
        if token != '':
            i = token.find('=')
            k = token[:i]
            v = token[i + 1:]
            attrs[k] = v
    return attrs


class GFF(object):

    '''TODO'''

    # TODO measure savings from slots
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

    class ParseError(Exception): pass

    @classmethod
    def from_string(cls, raw):
        '''Parse a GFF3 line.'''

        cols = raw.split('\t')

        # replace columns containing only '.' with None values
        cols = (col if col != '.' else None for col in cols)

        try:
            ref, source, ftype, start, end, score, strand, phase, attrs = cols
        except ValueError:
            raise GFF.ParseError('Invalid number of columns in raw GFF string')

        try:
            if start:
                start = int(start)
        except ValueError:
            raise GFF.ParseError("Couldn't parse start column as an integer")

        try:
            if end:
                end = int(end)
        except ValueError:
            raise GFF.ParseError("Couldn't parse end column as an integer")

        try:
            if score:
                score = float(score)
        except ValueError:
            raise GFF.ParseError("Couldn't parse score column as a float")

        try:
            if phase:
                phase = int(phase)
        except ValueError:
            raise GFF.ParseError("Couldn't parse phase column as a int")

        if attrs:
            attrs = parse_attributes_string(attrs)

        return cls(ref, source, ftype, start, end, score, strand, phase, attrs)

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

    # TODO full repr
    def __repr__(self):
        return 'GFF({}, {}, {}, {})'.format(self.seqid, self.type, self.start, self.end)

    @property
    def ID(self):
        return self.attributes.get('ID', '')

    @property
    def parent_IDs(self):
        parent_IDs = self.attributes.get('Parent')
        if parent_IDs:
            return parent_IDs.split(',')

    @property
    def parent_ID(self):
        parent_IDs = self.parent_IDs

        if parent_IDs:
            if len(parent_IDs) > 1:
                # Note: this library doesn't yet know how to handle multiple parents.
                #       We raise an exception to ensure the user knows that.
                raise MultipleParents()
            else:
                return parent_IDs[0]

    @classmethod
    def from_stream(cls, stream):
        '''Read a GFF3 stream, returning a GFF for every valid line.'''
        for line in stream:
            # skip GFF comment lines
            if line[0] != '#':
                yield cls.from_string(line.strip())



class TreeNode(object):
    __slots__ = ('_parent', 'children')

    def __init__(self):
        self._parent = None
        self.children = []

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


class GFFTreeNode(TreeNode):
    __slots__ = ('_parent', 'children', 'record')

    def __init__(self, record):
        super(GFFTreeNode, self).__init__()
        self.record = record

    def __repr__(self):
        return 'GFFTreeNode({})'.format(repr(self.record))

    @classmethod
    def from_records(cls, records):
        root = cls(None)

        # First we index records by their parent ID.
        by_parent_ID = defaultdict(list)

        # Orphans are records without a parent.
        orphans = []

        for record in records:

            # Some GFF files don't add a 'Parent' attribute for genes,
            # so we fallback to using the seqid.
            parent_IDs = record.parent_IDs or [record.seqid]

            if parent_IDs:
                for parent_ID in parent_IDs:
                    # NOTE: Records with multiple parents are duplicated!
                    node = cls(record)
                    by_parent_ID[parent_ID].append(node)
            else:
                node = cls(record)
                orphans.append(node)

        # Now we make a second pass, linking the nodes and building the tree.
        #
        # We do make two passes because we can't guarantee that
        # a parent is defined before it's children
        # (GFF files are frequently a mess).
        #
        # This two-phase process should be simple and robust,
        # although probably less efficient.
        #
        # TODO some GFF files are huge though, and take up tons of memory.
        #      how could we offer a more efficient (or possibly distributed) version?

        def link_children(node):
            ID = node.record.ID
            for child in by_parent_ID[ID]:
                child.parent = node
                link_children(child)

        for orphan in orphans:
            link_children(orphan)
            orphan.parent = root

        return root
