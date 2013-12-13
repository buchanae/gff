"""
Tools for working with a tree of GFF records.
"""

from collections import defaultdict


class GFFTreeNode(object):

    """
    Represents a node in a tree of GFF records.
    
    A tree node has a single parent, a list of children,
    and a reference to the record it represents.
    """

    __slots__ = ('_parent', 'children', 'record')

    def __init__(self, record):
        self._parent = None
        self.children = []
        self.record = record

    def __repr__(self):
        return 'GFFTreeNode({})'.format(repr(self.record))

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

    def walk(self):
        yield self
        for child in self.children:
            for node in child.walk():
                yield node
            

def build_tree(records, cls=GFFTreeNode):
    """
    Given GFF records, build a tree of instances.
    
    Usually GFF records represent a tree of genomic features:
    an annotation has chromosomes, chromosomes have genes,
    genes have transcripts, transcripts have exons, and so on.
    
    Often it's useful to work with the annotation in this tree form.
    For example, if you have an annotation like this:

    - Chromosome 1
    -- gene A
    --- transcript A.1
    --- transcript A.2
    - Chromosome 2
    ... etc ...

    You might want to do this:

    >>> records = GFF.from_file(open('path/to/annotation.gff'))
    >>> tree = build_tree(records)
    >>> tree.children[0].ID
    ... 'Chromosome 1'
    >>> chr1_genes = tree.children[0].children
    >>> chr1_genes[0].children[0].ID
    ... 'transcript A.1'


    Note, this tree building business can be tricky. GFF files often come
    in different forms: some put Parent attributes on genes, some don't,
    some put ID attributes on chromosomes, some don't. This method
    isn't intended to be robust enough to handle all those varieties.
    """

    root = cls(None)

    by_ID = {}

    def link(node, parent_ID):
        try:
            parent = by_ID[parent_ID]
            if parent is not node:
                node.parent = parent

        except KeyError:
            pass

    # Orphans are records without a parent.
    orphans = []

    for record in records:

        # Some GFF files don't add a 'Parent' attribute for genes,
        # so we fallback to using the seqid.
        parent_IDs = record.parent_IDs or [record.seqid]

        # Important! If a record has multiple parents,
        # multiple nodes are created.
        for parent_ID in parent_IDs:

            node = cls(record)

            ID = record.ID
            if ID:
                by_ID[ID] = node

            link(node, parent_ID)

            if not node.parent:
                orphans.append((node, parent_ID))

    # It's possible that the GFF records were defined out of order,
    # and the parent came after the child, so make a second pass
    # and attempt to link any orphans to their parents.
    #
    # If no parents were found for the orphan, link it to the root node.
    for orphan, parent_ID in orphans:
        link(orphan, parent_ID)
        if not orphan.parent:
            orphan.parent = root

    return root
