class ParseError(Exception): pass

class Parser(object):
    """
    Parse strings into GFF instances.

    Usually, you don't need to use Parser directly, you only need to use
    GFF.from_string().

    The "attributes_cls" argument allows you to define the data structure
    used to store attributes. By default, attributes are store in a dictionary.
    If you wanted to maintain the order of the attributes, you might use
    collections.OrderedDict instead:

    >>> from collections import OrderedDict
    >>> parser = Parser(GFF_class, attributes_cls=OrderedDict)
    >>> parser.parse(...)
    """

    def __init__(self, GFF_cls, attributes_cls=dict):
        self.GFF_cls = GFF_cls
        self.attributes_cls = attributes_cls

    def parse(self, raw):
        """Parse a single GFF3 record string."""

        cols = raw.split('\t')

        # replace columns containing only '.' with None values
        cols = (col if col != '.' else None for col in cols)

        try:
            ref, source, ftype, start, end, score, strand, phase, attrs = cols
        except ValueError:
            raise ParseError('Invalid number of columns in raw GFF string')

        try:
            if start:
                start = int(start)
        except ValueError:
            raise ParseError("Couldn't parse start column as an integer")

        try:
            if end:
                end = int(end)
        except ValueError:
            raise ParseError("Couldn't parse end column as an integer")

        try:
            if score:
                score = float(score)
        except ValueError:
            raise ParseError("Couldn't parse score column as a float")

        try:
            if phase:
                phase = int(phase)
        except ValueError:
            raise ParseError("Couldn't parse phase column as a int")

        if attrs:
            attrs = self.parse_attributes(attrs)

        return self.GFF_cls(ref, source, ftype, start, end, score, strand,
                            phase, attrs)

    def parse_many(self, lines):
        """Given GFF lines (e.g. from a file), generate GFF instances."""
        for line in lines:
            # skip GFF comment lines
            if line[0] != '#':
                yield self.GFF_cls.from_string(line.strip())

    def parse_attributes(self, raw):
        """
        A helper for parsing the GFF attributes column.

        >>> attributes_str = 'key_1=value_1;key_2=value_2a,value_2b'
        >>> parse_attributes(attributes_str)
        {'key_1': 'value_1', 'key_2': 'value_2a,value_2b'} 

        Note that the result is a dictionary, so the original ordering is lost.
        If you want to maintain ordering, or need some other kind of custom
        storage, use the "cls" keyword.

        >>> from collections import OrderedDict
        >>> parse_attributes
        """

        # Note that because we're storing attributes in a dictionary,
        # we lose the original ordering. We could use collections.OrderedDict,
        # but this is a huge performance hit.
        attrs = self.attributes_cls()
        for token in raw.split(';'):
            token = token.strip()
            if token != '':
                i = token.find('=')
                k = token[:i]
                v = token[i + 1:]
                attrs[k] = v
        return attrs
