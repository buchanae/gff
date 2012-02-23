from collections import OrderedDict


__version__ = '0.1'


class Attributes(OrderedDict):
    @classmethod
    def from_string(cls, raw):
        args = []
        for token in raw.split(';'):
            if token != '':
                k, v = token.strip().split("=")
                sp = v.split(',')
                if len(sp) > 1:
                    v = sp
                args.append((k, v))

        return cls(args)

    def __str__(self):
        return ';'.join([k + '=' + ','.join(self.as_list(k)) for k in self.keys()])

    def as_list(self, k):
        if type(self[k]) == list:
            return self[k]
        else:
            return [self[k]]


class Feature(object):

    """TODO"""

    class ParseError(Exception): pass

    @classmethod
    def from_file(cls, path):
        """Read a GFF3 file, returning a Feature for every valid line."""

        with open(path) as fh:
            for line in fh:
                if line[:2] != '##':
                    try:
                        yield cls.from_string(line.strip())
                    except cls.ParseError:
                        pass

    @classmethod
    def from_string(cls, raw):
        """Parse a GFF3 line."""

        cols = raw.split('\t')

        if len(cols) != 9:
            raise Feature.ParseError("invalid number of columns in raw GFF string")

        return cls(*cols)

    def __init__(self, seqid, source, feature_type, start, end, score, strand,
                 phase, raw_attributes):

        self.seqid = seqid
        self.source = source
        self.type = feature_type

        try:
            self.start = int(start)
            self.end = int(end)
        except ValueError:
            raise Feature.ParseError("couldn't parse start or end value")

        self.score = score
        self.strand = strand
        self.phase = phase

        self.attributes = Attributes.from_string(raw_attributes)

    @property
    def length(self):
        return self.end - self.start + 1

    def __str__(self):
        return '\t'.join([self.seqid, self.source, self.type, str(self.start), 
                          str(self.end), self.score, self.strand, self.phase, 
                          str(self.attributes)])
