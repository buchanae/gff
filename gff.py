class Feature(object):

    @classmethod
    def from_string(cls, raw):
        """Parse a GFF3 line."""

        cols = raw.split('\t')
        return cls(*cols)

    def __init__(self, seqid, source, feature_type, start, end, score, strand,
                 phase, raw_attributes):

        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = int(start)
        self.end = int(end)
        self.length = self.end - self.start + 1
        self.score = score
        self.strand = strand
        self.phase = phase
        self.raw_attributes = raw_attributes

        self.attributes = {}
        for token in raw_attributes.split(';'):
            if token != '':
                k, v = token.split("=")
                self.attributes[k] = v


class Reader(object):

    """Read a GFF3 file, returning a dictionary for each line."""

    def __init__(self, path):
        self.path = path

    def __iter__(self):
        def _reader():
            with open(self.path) as fh:
                for line in fh:
                    if line[:2] != '##':
                        yield line.strip()

        return _reader()
