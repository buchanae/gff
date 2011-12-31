fieldnames = ('seqid', 'source', 'type', 'start', 'end', 
              'score', 'strand', 'phase', 'attributes')


def parse_attributes(raw):
    """Parse a GFF3 attributes string."""

    att = {}
    for token in raw.split(';'):
        if token != '':
            key, value = token.split("=")
            att[key] = value
    return att

def parse(raw):
    """Parse a GFF3 line."""

    d = {}

    for i, val in enumerate(raw.strip().split('\t')):
        d[fieldnames[i]] = val

    for k, v in parse_attributes(d['attributes']).iteritems():
        d[k] = v

    d['start'] = int(d['start'])
    d['end'] = int(d['end'])

    return d


class Reader(object):

    """Read a GFF3 file, returning a dictionary for each line."""

    def __init__(self, path):
        self.path = path


    def __iter__(self):
        def _reader():
            with open(self.path) as fh:
                for line in fh:
                    if line[:2] != '##':
                        yield parse(line)

        return _reader()
