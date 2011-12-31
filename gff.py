from csv import DictReader


fieldnames = ('seqid', 'source', 'type', 'start', 'end', 
              'score', 'strand', 'phase', 'attributes')


def parse_attributes(raw):
    """Parse a GFF3 attribute string."""

    att = {}
    for token in raw.split(';'):
        if token != '':
            key, value = token.split("=")
            att[key] = value
    return att

def reader(path):
    """Read a GFF3 file, returning a dictionary for each line."""

    def _reader(path):
        with open(path) as fh:
            for line in fh:
                if line[:2] != '##':
                    yield line

    for row in DictReader(_reader(path), fieldnames, delimiter='\t'):
        for key, value in parse_attributes(row['attributes']).iteritems():
            row[key] = value

        row['start'] = int(row['start'])
        row['end'] = int(row['end'])
        
        yield row
