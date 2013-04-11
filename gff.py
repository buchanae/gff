from collections import OrderedDict


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


    def __init__(self, reference, source, feature_type, start, end,
                 score, strand, phase, attributes=()):

        self.reference = reference
        self.source = source
        self.type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = OrderedDict(attributes)


    def __str__(self):

        '''Return a GFF3 feature string.'''

        attributes = ';'.join([k + '=' + v for k, v in self.attributes.items()])

        cols = [self.reference, self.source, self.type, self.start, self.end,
                self.score, self.strand, self.phase, attributes]

        # If any of the columns are None value, repalce them with '.'
        cols = [col if col is not None else '.' for col in cols]

        # Convert all the columns to strings
        cols = [str(col) for col in cols]

        return '\t'.join(cols)


def Reader(stream):
    '''Read a GFF3 stream, returning a GFF for every valid line.'''
    for line in stream:
        # skip GFF comment lines
        if line[:2] != '##':
            yield GFF.from_string(line.strip())
