class Formatter(object):

    def format(self, record):
        """Given an instance of GFF, return a GFF-formatted string."""

        cols = [record.seqid, record.source, record.type, record.start,
                record.end, record.score, record.strand, record.phase,
                self.format_attributes(record)]

        # If any of the columns are None value, repalce them with '.'
        cols = [col if col is not None else '.' for col in cols]

        # Convert all the columns to strings
        cols = [str(col) for col in cols]

        return '\t'.join(cols)

    def format_attributes(self, record):
        """
        Given an instance of GFF,
        return it's attributes as a GFF-formatted string.

        GFF.attributes is a dictionary, so we've lost the original ordering,
        (see note in parser.Parser.parse_attributes).
        Here we sort the attributes so that there can be some sort of order.

        If you need custom ordering, override this method.
        """
        attributes = sorted(record.attributes.items())
        attributes = [k + '=' + v for k, v in attributes]
        return ';'.join(attributes)
