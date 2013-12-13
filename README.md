# gff

Tools for working with [GFF3](http://www.sequenceontology.org/gff3.shtml)-formatted files.

# Usage

    from gff import GFF

    fh = open('path/to/annotation.gff')
    records = GFF.from_file(fh)

    for record in records:
        print record.ID

# Known Issues

* Unquoting URL escaped values is not implemented yet.
