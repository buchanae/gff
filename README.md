# gff

A parser for [GFF3](http://www.sequenceontology.org/gff3.shtml)-formatted files.

# Usage

    from gff import Feature, Reader

    reader = Reader('/path/to/file.gff')

    for line in reader:
        feature = Feature.from_string(line)
