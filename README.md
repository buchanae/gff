# gff

A parser for [GFF3](http://www.sequenceontology.org/gff3.shtml)-formatted files.

NOTE: This is alpha quality.  It may not strictly implement every feature of the GFF3 format.  See the Known Issues section.

# Usage

    from gff import Feature, Reader

    reader = Reader('/path/to/file.gff')

    for line in reader:
        feature = Feature.from_string(line)

# Known Issues

* Unquoting URL escaped values is not implemented yet.
