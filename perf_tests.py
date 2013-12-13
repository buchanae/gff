from __future__ import print_function

from timeit import timeit

from gff import GFF


record = 'Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010'


if __name__ == '__main__':
    setup = 'from __main__ import GFF, record'
    test = 'GFF.from_string(record)'
    print(timeit(test, setup=setup))
