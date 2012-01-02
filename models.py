import random

import gff


class SuperReader(gff.Reader):
    type_to_raw = {
      'chromosome': ['chromosome'],
      'gene': ['gene', 'pseudogene', 'transposable_element_gene'],
      'transcript': ['miRNA', 'mRNA', 'ncRNA', 'pseudogenic_transcript', 
                     'rRNA', 'snoRNA', 'snRNA', 'tRNA'],
      'transcript part': ['exon', 'intron', 'pseudogenic_exon', 'pseudogenic_intron'],
    }
    raw_to_type = {}
    for t in type_to_raw:
        for r in type_to_raw[t]:
            raw_to_type[r] = t


class Genome(Entity):
    chromosomes = OneToMany('Chromosome')

    @classmethod
    def from_gff(cls, gff_file):
        setup_all()
        create_all()

        r = SuperReader(gff_file)
        genome = cls.get_by(digest=r.digest)

        if not genome:
            genome = cls(digest=r.digest)
            # TODO: this caching is pretty weak.  There's a better way
            chromosomes = {}
            genes = {}
            gff_file.seek(0)

            for row in sorted(r, key=lambda r: r['type']):
                if row['type'] == 'chromosome':
                    c = Chromosome(genome=genome, **row)
                    chromosomes[c.name] = c

                elif row['type'] == 'gene':
                    c = chromosomes[row['chrom']]
                    g = Gene(chromosome=c, **row)
                    genes[g.name] = g

                elif row['type'] == 'transcript':
                    g = genes[row['Parent']]
                    Transcript(gene=g, **row)

            session.commit()
            
            genome.delete_overlapping_genes(READ_LENGTH)
            genome.delete_genes_by_min_length(READ_LENGTH * 4)

        return genome

    @property
    def transcripts(self):
        def gen():
            for c in self.chromosomes:
                for g in c.genes:
                    for t in g.transcripts:
                        yield t
        return gen()

    @property
    def genes(self):
        def gen():
            for c in self.chromosomes:
                for g in c.genes:
                    yield g
        return gen()

    def delete_overlapping_genes(self, distance):
        session.autoflush = False
        for c in self.chromosomes:
            s = sorted(c.genes, key=lambda g: g.start)
            for i in xrange(len(s)):
                cur = s[i]
                i -= 1
                while i >= 0:
                    if cur.start - s[i].end < distance:
                        session.delete(cur)
                        session.delete(s[i])
                    i -= 1
        session.commit()
        session.autoflush = True

    def delete_genes_by_min_length(self, length):
        for c in self.chromosomes:
            for g in c.genes:
                if g.length < length:
                    g.delete()
        session.commit()

    def random_gene(self):
        c = random.choice(self.chromosomes)
        return random.choice(c.genes)

    def random_transcript(self):
        gene = self.random_gene()
        #TODO if there are no transcripts, this is an infinite loop
        while True:
            try:
                return gene.random_transcript()
            except NotEnoughTranscriptsError:
                continue


class Chromosome(Feature):
    genome = ManyToOne('Genome')
    genes = OneToMany('Gene', cascade='all,delete-orphan')

    @property
    def sequence(self):
        seqs = []
        offset = 12

        for gene in self.genes:
            pos = gene.start - offset
            if pos < 1:
                pos = 1

            seqs.append((pos, gene.head))
            pos = pos + READ_LENGTH + PAD_LENGTH
            seqs.append((pos, gene.first))
            pos = pos + READ_LENGTH + PAD_LENGTH

            for read in gene.reads:
                seqs.append((pos, read))
                pos += len(read) + PAD_LENGTH

            seqs.append((gene.end - offset - PAD_LENGTH - READ_LENGTH, gene.last))
            seqs.append((gene.end - offset, gene.tail))

        return build_seq(seqs, self.length)


class HitFeature(Feature):
    def __init__(self, *args, **kwargs):
        super(Feature, self).__init__(*args, **kwargs)
        self.init_on_load()

    @orm.reconstructor
    def init_on_load(self):
        self.reads = []
        self.first = seq_gen.next()
        self.last = seq_gen.next()

    def has_room(self, num):
        # does this feature has room for more reads?
        max_len = (self.length - READ_LENGTH * 4) / (READ_LENGTH + PAD_LENGTH)
        return len(self.reads) + num <= max_len

    def add_reads(self, *args):
        if not self.has_room(len(args)):
            raise MaxReads()
        self.reads.extend(args)

    def hit(self):
        r = seq_gen.next()
        self.add_reads(r)
        return r
    

class Gene(HitFeature):
    chromosome = ManyToOne('Chromosome')
    transcripts = OneToMany('Transcript', cascade='all,delete-orphan')

    def __init__(self, *args, **kwargs):
        super(Gene, self).__init__(*args, **kwargs)
        self.init_on_load()

    @orm.reconstructor
    def init_on_load(self):
        super(Gene, self).init_on_load()
        self.head = seq_gen.next()
        self.tail = seq_gen.next()

    def random_transcript(self):
        if len(self.transcripts) == 0:
            raise NotEnoughTranscriptsError()
        return random.choice(self.transcripts)


class Transcript(HitFeature):
    gene = ManyToOne('Gene')
