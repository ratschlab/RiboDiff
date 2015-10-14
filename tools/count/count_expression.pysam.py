import sys
import pdb
import pysam
import time
import re
import array
import cPickle
import os

def parse_options(argv):

    """Parses options from the command line """

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()
    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-a', '--annotation', dest='anno', metavar='FILE', help='annotation file in GTF/GFF3 format', default='-')
    required.add_option('-o', '--outfile', dest='outfile', metavar='FILE', help='outfile to store counts in tab delimited format [stdin]', default='-')
    required.add_option('-A', '--alignment', dest='alignment', metavar='FILE', help='alignment in sam or bam format [stdin - sam]', default='-')
    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-f', '--fields', dest='fields', metavar='STRING', help='annotation fields [exon], comma separated', default='exon')
    optional.add_option('-b', '--bam_force', dest='bam_force', action='store_true', help='force BAM as input even if file ending is different from .bam - does not work for STDIN', default=False)
    optional.add_option('-B', '--best_only', dest='best_only', action='store_true', help='count only the best alignment per read [off]', default=False)
    optional.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbosity', default=False)
    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args()
    
    if len(argv) < 2:
        parser.print_help()
        sys.exit(2)

    return options

def parse_anno_from_gff3(options, contigs):
    """This function reads the gff3 input file and returns the information in an
       internal data structure"""

    anno = dict()
    idx2gene = dict()
    gene2idx = dict()

    if options.verbose:
        print >> sys.stderr, "Parsing annotation from %s ..." % options.anno
    
    ### initial run to get the transcript to gene mapping
    if options.verbose:
        print >> sys.stderr, "... init structure"

    trans2gene = dict()
    for line in open(options.anno, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        if sl[2] in ['mRNA', 'transcript']:
            tags=sl[8].split(';')
            assert(tags[0][:2] == 'ID')
            assert(tags[1][:6] == 'Parent')
            tags[0] = tags[0].replace("transcript", "")
            tags[0] = tags[0].replace("mRNA", "")
            tags[1] = tags[1].replace("gene", "")
            trans2gene[tags[0][3:]] = tags[1][7:]

    ### init genome structure
    for c in contigs:
        if options.verbose:
            print 'reserving memory for chr %s of len %s' % (c, contigs[c])
        anno[c] = array.array('H', '\x00\x00' * (contigs[c] + 1))

    ### init list of  considered GFF fields
    fields = options.fields.split(',')

    ### generate a list of exons with attached gene/transcript information
    ### one list per chromsome
    counter = 1
    gene_counter = 1

    t0 = time.time()
    for line in open(options.anno, 'r'):
        if options.verbose and counter % 10000 == 0:
            print >> sys.stderr, '.',
            if counter % 100000 == 0:
                t1 = time.time() - t0
                print >> sys.stderr, "%i - took %.2f secs" % (counter, t1)
                t0 = time.time()
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        if not sl[2] in fields:
            continue

        tags = sl[8].split(';')
        if sl[2] == 'exon':
            tags[0] = tags[0].replace("exon", "")
            trans_id = tags[0][3:]
			
            gene_id = trans2gene[trans_id]
        else:
            print >> sys.stderr, 'Currently only >exon< is supported'
            sys.exit(1)

        if not gene2idx.has_key(gene_id):
            gene2idx[gene_id] = gene_counter
            idx2gene[gene_counter] = gene_id
            gene_counter += 1

        #if sl[0] == 'M':
        #    sl[0] = 'MT'
        ### store for each position of the transcriptome one gene id
        anno[sl[0]][int(sl[3]):int(sl[4]) + 1] = array.array('H', [gene2idx[gene_id]] * (int(sl[4]) + 1 - int(sl[3])))

    if options.verbose:
        print >> sys.stderr, "... done"
    
    if options.verbose:
        print >> sys.stderr, "Dumping exon array ..."

    ### sparsify and dump annotation
    dump_info = open(options.anno + '.dump.info', 'w')
    for k in anno.keys():
        if options.verbose:
            print >> sys.stderr, "... %s" % k
        out_fn = options.anno + '.' + k + '.dump'
        anno[k].tofile(open(out_fn, 'w'))
        print >> dump_info, "%s\t%s\t%i" % (k, out_fn, len(anno[k]))
    dump_info.close()
        
    if options.verbose:
        print >> sys.stderr, "... pickling gene ID map"

    cPickle.dump(idx2gene, open(options.anno + '.pickle', 'w'))

    if options.verbose:
        print >> sys.stderr, "... done"

    return (anno, idx2gene)

def parse_anno_from_gtf(options, contigs):
    """This function reads the gtf input file and returns the information in an
       internal data structure"""

    anno = dict()
    idx2gene = dict()
    gene2idx = dict()

    if options.verbose:
        print >> sys.stderr, "Parsing annotation from %s ..." % options.anno
    
    ### init genome structure
    for c in contigs:
        if options.verbose:
            print 'reserving memory for chr %s of len %s' % (c, contigs[c])
        anno[c] = array.array('H', '\x00\x00' * (contigs[c] + 1))

    ### init list of  considered GFF fields
    fields = options.fields.split(',')

    ### generate a list of exons with attached gene/transcript information
    ### one list per chromsome
    counter = 1
    gene_counter = 1

    t0 = time.time()
    for line in open(options.anno, 'r'):
        if options.verbose and counter % 10000 == 0:
            print >> sys.stderr, '.',
            if counter % 100000 == 0:
                t1 = time.time() - t0
                print >> sys.stderr, "%i - took %.2f secs" % (counter, t1)
                t0 = time.time()
        counter += 1        

        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        
        if not sl[2] in fields:
            continue
        tags = sl[8].split(';')
        gene_id = tags[0][9:-1]
        #transcript_id = '' # tags[1][16:-1]
        if not gene2idx.has_key(gene_id):
            gene2idx[gene_id] = gene_counter
            idx2gene[gene_counter] = gene_id
            gene_counter += 1

        #if sl[0] == 'chrM':
        #    sl[0] = 'chrM_rCRS'
        ### store for each position of the transcriptome one gene id
        anno[sl[0]][int(sl[3]):int(sl[4]) + 1] = array.array('H', [gene2idx[gene_id]] * (int(sl[4]) + 1 - int(sl[3])))

    if options.verbose:
        print >> sys.stderr, "... done"
    
    if options.verbose:
        print >> sys.stderr, "Dumping exon array ..."

    ### sparsify and dump annotation
    dump_info = open(options.anno + '.dump.info', 'w')
    for k in anno.keys():
        if options.verbose:
            print >> sys.stderr, "... %s" % k
        out_fn = options.anno + '.' + k + '.dump'
        anno[k].tofile(open(out_fn, 'w'))
        print >> dump_info, "%s\t%s\t%i" % (k, out_fn, len(anno[k]))
    dump_info.close()
        
    if options.verbose:
        print >> sys.stderr, "... pickling gene ID map"

    cPickle.dump(idx2gene, open(options.anno + '.pickle', 'w'))

    if options.verbose:
        print >> sys.stderr, "... done"

    return (anno, idx2gene)
    #pdb.set_trace()

def read_header(options, infile):
    """Parses the alignment header and extracts contig information"""

    contigs = dict()
    line = ''
    if options.is_bam:
        #chrm = infile.getrname(line.tid).replace('chr', '')
        for i in range(len(infile.references)):
            if contigs.has_key(infile.references[i]):
                if not contigs[infile.references[i]] == infile.lengths[i]:
                    print >> sys.stderr, "Headers in BAM files have inconsistent contig lengths. Stopping ..."
                    sys.exit(1)
            else:
                contigs[infile.references[i]] = infile.lengths[i]

    else:
        for line in infile:
            if not line[0] == '@':
                if len(contigs) == 0:
                    print >> sys.stderr, "No header found in %s. Stopping." % file
                    sys.exit(1)
                else:
                    break

            sl = line.strip().split('\t')

            if not sl[0] == '@SQ':
                continue

            if contigs.has_key(sl[1][3:]):
                if not contigs[sl[1][3:]] == int(sl[2][3:]):
                    print >> sys.stderr, "Headers in BAM files have inconsistent contig lengths. Stopping ..."
                    sys.exit(1)
            else:
                contigs[sl[1][3:]] = int(sl[2][3:])
                        
    return (contigs, line)

def compress_counts(count_list, genes):
    """Takes a list of gene IDs and compresses them to a list of tuples"""

    a = 0
    g = 0
    compressed_list = []

    print >> sys.stderr, " [compressing gene list] ",

    while g < len(genes):
        while g < len(genes) and (a == len(count_list) or genes[g] < count_list[a]):
            g += 1
        if g < len(genes):
            b = a
            while a < len(count_list) and genes[g] == count_list[a]:
                a += 1
            compressed_list.append([genes[g], a - b])
            g += 1
    return compressed_list
    

def main():
    """Main Program Procedure"""
    
    options = parse_options(sys.argv)
    contigs = dict()

    time_total = time.time()

    ### iterate over alignment file(s)
    for file in options.alignment.split(','):
        options.is_bam = False
        ### open file stream
        if file == '-':
            infile = sys.stdin
        elif (len(file) > 3 and file[-3:] == 'bam') or options.bam_force:
            infile = pysam.Samfile(file, 'rb')
            options.is_bam = True
        else:
            infile = open(file, 'r')

        if options.verbose:
            if options.alignment == '-':
                print >> sys.stderr, "Reading alignment from stdin\n"
            else:
                print >> sys.stderr, "Reading alignment from %s\n" % options.alignment

        ### get contigs from alignment data
        if len(contigs) == 0:
            (contigs, lastline) = read_header(options, infile)
            ### TODO handle lastline (line after header) for SAM input

            ### check if we have a version on disk
            if os.path.isfile(options.anno + '.pickle') and os.path.isfile(options.anno + '.dump.info'):
                if options.verbose:
                    t0 = time.time()
                    print >> sys.stderr, 'Loading annotation from pickle/dump files ...'
                idx2gene = cPickle.load(open(options.anno + '.pickle', 'r'))
                anno = dict()
                info_file = open(options.anno + '.dump.info', 'r')
                for line in info_file:
                    sl = line.strip().split('\t')
                    anno[sl[0]] = array.array('H')
                    anno[sl[0]].fromfile(open(sl[1], 'r'), int(sl[2])) 
                    if options.verbose:
                        t1 = time.time() - t0
                        print >>sys.stderr, "... %s took %i secs" % (sl[0], t1)
                        t0 = time.time()
                info_file.close()
                if options.verbose:
                    t1 = time.time() - t0
                    print >> sys.stderr, "... done - last took %i secs" % t1
            else:
                if options.anno[-4:] == 'gff3':
                    ### read annotation from GFF3
                    (anno, idx2gene) = parse_anno_from_gff3(options, contigs)
                else:
                    ### read annotation from GTF
                    (anno, idx2gene) = parse_anno_from_gtf(options, contigs)

        ### count reads
        counter = 1
        t0 = time.time()
        tmp_count = []
        compressed_counts = []
        genes = sorted(idx2gene.keys())
        for line in infile:
            if counter % 100000 == 0:
                print >> sys.stderr, '.',
                if counter % 1000000 == 0:
                    if len(tmp_count) > 5000000:
                        compressed_counts.extend(compress_counts(sorted(tmp_count), genes))
                        tmp_count = []
                    t1 = time.time() - t0
                    print >> sys.stderr, '%i (last 1000000 took %.2f secs)' % (counter, t1) 
                    t0 = time.time()

            counter += 1

            if options.is_bam:
                #chrm = infile.getrname(line.tid).replace('chr', '')
                chrm = infile.getrname(line.tid)
                pos = line.pos
                broken = False

                if line.is_unmapped:
                    continue

                if options.best_only and line.is_secondary:
                    continue

                for o in line.cigar:
                    if o[0] in [0, 2]:
                        for p in range(o[1]):
                            try:
                                g = anno[chrm][pos + p]
                                if g > 0:
                                    tmp_count.append(g)
                                    break
                            except KeyError:
                                continue
                            except IndexError:
                                if chrm in ['chrM', 'M', 'chrM_rCRS']:
                                    continue
                                else:
                                    print >> sys.stderr, 'ERROR: %i exceeds length of %s' % (pos + p, chrm)
                    if broken:
                        break
                    if not o[0] in [1, 5]:
                        pos += o[1]
                
            else:
                sl = line.strip().split('\t')
                if len(sl) < 9:
                    print >> sys.stderr, "ERROR: invalid SAM line\n%s" % line
                    sys.exit(1)

                (size, op) = (re.split('[^0-9]', sl[5])[:-1], re.split('[0-9]*', sl[5])[1:])
                size = [int(i) for i in size]

                #chrm = sl[2].replace('chr', '')
                chrm = sl[2]
                pos = int(sl[3]) - 1
                broken = False

                ## is unmapped ?
                if (int(sl[1]) & 4) == 4:
                    continue

                ## is secondary ?
                if options.best_only and (int(sl[1]) & 256 == 256):
                    continue

                for o in range(len(op)):
                    if op[o] in ['M', 'D']:
                        for p in range(size[o]):
                            try:
                                g = anno[chrm][pos + p]
                                if g > 0:
                                    tmp_count.append(g)
                                    break
                            except KeyError:
                                continue
                            except IndexError:
                                if chrm in ['chrM', 'M', 'chrM_rCRS']:
                                    continue
                                else:
                                    print >> sys.stderr, 'ERROR: %i exceeds length of %s' % (pos + p, chrm)
                    if broken:
                        break
                    if not op[o] in ['H', 'I']:
                        pos += size[o]
                           
        ### close file stream
        if not file == '-':
            infile.close()

        ### compress remaining counts
        compressed_counts.extend(compress_counts(sorted(tmp_count), genes))
        tmp_count = []

        ### report counts to outfile
        outfile = open(options.outfile, 'w')
        print >> sys.stderr, "Sorting and condensing compressed list ..."
        t0 = time.time()
        compressed_counts = sorted(compressed_counts, key = lambda x: x[0])
        for i in range(1, len(compressed_counts)):
            if compressed_counts[i-1][0] == compressed_counts[i][0]:
                compressed_counts[i][1] += compressed_counts[i-1][1]
                compressed_counts[i-1][1] = -1
        compressed_counts = filter(lambda x: x[1] >= 0, compressed_counts)
        t1 = time.time() - t0
        print >> sys.stderr, "... done. took %.2f secs" % t1

        if options.verbose:
            print >> sys.stderr, "Summarizing gene counts ..."

        a = 0
        g = 0
        ### seek to first position that mapped to gene (0 means not gene found)
        while g < len(genes):
            while g < len(genes) and (a == len(compressed_counts) or genes[g] < compressed_counts[a][0]):
                print >> outfile, '%s\t0' % idx2gene[genes[g]]
                if options.verbose and g % 100 == 0:
                    print >> sys.stderr, "%.2f / 100 percent \r" % (float(g) / len(genes) * 100),
                g += 1
            while a < len(compressed_counts) and g < len(genes) and genes[g] == compressed_counts[a][0]:
                print >> outfile, '%s\t%i' % (idx2gene[genes[g]], compressed_counts[a][1])
                a += 1
                g += 1
                if options.verbose and g % 100 == 0:
                    print >> sys.stderr, "%.2f / 100 percent \r" % (float(g) / len(genes) * 100),

        if options.verbose:
            t1 = time.time() - time_total
            print >> sys.stderr, "\n... done - total run took %i secs." % t1
        outfile.close()

if __name__ == "__main__":
    main()
