#!/usr/bin/env python

import sys
import cPickle as pickle
import gzip
import string

def read_genome(genomeFileIn):

    FileIn = open(genomeFileIn, 'r')

    genomeSeq = {}
    id = ''
    seq = ''
    n = 0
    for line in FileIn:
        sys.stdout.flush()
        if line.lstrip().startswith('>'):
            n += 1
            print '\rStart reading chromosome: %i...' % n,
            if id and seq:
                genomeSeq[id] = seq
            id = line.strip()[1:]
            seq = ''
        else:
            seq = seq + line.strip()

    FileIn.close()

    print '\nRead genome sequence done.'

    return genomeSeq

def find_len_extreme(Obj, seqFeature):

    """
    If there are two or more longest/shortest sequences, this function only returns the first one.
    The "shortest sequences" means the shortest one that longer than 0 bp. If all the isoforms are 
    0 bp, it returns idx as [0].
    """

    objLength = []
    extremeIdx = []
    for i in range(len(Obj)):
        if Obj[i].shape[0] == 0:
            eachLength = 0
        else:
            eachLength = sum(Obj[i][:,1] - Obj[i][:,0] + 1)
        objLength.extend([eachLength])

    if seqFeature == 'LONGEST':
        extremeIdx.append(objLength.index(max(objLength)))
    elif seqFeature == 'SHORTEST':
        if min(objLength) != 0:
            extremeIdx.append(objLength.index(min(objLength)))
        else:
            if sum(objLength) > 0:
                extremeIdx.append(objLength.index(sorted(each for each in objLength if each > 0)[0]))
            else:
                extremeIdx.append(objLength.index(0))
    else:
        pass

    return extremeIdx

def get_sequences(genome, annoFile, geneList, seqType, seqFeature, outputFasta):

    FileIn = gzip.GzipFile(annoFile, 'rb')
    (gene_id, gene_name, gene_biotype, chr, gene_region, gene_strand, transc_id, transc_type, transc, cds, utr5, utr3) = pickle.load(FileIn)
    FileIn.close()
    print 'Load annotation done.'

    FileIn = open(geneList, 'r')
    FileOut = open(outputFasta, 'w')
    for line in FileIn:
        if not line.lstrip().startswith('ENSG'):
            continue
        geneID = line.strip().split('\t')[0]

        if geneID not in chr:
            print '%s dose not exists in annotation file.' % geneID
            #sys.exit()
            continue

        Chr = chr[geneID]
        Strand = gene_strand[geneID]

        if seqType != 'GENE':
            if seqType == 'UTR5':
                Obj = utr5[geneID]
            elif seqType == 'UTR3':
                Obj = utr3[geneID]
            elif seqType == 'CDS':
                Obj = cds[geneID]
            elif seqType == 'TRANSCRIPT':
                Obj = transc[geneID]
            else:
                pass

            if seqFeature == 'LONGEST' or seqFeature == 'SHORTEST':
                idx = find_len_extreme(Obj, seqFeature)
            elif seqFeature == 'ALL':
                idx = range(len(Obj))
            else:
                pass

            for i in idx:
                transcID = transc_id[geneID][i]
                seq = ''
                if Obj[i].shape[0] == 0:
                    continue
                for j in range(Obj[i].shape[0]):
                    left = Obj[i][j, 0]
                    right = Obj[i][j, 1]
                    seqExon = genome[Chr][left-1:right]
                    seq = seq + seqExon

                if Strand == '-':
                    seq = seq[::-1].translate(string.maketrans('ATGCNatgcn', 'TACGNtacgn'))

                FileOut.write('>%s|%s | %s\n' % (geneID, transcID, seqType))
                #FileOut.write('>%s|%s|%s\n' % (geneID, transcID, gene_name[geneID]))
                FileOut.write('%s\n' % seq)

        if seqType == 'GENE':
            left = gene_region[geneID][0]
            right = gene_region[geneID][1]
            seq = genome[Chr][left-1:right]

            if Strand == '-':
                seq = seq[::-1].translate(string.maketrans('ATGCNatgcn', 'TACGNtacgn'))

            FileOut.write('>%s\n' % geneID)
            FileOut.write('%s\n' % seq)

    FileIn.close()
    FileOut.close()

    print 'Get sequences done!'

if __name__ == '__main__':

    if len(sys.argv) == 7:
        genomeFile = sys.argv[1]
        annoPklz = sys.argv[2]
        geneList = sys.argv[3]
        seqType = sys.argv[4]       # 'UTR5', 'UTR3', 'CDS', 'TRANSCRIPT', 'GENE' #
        seqFeature = sys.argv[5]    # 'LONGEST', 'SHORTEST', or 'ALL' #
        outputFasta = sys.argv[6]
    elif len(sys.argv) == 1:
        genomeFile = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.STAR.fa'
        annoPklz = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/annotation/Homo_sapiens.GRCh37.75.primary_assembly.anno.pklz'
        geneList = '/cbio/grlab/projects/RibosomeFootprint/GQ/Lists/AllProtGenes.txt'
        seqType = 'UTR5'
        seqFeature = 'ALL'
        outputFasta = '/cbio/grlab/projects/RibosomeFootprint/GQ/Seq/AllProtGenes.all.fasta'
    else:
        sys.stderr.write('\nmissing input or output arguments.\n\n')

    genome = read_genome(genomeFile)

    get_sequences(genome, annoPklz, geneList, seqType, seqFeature, outputFasta)

