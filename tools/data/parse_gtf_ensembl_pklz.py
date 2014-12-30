#!/usr/bin/env python

import sys
import numpy as np
import cPickle as pickle 
import gzip

def get_utr(transcript, codingsqc, strand):

    """ Assume protein coding exon includes stop codon AND the coordinate of transcript / cds exons are sorted from small to large """

    if strand == '+':
        if codingsqc[0, 0] > transcript[0, 0]:
            idxUTR5 = (codingsqc[0, 0] > transcript[:, 0]).nonzero()[0]
            UTR5 = transcript[0:idxUTR5[-1]+1, :].copy()
            if codingsqc[0, 0] > transcript[idxUTR5[-1], 0] and codingsqc[0, 0] <= transcript[idxUTR5[-1], 1]:
                UTR5[-1, 1] = codingsqc[0, 0] - 1
        else:
            UTR5 = np.array([])

        if codingsqc[-1, 1] < transcript[-1, 1]:
            idxUTR3 = (codingsqc[-1, 1] < transcript[:, 1]).nonzero()[0]
            UTR3 = transcript[idxUTR3[0]:, :].copy()
            if codingsqc[-1, 1] < transcript[idxUTR3[0], 1] and codingsqc[-1, 1] >= transcript[idxUTR3[0], 0]:
                UTR3[0, 0] = codingsqc[-1, 1] + 1
        else:
            UTR3 = np.array([])
    elif strand == '-':
        if codingsqc[-1, 1] < transcript[-1, 1]:
            idxUTR5 = (codingsqc[-1, 1] < transcript[:, 1]).nonzero()[0]
            UTR5 = transcript[idxUTR5[0]:, :].copy()
            if codingsqc[-1, 1] < transcript[idxUTR5[0], 1] and codingsqc[-1, 1] >= transcript[idxUTR5[0], 0]:
                UTR5[0, 0] = codingsqc[-1, 1] + 1
        else:
            UTR5 = np.array([])

        if codingsqc[0, 0] > transcript[0, 0]:
            idxUTR3 = (codingsqc[0, 0] > transcript[:, 0]).nonzero()[0]
            UTR3 = transcript[0:idxUTR3[-1]+1, :].copy()
            if codingsqc[0, 0] > transcript[idxUTR3[-1], 0] and codingsqc[0, 0] <= transcript[idxUTR3[-1], 1]:
                UTR3[-1, 1] = codingsqc[0, 0] - 1
        else:
            UTR3 = np.array([])
    else:
        pass

    return (UTR5, UTR3)

def parse_gtf(annotationFile, pklzFile):

    gene_id = []
    gene_name = {}
    gene_biotype = {}
    chr = {}
    gene_region = {}
    gene_strand = {}
    transc_id = {}
    transc_type = {}
    transc = {}
    cds = {}
    utr5 = {}
    utr3 = {}

    transcript = {}
    codingsqc = {}
    proteinstop = {}

    with open(annotationFile, 'r') as FileIn:

        for line in FileIn:
            if line.lstrip().startswith('#'):
                continue

            lineList = line.rstrip().split('\t')

            if lineList[2] == 'gene':
                geneID = lineList[8].split('"')[1]
                gene_id.extend([geneID])

                geneName = lineList[8].split('"')[3]
                gene_name[geneID] = geneName

                geneBioType = lineList[8].split('"')[7]
                gene_biotype[geneID] = geneBioType

                chr[geneID] = lineList[0]

                gene_region[geneID] = [int(lineList[3]), int(lineList[4])]

                gene_strand[geneID] = lineList[6]

            if lineList[2] == 'exon':
                geneID = lineList[8].split('"')[1]
                transcID = lineList[8].split('"')[3]
                transcType = lineList[1]
                if transcID not in transcript:
                    transcript[transcID] = np.array([[int(lineList[3]), int(lineList[4])]])
                    if geneID not in transc_id:
                        transc_id[geneID] = [transcID]
                        transc_type[geneID] = [transcType]
                    else:
                        transc_id[geneID].extend([transcID])
                        transc_type[geneID].extend([transcType])
                else:
                    transcript[transcID] = np.vstack([transcript[transcID], np.array([int(lineList[3]), int(lineList[4])])])

            if lineList[2] == 'CDS':
                transcID = lineList[8].split('"')[3]
                if transcID not in codingsqc:
                    codingsqc[transcID] = np.array([[int(lineList[3]), int(lineList[4])]])
                else:
                    codingsqc[transcID] = np.vstack([codingsqc[transcID], np.array([int(lineList[3]), int(lineList[4])])])

            if lineList[2] == 'stop_codon':
                transcID = lineList[8].split('"')[3]
                if transcID not in proteinstop:
                    proteinstop[transcID] = np.array([[int(lineList[3]), int(lineList[4])]])
                else:
                    sys.stderr.write('More than one stop codon: %s, %s, %s.\n' % (geneID, transcID, transcType))

    for key in transc_id:
        transc[key] = []
        cds[key] = []
        utr5[key] = []
        utr3[key] = []
        for each_transcID in transc_id[key]:
            # transcript exons' coordinate order: small to large
            if gene_strand[key] == '+':
                TRANSC = transcript[each_transcID]
            else:
                TRANSC = np.flipud(transcript[each_transcID])
            transc[key].extend([TRANSC])

            if each_transcID in codingsqc:
                # cds exons' coordinate order: small to large
                if gene_strand[key] == '+':
                    CDS = codingsqc[each_transcID]
                    # cds exon includes the stop codon
                    if each_transcID in proteinstop:
                        CDS[-1, 1] = CDS[-1, 1] + 3
                else:
                    CDS = np.flipud(codingsqc[each_transcID])
                    if each_transcID in proteinstop:
                        CDS[0, 0] = CDS[0, 0] - 3
                cds[key].extend([CDS])

                (UTR5, UTR3) = get_utr(TRANSC, CDS, gene_strand[key])
                utr5[key].extend([UTR5])
                utr3[key].extend([UTR3])
            else:
                cds[key].extend([np.array([])])
                utr5[key].extend([np.array([])])
                utr3[key].extend([np.array([])])

    genes = (gene_id, gene_name, gene_biotype, chr, gene_region, gene_strand, transc_id, transc_type, transc, cds, utr5, utr3)

    FileOut = gzip.GzipFile(pklzFile, 'wb')
    pickle.dump(genes, FileOut, pickle.HIGHEST_PROTOCOL)
    FileOut.close()

if __name__ == '__main__':

    if len(sys.argv) == 1:
        fileInName = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/annotation/Homo_sapiens.GRCh37.75.primary_assembly.gtf'
        fileOutName = '/cbio/grlab/projects/Ivashkiv/MacrophageRpfRnaSeq/results/motif/annotation/Homo_sapiens.GRCh37.75.primary_assembly.anno.pklz'
    elif len(sys.argv) == 3:
        fileInName = sys.argv[1]
        fileOutName = sys.argv[2]
    else:
        sys.stderr.write('\nmissing input GTF file or output pklz file.\n\n')

    parse_gtf(fileInName, fileOutName)

