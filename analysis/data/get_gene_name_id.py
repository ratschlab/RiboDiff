#!/usr/bin/env python

import sys

def get_gene_name(annotationFile, inputFile, outputFile):

    GeneName = {}

    with open(annotationFile, 'r') as FileIn:
        for line in FileIn:
            if line.lstrip().startswith('#'):
                continue
            lineList = line.rstrip().split('\t')
            if lineList[2] == 'gene':
                gene_id = lineList[8].split('"')[1]
                gene_name = lineList[8].split('"')[3]
                GeneName[gene_id] = gene_name

    FileOut = open(outputFile, 'w')
    with open(inputFile, 'r') as FileIn:
        for line in FileIn:
            if line.strip().split('\t')[0] in GeneName:
                FileOut.write('%s\t%s\n' % (line.strip().split('\t')[0], GeneName[line.strip().split('\t')[0]]))
            else:
                print '%s dose not exists in annotation file.' % line.strip().split('\t')[0]

    FileOut.close()

def get_gene_id(annotationFile, inputFile, outputFile):

    GeneId = {}

    with open(annotationFile, 'r') as FileIn:
        for line in FileIn:
            if line.lstrip().startswith('#'):
                continue
            lineList = line.rstrip().split('\t')
            if lineList[2] == 'gene':
                gene_id = lineList[8].split('"')[1]
                gene_name = lineList[8].split('"')[3]
                GeneId[gene_name] = gene_id

    FileOut = open(outputFile, 'w')
    with open(inputFile, 'r') as FileIn:
        for line in FileIn:
            if line.strip().split('\t')[0] in GeneId:
                FileOut.write('%s\t%s\n' % (GeneId[line.strip().split('\t')[0]], line.strip().split('\t')[0]))
            else:
                print '%s dose not exists in annotation file.' % line.strip().split('\t')[0]

    FileOut.close()

if __name__ == '__main__':

    if sys.argv[4] == 'ForGeneName':
        get_gene_name(sys.argv[1], sys.argv[2], sys.argv[3])
    if sys.argv[4] == 'ForGeneID':
        get_gene_id(sys.argv[1], sys.argv[2], sys.argv[3])
