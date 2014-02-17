import sys
import ReadInputs as RI

def usage():
    sys.stderr.write('Usage:' + '\n' + 'python ReadInputs.py Experiment_Outline_File Gene_Count_File' + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    else:
        FileIn = RI.ReadInputs(sys.argv[1], sys.argv[2])
        experRF, experRNA = FileIn.ParseExper()
        Genes = FileIn.ReadCount(experRF, experRNA)
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % Genes.geneIDs.size
