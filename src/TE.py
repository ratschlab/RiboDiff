import sys
import readInputs as ri

def usage():
    sys.stderr.write('Usage:' + '\n' + 'python ReadInputs.py Experiment_Outline_File Gene_Count_File' + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    else:
        print '*'*25
        FileIn = ri.readInputs(sys.argv[1], sys.argv[2])
        experRF, experRNA = FileIn.ParseExper()
        data = FileIn.ReadCount(experRF, experRNA)
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size

        print '*'*25
        data = XXXXX(data)
