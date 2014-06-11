import sys
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python plot.py data.pkl FigDispersion.pdf FigTEchange.pdf' + '\n')

def plot_dispersion(data, fileName):

    countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
    countRiboMean=np.reshape(countRiboMean, (len(countRiboMean), 1))
    disperRaw = data.disperRaw
    disperAdj = data.disperAdj
    beta = data.beta

    fig, ax = plt.subplots()
    index = ~np.isnan(disperRaw)
    ax.scatter(countRiboMean[index], disperRaw[index], marker='o', color='#C0C0C0', s=5, lw=0, label='Raw')
    ax.scatter(countRiboMean[index], disperAdj[index], marker='o', color='#6666FF', s=5, lw=0, label='Adjusted')
    countRiboMeanSorted = np.sort(countRiboMean[index])
    ax.plot(countRiboMeanSorted, beta[0]/countRiboMeanSorted + beta[1], color='#FF3333', label='Fitted')

    ax.legend(loc='upper right')

    xLowBound = 0
    xUpBound  = np.percentile(countRiboMean[index], 90)
    yLowBound = np.percentile(disperRaw[index], 5)
    yUpBound  = np.percentile(disperRaw[index], 95)
    ax.set_xlim(xLowBound, xUpBound)
    ax.set_ylim(yLowBound, yUpBound)
    #ax.set_xlim(0, 2000)
    #ax.set_ylim(0, 0.1)

    plt.xlabel('Mean count of Ribosome Footprint')
    plt.ylabel('Dispersion')

    plt.savefig(fileName)

def plot_TEchange(data, fileName, threshold):

    countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
    countRiboMean=np.reshape(countRiboMean, (len(countRiboMean), 1))
    foldChange = data.logFoldChange

    index = ~np.isnan(data.disperRaw)
    countRiboMeanEff = countRiboMean[index]
    foldChangeEff = foldChange[index]
    padjEff = data.padj[index]

    idx = padjEff < threshold
    countRiboMeanSig = countRiboMeanEff[idx]
    foldChangeSig = foldChangeEff[idx]

    fig, ax = plt.subplots()
    ax.scatter(countRiboMeanEff, foldChangeEff, marker='o', color='#C0C0C0', s=3, lw=0, label='Tested genes')
    ax.scatter(countRiboMeanSig, foldChangeSig, marker='o', color='#FF9933', s=3, lw=0, label='Significant genes')
    ax.legend(loc='upper right')
    plt.xlabel('Mean count of Ribosome Footprint')
    plt.ylabel('log TE Fold Change')

    xLowBound = (np.percentile(countRiboMeanSig, 97.5) - np.min(countRiboMeanSig)) * -0.02
    xUpBound  = np.percentile(countRiboMeanSig, 97.5)
    yLowBound = np.percentile(foldChangeSig, 2.5)
    yUpBound  = np.percentile(foldChangeSig, 97.5)

    ax.set_xlim(xLowBound, xUpBound)
    ax.set_ylim(yLowBound, yUpBound)

    plt.savefig(fileName)

if __name__ == '__main__':

    if len(sys.argv) != 4:
        usage()
    else:
        with open(sys.argv[1], 'rb') as FileIn:
            data = pickle.load(FileIn)

        plot_dispersion(data, sys.argv[2])
        plot_TEchange(data, sys.argv[3], threshold=0.05)
