#!/usr/bin/env python

import sys
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

def qq_plot(data1, data2, fileOutName):

    idx = np.logical_and(~np.isnan(data1.pval), ~np.isnan(data2.pval)).nonzero()[0]

    prob1 = -np.log10(np.sort(data1.pval[idx], axis=None))
    prob2 = -np.log10(np.sort(data2.pval[idx], axis=None))

    probX = -np.log10(np.sort(np.random.random_sample(idx.size)))

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(10, 3))
    ax1.hist(data1.pval[idx], np.arange(0, 1, 0.05), histtype='bar', align='right', color='tomato', lw=0.5, edgecolor='white')
    ax2.hist(data2.pval[idx], np.arange(0, 1, 0.05), histtype='bar', align='right', color='blue', lw=0.5, edgecolor='white')
    ax3.scatter(probX, prob1, marker='o', color='tomato', s=2, lw=0, label='single disp')
    ax3.scatter(probX, prob2, marker='o', color='blue', s=2, lw=0, label='two disp')
    ax3.plot([min(probX), max(probX)], [min(probX), max(probX)], color='black', linestyle='-')
    ax3.legend(loc='upper left', prop={'size':8})

    ax1.set_ylim(0, 1300)
    ax2.set_ylim(0, 1300)

    ax1.tick_params(labelsize=9)
    ax2.tick_params(labelsize=9)
    ax3.tick_params(labelsize=9)

    plt.savefig(fileOutName, format='pdf', bbox_inches='tight')

if __name__ == '__main__':

    with open(sys.argv[1], 'rb') as FileIn:
        data1 = pickle.load(FileIn)

    with open(sys.argv[2], 'rb') as FileIn:
        data2 = pickle.load(FileIn)

    fileOutName = sys.argv[3]

    qq_plot(data1, data2, fileOutName)
