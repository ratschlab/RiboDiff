#!/bin/bash

set -e

#time python -B ../src/TE.py -e ../exp/Sim/ExpOutlineRep3.txt -c ../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.cnt.txt -o ../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.res0.txt -d 0 -r 1 -p 0

#sleep 300s

time python -B ../src/TE.py -e ../exp/Sim/ExpOutlineRep3.txt -c ../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.cnt.txt -o ../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.res1.txt -d 1 -r 1 -p 0

