#!/bin/bash

set -e

time python -B ../src/TE.py -e ../exp/Sim/ExpOutlineRep3.txt -c ../exp/Sim/Sim.Merged.Rep3.G10K3.Diff1KRb.Sh1.5.Sc0.5.cnt.txt -o ../exp/Sim/Sim.Merged.Rep3.G10K3.Diff1KRb.Sh1.5.Sc0.5.res0.txt -d 0 -r 1 -p 0

sleep 300s

time python -B ../src/TE.py -e ../exp/Sim/ExpOutlineRep3.txt -c ../exp/Sim/Sim.Merged.Rep3.G10K3.Diff1KRb.Sh1.5.Sc0.5.cnt.txt -o ../exp/Sim/Sim.Merged.Rep3.G10K3.Diff1KRb.Sh1.5.Sc0.5.res1.txt -d 1 -r 1 -p 0

