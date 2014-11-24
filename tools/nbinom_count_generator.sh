#!/bin/bash

set -e

#python nbinom_count_generator.py --numEntry 10000  --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G10K.Diff0K.cnt.txt --pParamNB 0.0008 --beta1 0.1 --beta2 0.0001

#python nbinom_count_generator.py --numEntry 20000  --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G20K.Diff0K.cnt.txt --pParamNB 0.0008 --beta1 0.1 --beta2 0.0001

python nbinom_count_generator.py --numEntry 10000  --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Ribo.Rep3.G10K3.Diff1K.Sh1.5.Sc0.5.cnt.txt --pParamNB 0.008 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --shapeGamma 1.5 --scaleGamma 0.5


