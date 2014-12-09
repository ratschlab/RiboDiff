#!/bin/bash

set -e

#Before Dec4
#python nbinom_count_generator.py --numEntry 10000  --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G10K.Diff0K.cnt.txt --pParamNB 0.0008 --beta1 0.1 --beta2 0.0001
#python nbinom_count_generator.py --numEntry 10000  --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Ribo.Rep3.G10K3.Diff1K.Sh1.5.Sc0.5.cnt.txt --pParamNB 0.008 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --shapeGamma 1.5 --scaleGamma 0.5


#python nbinom_count_generator.py --numEntry 10000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G10K.Diff0K.cnt.txt --nParamNB 1.0 --pParamNB 0.0002 --beta1 0.1 --beta2 0.0001
#python nbinom_count_generator.py --numEntry 10000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Ribo.Rep3.G10K.Diff5K.Sh1.5.Sc0.5.cnt2.txt --nParamNB 1.0 --pParamNB 0.008 --beta1 0.1 --beta2 0.0001 --numDiff 5000 --shapeGamma 1.5 --scaleGamma 0.5

#----------------

#python nbinom_count_generator.py --numEntry 2000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G2K.Diff0K.cnt.txt --nParamNB 1.0 --pParamNB 0.0002 --beta1 0.1 --beta2 0.0001
#python nbinom_count_generator.py --numEntry 2000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Ribo.Rep3.G2K.Diff1K.Sh1.5.Sc0.5.cnt.txt --nParamNB 1.0 --pParamNB 0.008 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --shapeGamma 1.5 --scaleGamma 0.5

#python nbinom_count_generator.py --numEntry 2000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G2K.Diff1K.Sh1.5.Sc0.5.cnt.txt --nParamNB 1.0 --pParamNB 0.0002 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --shapeGamma 1.5 --scaleGamma 0.5
#python nbinom_count_generator.py --numEntry 2000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Ribo.Rep3.G2K.Diff0K.cnt.txt --nParamNB 1.0 --pParamNB 0.008 --beta1 0.1 --beta2 0.0001

#----------------

#python nbinom_count_generator.py --numEntry 2000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Ribo.Rep3.G2K.Diff1K.Sh0.8.Sc0.5.cnt.txt --nParamNB 1.0 --pParamNB 0.008 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --shapeGamma 0.8 --scaleGamma 0.5
python nbinom_count_generator.py --numEntry 2000 --numSampleConA 3 --numSampleConB 3 --output ../exp/Sim/Sim.Rna.Rep3.G2K.Diff1K.Sh0.6.Sc0.5.cnt.txt --nParamNB 1.0 --pParamNB 0.0002 --beta1 0.1 --beta2 0.0001 --diffFile ../exp/Sim/Sim.diffFile --shapeGamma 0.6 --scaleGamma 0.5
