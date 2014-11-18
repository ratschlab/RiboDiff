#!/bin/bash

set -e

#python nbinom_count_generator.py --numEntry 2000  --numSampleConA 3 --numSampleConB 3 --output ./test/disp012/Sim/RiboCnt1.0Fold.Beta1.0.txt --pParamNB 0.01 --beta1 1.0 --beta2 0.01
#python nbinom_count_generator.py --numEntry 2000  --numSampleConA 3 --numSampleConB 3 --output ./test/disp012/Sim/RnaCnt1.0Fold.Beta0.1.txt --pParamNB 0.006 --beta1 0.1 --beta2 0.001
#
#python nbinom_count_generator.py --numEntry 2000  --numSampleConA 3 --numSampleConB 3 --output ./test/disp012/Sim/RiboCnt2.0Fold.Beta0.5.txt --pParamNB 0.01 --beta1 0.5 --beta2 0.001 --numDiff 1000 --meanFoldDiff 2.0
#cat ./test/disp012/Sim/RiboCnt2.0Fold.Beta0.5.txt | tail +2 | cut -f 8 > ./test/disp012/Sim/dispersions.txt
#python nbinom_count_generator.py --numEntry 2000  --numSampleConA 3 --numSampleConB 3 --output ./test/disp012/Sim/RnaCnt1.0Fold.useDisp.txt --pParamNB 0.006 --dispFile ./test/disp012/Sim/dispersions.txt


#python nbinom_count_generator.py --numEntry 20000  --numSampleConA 10 --numSampleConB 10 --output ./test/disp012/Sim/RiboCnt2.0Fold.Beta0.1.Rep10.20K.txt --pParamNB 0.01 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --meanFoldDiff 2.0
#python nbinom_count_generator.py --numEntry 20000  --numSampleConA 10 --numSampleConB 10 --output ./test/disp012/Sim/RnaCnt1.0Fold.Beta0.1.Rep10.20K.txt --pParamNB 0.001 --beta1 0.1 --beta2 0.0001

python nbinom_count_generator.py --numEntry 20000  --numSampleConA 3 --numSampleConB 3 --output ./test/disp012/Sim/RiboCnt.Beta0.1.Rep3.G20K.Diff1K.Sh2.0.Sc0.5.txt --pParamNB 0.01 --beta1 0.1 --beta2 0.0001 --numDiff 1000 --shapeGamma 2.0 --scaleGamma 0.5
python nbinom_count_generator.py --numEntry 20000  --numSampleConA 3 --numSampleConB 3 --output ./test/disp012/Sim/RnaCnt.Beta0.1.Rep3.G20K.Diff0K.txt --pParamNB 0.001 --beta1 0.1 --beta2 0.0001
