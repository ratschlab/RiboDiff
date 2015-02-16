#!/bin/bash

set -e

python ../src/TE.py -e exp_outline.txt -c cnt_table.txt -o test_result.txt -d 0 -r 1 -p 1 -q 0.1

echo 'Test run finished. Please check the results in the current directory.'

