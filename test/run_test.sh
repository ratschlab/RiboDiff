#!/bin/bash

set -e

python ../src/TE.py -e ./exp_outline.txt -c ./cnt_table.txt -o test_res.txt -d 0 -r 1 -p 1
