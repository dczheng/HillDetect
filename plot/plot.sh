#!/bin/bash

for i in `cat ./filelist.txt | head 100`
    do
        ./plot_lines.py $i  ./data_fits ./outputs & 
    done
