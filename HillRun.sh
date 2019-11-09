#!/bin/bash

make clean
make 

./HillDetect hill.in
./plot.py hill.in 15 
#./gen_catalog.py hill.in

./HillDetect hill_mask.in
./plot.py hill_mask.in 15 
#./gen_catalog.py hill_mask.in

add_hills.py hill.in 

./HillDetect hill_add.in
./plot.py hill_add.in fits/test.fits
#./gen_catalog.py hill_add.in
