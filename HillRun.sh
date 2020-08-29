#!/bin/bash

make clean
make 

#1st run
./HillDetect hill.in
./plot.py hill.in 15 
#./gen_catalog.py hill.in

#1st mask run
#./HillDetect hill_mask.in
#./plot.py hill_mask.in 15 
#./gen_catalog.py hill_mask.in

#1st add hills
#./add_hills.py hill.in 

#1st add run
#./HillDetect hill_add.in
#./plot.py hill_add.in fits/ll2000.fits
#./gen_catalog.py hill_add.in
