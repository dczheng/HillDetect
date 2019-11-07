#!/bin/bash

#make clean
make 

./HillDetect hill.in
./plot.py hill.in 15 
./gen_csv_catalog.py hill.in
