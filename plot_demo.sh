#!/bin/bash

./plot_lset0.py new_outputs sub40Kv1.fits  20 
./plot_fof.py new_outputs sub40Kv1.fits 20  

./plot_final.py ./test_fits/sub40Kv1.fits new_outputs
./gen_catalog.py  ./test_fits/sub40Kv1.fits new_outputs
