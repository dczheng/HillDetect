#!/bin/bash
# run next should indentify

./gen_catalog.py ./test_fits/$1 new_outputs
./gen_skymodel.py ./test_fits/$1 new_outputs
./gen_aomodel_img.py ./test_fits/$1 new_outputs
./plot_final.py ./test_fits/$1 new_outputs
./plot_lset0.py new_outputs $1 $2
./plot_fof.py new_outputs $1 $3
./plot_mask.py test_fits/$1 new_outputs/
