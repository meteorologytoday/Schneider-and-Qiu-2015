#!/bin/bash

selected_dSST=1.0
RH=0.9


for DTheta in 1 5 10 ; do

    python3 plot_param_space.py  \
        --input-dir output_strong_E0 \
        --selected-dSST $selected_dSST \
        --Ugs 0 5 10 15 20 \
        --wvlens 20 40 60 80 100 120 140 160 180 200 \
        --output phase_diagram-DTheta_${DTheta}.png \
        --RH ${RH} \
        --no-display
