#!/bin/bash


python3 plot_2D_SST_bump.py --input-file hycom_GLBv0.08_572_2017010112_t000.nc --lat-rng 25 50 --lon-rng -180 -140 --cutoff-wvlen 2.0 --output-SSTmap sst_analysis_map.png --output-SSTspec sst_analysis_spec.png --no-display


