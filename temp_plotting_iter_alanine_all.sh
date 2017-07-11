#!/usr/bin/env bash

for item in {1..10}; do
    python temp_plotting_iter_alanine.py ${item} 0
    python temp_plotting_iter_alanine.py ${item} 1
done
