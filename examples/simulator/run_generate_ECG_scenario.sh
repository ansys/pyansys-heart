#!/bin/bash
"""
Execute this file in GIT BASH:
xuhu@LYOTRAINEE14 MINGW64 ~/pyheart-lib/examples/simulator
$ bash ./run_ECG_scenario.sh
"""
# params=(2 1.5 1 0.5 0.25)
param=(2 1.5 1 0.5)

for param in "${params[@]}"
do
    python doc_ECG_scenario.py "$param"
done