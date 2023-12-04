#!/bin/bash
"""

(env38) PS C:\Users\xuhu\pyheart-lib\examples\simulator> python doc_ECG_scenario.py 0.25

GIT BASH:
xuhu@LYOTRAINEE14 MINGW64 ~/pyheart-lib/examples/simulator (438-AI-modeling-for-ECG)
$ bash ./run_ECG_scenario.sh
"""
params=(2 1.5 1 0.5 0.25)

for param in "${params[@]}"
do
    python doc_ECG_scenario.py "$param"
done
