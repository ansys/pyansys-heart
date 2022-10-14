#!/usr/bin/env sh

echo start lsdyna...
mpirun -np 4 mppdyna_d_sse2_linux86_64_intelmmpi i=main.k case
l2a iter3.binout0000