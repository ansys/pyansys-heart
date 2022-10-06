#!/usr/bin/env sh
# export LSTC_FILE=ansys
# export LSTC_LICENSE=ansys
# export ANSYSLMD_LICENSE_FILE=1055@10.110.4.5
# export ANSYSLI_SERVERS=2325@10.110.4.5
echo xxx
mpirun -np 4 mppdyna_d_sse2_linux86_64_intelmmpi i=main.k case