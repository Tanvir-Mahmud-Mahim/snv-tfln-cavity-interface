#!/bin/bash
# FDTD runs of the loaded one-sided cavity used in the revised article,
# in order (see BUILD_MEEP.md).  Arguments of cavity_fdtd.py:
#   N_out  resolution(px/um)  ring-down time after the source (Meep units)
# N_out = 0 runs until the field has decayed (LDOS check, V17); the other
# runs need only the ring-down for Q, the mode volume and the fractions.
# Set PYTHONPATH and LD_LIBRARY_PATH to your Meep installation first.
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
mkdir -p results
for job in "bulk 40" "0 40 1000" "1 40 300" "2 40 300" "3 40 300" "0 50 300"; do
  python cavity_fdtd.py $job > results/fdtd_$(echo $job | tr " " "_").log 2>&1
  echo "finished $job rc=$? $(date)" >> results/fdtd_batch.log
done
