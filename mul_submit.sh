#!/bin/bash
for LATO in 1024 2048 4096 8192 16384 32768 
do
  for NPROC in 1 2 4 8 16 21 32 48
  do
    llsubmit "./jobs/FFT_Parallel_${NPROC}_${LATO}.job"
  done
done
