#!/bin/bash

vfrac=$1
run_start=$2
run_end=$3
run_inc=$4
dir=$5

vfrac=$(python -c "print '%.2f' % ($vfrac)")
run=$run_start

for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
do
    echo "Creating files for vfrac = $vfrac run = $run"
    bash drosophila_equil.sh $vfrac $run $dir
done

