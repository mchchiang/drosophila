#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

# Read in and set parameters
fraction=$1       # volume fraction of chromatin
run=$2            # trial number
run_dir=$3        # run directory

ploidy=1          # 1 = haploid, 2 = duploid
buffer=5.0        # buffer from each side of the box

# Format variables
fraction=$(python -c "print '%.2f'%($fraction)")
fraction2=$(python -c "print $fraction/20.0")

# Init config (default is random walk)
gen_chromo_exe="../bin/exe/Gen_Drosophila"

# Data files (chromosome lengths and chromatin states)
chromo_file="../../data/drosophila_chromo_length.dat"
hmm_file="../../data/9state_genome_mod.dat"

# Make execution directory
sim_name="dros_vfrac_${fraction}_run_${run}"
run_dir="${run_dir}/${sim_name}"
mkdir $run_dir

# Set output file names
init_file="init_${sim_name}.in"
restart_file="restart_${sim_name}"
prep1_outfile="prep1_${sim_name}.lammpstrj"
prep1_endfile="prep1_${sim_name}.out"
prep2_outfile="prep2_${sim_name}.lammpstrj"
prep2_endfile="prep2_${sim_name}.out"
prep3_outfile="prep3_${sim_name}.lammpstrj"
prep3_endfile="prep3_${sim_name}.out"
pos_file="pos_${sim_name}.dat"
map_file="${sim_name}.lammpsmap"

# Generate chromatin
${gen_chromo_exe} $chromo_file $hmm_file $ploidy $fraction2 $buffer "${run_dir}/${init_file}" "${run_dir}/${map_file}"

# Determine the box size
# Initial box size would be such that the volume fraction is doubled
# that of the main run (obtained from the generated input file)
nbeads=$(awk '{if(NR==3){print $1}}' ${run_dir}/${init_file})
Lhalf=$(awk '{if(NR==11){print $2}}' ${run_dir}/${init_file})
echo "nbeads = $nbeads"
echo "Lhalf = $Lhalf"
L=$(python -c "print $Lhalf*2.0")
r_start=$(python -c "print $Lhalf*(3.0**0.5)") # Initial radius
r_end=$(python -c "print ($nbeads/$fraction)**(1.0/3.0)/2.0") # Final radius
init_box_size=$L
box_size=$(python -c "print $r_end*2.0")

# Initial box size
ilo=$(python -c "print -$init_box_size/2.0")
ihi=$(python -c "print $init_box_size/2.0")

# Final box size
lo=$(python -c "print -$box_size/2.0")
hi=$(python -c "print $box_size/2.0")

max_seed=100000

restart_freq=1000

# Prep/equilibration times
prep1_printfreq=10
prep1_seed=$(python GetRandom.py $max_seed)
prep1_time_1=40 # 4000
prep1_time_2=40 # 4000
prep1_time_3=20 # 2000

prep2_printfreq=10
prep2_seed=$(python GetRandom.py $max_seed)
prep2_time=50 # 5000

prep3_printfreq=10
prep3_seed=$(python GetRandom.py $max_seed)
prep3_time=50 # 5000

lam_atoms=2500
lam_seed=$(python GetRandom.py $max_seed)

delta_t=0.01       # time step size in Brownian time units

# Wall type
wall="bead"        # choose between "ljwall" (default) or "bead"

# Interaction energy parameters
# Harmonic potential
e_harm=100.0

# Soft potential
e_soft=100.0

# LJ potentials
sigma=1.0
cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")

# Normalisation (ensure minimum of potential is actually epsilon)
norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$cutoff)**12-($sigma/$cutoff)**6))")

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
prep1_time_1=$(bc <<< "$prep1_time_1/$delta_t")
prep1_time_2=$(bc <<< "$prep1_time_2/$delta_t")
prep1_time_3=$(bc <<< "$prep1_time_3/$delta_t")
prep1_printfreq=$(bc <<< "$prep1_printfreq/$delta_t")
prep2_time=$(bc <<< "$prep2_time/$delta_t")
prep2_printfreq=$(bc <<< "$prep2_printfreq/$delta_t")
prep3_time=$(bc <<< "$prep3_time/$delta_t")
prep3_printfreq=$(bc <<< "$prep3_printfreq/$delta_t")

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
file="${run_dir}/${lammps_file}"

# Choose template depending on the type of wall used
cp Drosophila_equil.lam $file

# Replace macros in template with input values
sed -i -- "s/INIT_FILE/${init_file}/g" $file
sed -i -- "s/RESTART_FILE/${restart_file}/g" $file

sed -i -- "s/IXLO/${ilo}/g" $file
sed -i -- "s/IXHI/${ihi}/g" $file
sed -i -- "s/IYLO/${ilo}/g" $file
sed -i -- "s/IYHI/${ihi}/g" $file
sed -i -- "s/IZLO/${ilo}/g" $file
sed -i -- "s/IZHI/${ihi}/g" $file

sed -i -- "s/XLO/${lo}/g" $file
sed -i -- "s/XHI/${hi}/g" $file
sed -i -- "s/YLO/${lo}/g" $file
sed -i -- "s/YHI/${hi}/g" $file
sed -i -- "s/ZLO/${lo}/g" $file
sed -i -- "s/ZHI/${hi}/g" $file

sed -i -- "s/RESTART_FREQ/${restart_freq}/g" $file

sed -i -- "s/PREP1_PRINTFREQ/${prep1_printfreq}/g" $file
sed -i -- "s/PREP1_OUTFILE/${prep1_outfile}/g" $file
sed -i -- "s/PREP1_ENDFILE/${prep1_endfile}/g" $file
sed -i -- "s/PREP1_SEED/${prep1_seed}/g" $file
sed -i -- "s/PREP1_TIME_1/${prep1_time_1}/g" $file
sed -i -- "s/PREP1_TIME_2/${prep1_time_2}/g" $file
sed -i -- "s/PREP1_TIME_3/${prep1_time_3}/g" $file

sed -i -- "s/PREP2_PRINTFREQ/${prep2_printfreq}/g" $file
sed -i -- "s/PREP2_OUTFILE/${prep2_outfile}/g" $file
sed -i -- "s/PREP2_ENDFILE/${prep2_endfile}/g" $file
sed -i -- "s/PREP2_SEED/${prep2_seed}/g" $file
sed -i -- "s/PREP2_TIME/${prep2_time}/g" $file

sed -i -- "s/PREP3_PRINTFREQ/${prep3_printfreq}/g" $file
sed -i -- "s/PREP3_OUTFILE/${prep3_outfile}/g" $file
sed -i -- "s/PREP3_ENDFILE/${prep3_endfile}/g" $file
sed -i -- "s/PREP3_SEED/${prep3_seed}/g" $file
sed -i -- "s/PREP3_TIME/${prep3_time}/g" $file

sed -i -- "s/LAM_ATOMS/${lam_atoms}/g" $file
sed -i -- "s/LAM_SEED/${lam_seed}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

sed -i -- "s/EHARM/${e_harm}/g" $file
sed -i -- "s/ESOFT/${e_soft}/g" $file

sed -i -- "s/SIGMA/${sigma}/g" $file

sed -i -- "s/CUTOFF/${cutoff}/g" $file

sed -i -- "s/RSTART/${r_start}/g" $file
sed -i -- "s/REND/${r_end}/g" $file
