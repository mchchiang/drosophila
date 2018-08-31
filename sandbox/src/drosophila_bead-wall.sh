#!/bin/bash

# An input script for setting parameters in the LAMMPS config file

# Read in and set parameters
fraction=$1       # volume fraction of chromatin
run=$2            # trial number
run_dir=$3        # run directory
equil_end_file=$4 # end file from equilibration

# Interaction energy parameters
# LJ potentials
e_hethet=1.5
e_hetlam=2.0
e_eueu=0.0
e_tetlam=2.0 # interaction between the tethering bead and the lamina

# Format variables for volume fraction
fraction=$(python -c "print '%.2f'%($fraction)")

# Make execution directory
sim_name="dros-lam-bead_vfrac_${fraction}_run_${run}"
run_dir="${run_dir}/${sim_name}"

if [ ! -d $run_dir ]; then
mkdir $run_dir
fi

# Set output file names
init_file="init_${sim_name}.in"
restart_file="restart_${sim_name}"
prep0_outfile="prep0_${sim_name}.lammpstrj"
prep0_endfile="prep0_${sim_name}.out"
prep1_outfile="prep1_${sim_name}.lammpstrj"
prep1_endfile="prep1_${sim_name}.out"
prep2_outfile="prep2_${sim_name}.lammpstrj"
prep2_endfile="prep2_${sim_name}.out"
run_outfile="run_${sim_name}.lammpstrj"
run_endfile="run_${sim_name}.out"

# Determine the box size
nbeads=$(awk '{if(NR==3){print $1}}' $equil_end_file)
radius=$(python -c "print ($nbeads/$fraction)**(1.0/3.0)/2.0")
box_size=$(python -c "print ($radius+2.5)*2.0")

# Final box size
lo=$(python -c "print -$box_size/2.0")
hi=$(python -c "print $box_size/2.0")

echo "box radius = $radius"
echo "box range = [$lo,$hi]"

# Modify the end file from the equilibration simulation 
# for incorporating the lamina
echo "Modifying input scripts ... "
file=${run_dir}/${init_file}
cp $equil_end_file $file

# Reset the timestep to zero
sed -i -- "1,1s/timestep = [0-9]*/timestep = 0/" $file

# Change the number of atom types
sed -i -- "4,4s/3 atom types/4 atom types/" $file

# Change the box size
sed -i -- "/xlo xhi/c\\${lo} ${hi} xlo xhi" $file
sed -i -- "/ylo yhi/c\\${lo} ${hi} ylo yhi" $file
sed -i -- "/zlo zhi/c\\${lo} ${hi} zlo zhi" $file

# Add a new atom type
sed -i -- "19i4 1" $file

max_seed=100000

# Determine the number of lamina beads
lam_atoms=20000
lam_seed=$(python GetRandom.py $max_seed)

restart_freq=1000

# Prep and run times
thermo_printfreq=10

prep0_printfreq=10
prep0_seed=$(python GetRandom.py $max_seed)
prep0_time=100

prep1_printfreq=10
prep1_seed=$(python GetRandom.py $max_seed)
prep1_time=100

prep2_printfreq=10
prep2_seed=$(python GetRandom.py $max_seed)
prep2_time=500

run_printfreq=10
run_seed=$(python GetRandom.py $max_seed)
run_time=200

delta_t=0.01       # time step size in Brownian time units

# Normalisation (ensure minimum of potential is actually epsilon)
sigma=1.0
cutoff=$(python -c "print '%.13f' % (1.8*$sigma)")
norm=$(python -c "print '%.13f' % (1.0 + 4.0*(($sigma/$cutoff)**12-($sigma/$cutoff)**6))")

e_hethet_norm=$(python -c "print '%.13f' % ($e_hethet/$norm)")
e_hetlam_norm=$(python -c "print '%.13f' % ($e_hetlam/$norm)")
e_tetlam_norm=$(python -c "print '%.13f' % ($e_tetlam/$norm)")
e_eueu_norm=$(python -c "print '%.13f' % ($e_eueu/$norm)")

# Convert all time values to simulation time (i.e. rescale by delta t)
restart_freq=$(bc <<< "$restart_freq/$delta_t")
thermo_printfreq=$(bc <<< "$thermo_printfreq/$delta_t")
prep0_time=$(bc <<< "$prep0_time/$delta_t")
prep0_printfreq=$(bc <<< "$prep0_printfreq/$delta_t")
prep1_time=$(bc <<< "$prep1_time/$delta_t")
prep1_printfreq=$(bc <<< "$prep1_printfreq/$delta_t")
prep2_time=$(bc <<< "$prep2_time/$delta_t")
prep2_printfreq=$(bc <<< "$prep2_printfreq/$delta_t")
run_time=$(bc <<< "$run_time/$delta_t")
run_printfreq=$(bc <<< "$run_printfreq/$delta_t")

# Create the lammps command file based on template
lammps_file="${sim_name}.lam"
qsub_file="qsub_${sim_name}.sh"
file="${run_dir}/${lammps_file}"
qsub_file="${run_dir}/${qsub_file}"
log_file="${sim_name}.log"

# Copy template files
cp Drosophila_bead-wall.lam $file
cp qsub_script.sh $qsub_file

echo "Creating lammps script ... "
# Replace macros in template with input values
sed -i -- "s/INIT_FILE/${init_file}/g" $file
sed -i -- "s/RESTART_FILE/${restart_file}/g" $file

sed -i -- "s/RADIUS/${radius}/g" $file

sed -i -- "s/XLO/${lo}/g" $file
sed -i -- "s/XHI/${hi}/g" $file
sed -i -- "s/YLO/${lo}/g" $file
sed -i -- "s/YHI/${hi}/g" $file
sed -i -- "s/ZLO/${lo}/g" $file
sed -i -- "s/ZHI/${hi}/g" $file

sed -i -- "s/RESTART_FREQ/${restart_freq}/g" $file

sed -i -- "s/THERMO_PRINTFREQ/${thermo_printfreq}/g" $file

sed -i -- "s/PREP0_PRINTFREQ/${prep1_printfreq}/g" $file
sed -i -- "s/PREP0_OUTFILE/${prep1_outfile}/g" $file
sed -i -- "s/PREP0_ENDFILE/${prep1_endfile}/g" $file
sed -i -- "s/PREP0_SEED/${prep1_seed}/g" $file
sed -i -- "s/PREP0_TIME/${prep1_time}/g" $file

sed -i -- "s/PREP1_PRINTFREQ/${prep1_printfreq}/g" $file
sed -i -- "s/PREP1_OUTFILE/${prep1_outfile}/g" $file
sed -i -- "s/PREP1_ENDFILE/${prep1_endfile}/g" $file
sed -i -- "s/PREP1_SEED/${prep1_seed}/g" $file
sed -i -- "s/PREP1_TIME/${prep1_time}/g" $file

sed -i -- "s/PREP2_PRINTFREQ/${prep2_printfreq}/g" $file
sed -i -- "s/PREP2_OUTFILE/${prep2_outfile}/g" $file
sed -i -- "s/PREP2_ENDFILE/${prep2_endfile}/g" $file
sed -i -- "s/PREP2_SEED/${prep2_seed}/g" $file
sed -i -- "s/PREP2_TIME/${prep2_time}/g" $file

sed -i -- "s/RUN_PRINTFREQ/${run_printfreq}/g" $file
sed -i -- "s/RUN_OUTFILE/${run_outfile}/g" $file
sed -i -- "s/RUN_ENDFILE/${run_endfile}/g" $file
sed -i -- "s/RUN_SEED/${run_seed}/g" $file
sed -i -- "s/RUN_TIME/${run_time}/g" $file

sed -i -- "s/DELTA_T/${delta_t}/g" $file

sed -i -- "s/SIGMA/${sigma}/g" $file
sed -i -- "s/CUTOFF/${cutoff}/g" $file

sed -i -- "s/EHETHET/${e_hethet_norm}/g" $file
sed -i -- "s/EHETLAM/${e_hetlam_norm}/g" $file
sed -i -- "s/ETETLAM/${e_tetlam_norm}/g" $file
sed -i -- "s/EEUEU/${e_eueu_norm}/g" $file

sed -i -- "s/LAM_ATOMS/${lam_atoms}/g" $file
sed -i -- "s/LAM_SEED/${lam_seed}/g" $file

# For qsub file
sed -i -- "s/RUNNAME/${sim_name}/g" $qsub_file
sed -i -- "s/LAMMPSSCRIPT/${lammps_file}/g" $qsub_file
sed -i -- "s/LOGFILE/${log_file}/g" $qsub_file
