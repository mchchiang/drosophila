#!/bin/bash

num_beads=$1
in_file=$2
out_file=$3

awk -v nbeads=$num_beads '{a=(NR-1)%(nbeads+2)-2;if (a<0){} else if (a<7662){$1="O"} else if (a<14711){$1="N"} else if (a<22893){$1="C"} else if (a<32195){$1="B"} else if (a<32646){$1="F"} else if (a<40121){$1="S"} else if (a==40121){$1="U"} else if (a<nbeads){$1="H"}; print}' $in_file > $out_file 
