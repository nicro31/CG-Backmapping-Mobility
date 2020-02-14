#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -A snic2018-2-17
#SBATCH -n 8


# Load LAMMPS
module load LAMMPS/22Aug18-nsc1-intel-2018b-eb

declare -a frameName=("0ps.data" "50ps.data" "100ps.data" "150ps.data" "200ps.data" "250ps.data" "300ps.data" "350ps.data" "400ps.data" "450ps.data" "500ps.data" "550ps.data" "600ps.data" "650ps.data" "700ps.data" "750ps.data" "800ps.data")


cd XRD
	
for f in "${frameName[@]}"
do

lmp < $f".in"


done

cd ..

