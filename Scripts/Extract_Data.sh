module load vmd/1.9.3-nsc1

declare -a frame=("0ps" "50ps" "100ps" "150ps" "200ps" "250ps" "300ps" "350ps" "400ps" "450ps" "500ps" "550ps" "600ps" "650ps" "700ps" "750ps" "800ps")

count=0

cd XRD

for f in "${frame[@]}"
do

# Create the VMD script to extract frames of interest
scriptName="extract.tcl"

cat << __SCRIPT__ > $scriptName

package require topotools


# Load a molecule
mol new strides.gro first $count last $count

# Select non water molecules
atomselect top "not resname SOL"
animate write gro temp.gro sel atomselect0

# Load file without solvent
mol delete 0
mol new temp.gro

# Save frame as lammps datafile
topo writelammpsdata $f.data full

exit

__SCRIPT__

    count=$(($count+1))

vmd -dispdev text -e extract.tcl

done

cd ..

