# Load VMD
module load vmd/1.9.3-nsc1

# Filename of the trajectory file from which frames are extracted
fileName="6-mdnopr-0.002"

# Create the VMD script to extract frames of interest
scriptName="stride.tcl"

cat << __SCRIPT__ > $scriptName

package require topotools

# Load a molecule
mol new $fileName.xtc first 0 last 800 step 50 waitfor all
mol addfile $fileName.gro waitfor all

# Save frames as lammps datafile
animate goto start
animate write gro strides.gro

exit


__SCRIPT__



# Launch VMD script to stride frames
vmd -dispdev text -e $scriptName

# Put results in a new folder
mkdir XRD
mv strides.gro XRD
