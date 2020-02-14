
package require topotools

# Load a molecule
mol new 6-mdnopr-0.002.xtc first 0 last 800 step 50 waitfor all
mol addfile 6-mdnopr-0.002.gro waitfor all

# Save frames as lammps datafile
animate goto start
animate write gro strides.gro

exit


