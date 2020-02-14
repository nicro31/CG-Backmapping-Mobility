#############################################################
#############################################################
#########             LOE-CTP-FRAG                   ########
#########  N.Rolland, LOE Linköping University/KTH   ########
#########      A code for mobility calculation       ########
#############################################################
#############################################################

A C++ code, with dependancies on ApproxMVBB and QC_Tools (source files/compilation included in the sofware installation).
External python scripts are required for a few operations, therefore you need a python3 installation. You also need VMD.
Finally, you need a suitable compiler suite (note LOE-CTP-FRAG has been compiled with intel compiler suite until now).


################
# INSTALLATION #
################
- In the folder where you want to install the software do: 
git clone https://bitbucket.org/nicro31/loe-ctp-frag
- Then usual cmake and make procedure
- In the scripts folder, there is Load_LOE-CTP-FRAG.sh that is useful to load the software when you want to use it. Note intel, python and vmd modules should be adapted to your environment (present commands are suitable for NSC tetralith/gamma)
- Enjoy !


###############
#### RUN ######
###############

0) Take a look in the example folder to follow this quick introduction.

1) Let's assume you have a morphology ready in a gromacs file morphology.gro. You first need to create a .map file (or several in case you have different molecules that participate to the transport) that describes how the transport units are built. This file has the same name than the residues in the .gro file. In this .map file, you define a set of rigid fragments (e.g. [NDI1], [NDI2]...) which contain a list of atom with their gromacs serial number and name (e.g. C,H,O,S...). You can also indicate if you want to add some hydrogens to some atoms in the system (useful to saturate bonds if a transport unit is defined as part of a molecule only). Finally you define a tranport unit (e.g. {TRA}) by indicating the rigid fragments that formed it. Note the way you define rigid fragments/transport units is very important because for each transport unit, a position and orientation will be calculated. From these positions, an average position will be calculated for each transport unit, and two transport units will be linked according to the distance between their two closest rigid fragments. Finally, you can prepare the system with the following command in interactive mode (or even in the logging node but in this case please do: export OMP_NUM_THREADS= a reasonnable number of processors to not disturb other users on the login node, e.g. 4):

LOE-CTP-FRAG prep morphology.gro config.txt morpho.pdb

The .map and config.txt file have to be in the same folder than the morphology.gro file. morpho.pdb is the output file that will be used in the following calculations. Note you have to set the parameters that you want in the config.txt file. At this stage, the most important are probably the distance threshold between rigid fragments to link the transport units, and the number of shared processors to use for the Gaussian calculations. When this step is done, you have the morpho.pdb file that contains the morphology (note the morphology contains 7 replicants of the original simulation box to take into account pbc) and an orientation.tcl file that contains position and orientations of the rigid fragments (do not delete it because it is used in the following calculations!!!). You can visualize position and orientations of the rigid fragments by loading the orientation.tcl file in the vmd Tk console mode. Also, you have now a bunch of folders with names "site...". These contains the .com files to lauch Gaussian calculations. In case you realize that you need to change the distance threshold for instance, note you can run:

LOE-CTP-FRAG prep_gaussian morpho.pdb config.txt

This command will generate the new Gaussian files with increased distance threshold, without re-calculating the orientations/positions of the rigid fragments.

2) Launch the Gaussian calculations. To check the number of calculations that you need to do, and the number of calculations already performed, you can do Check_GProgress.sh. To launch the calculations, you can use the gaussian.py python script provided with LOE-CTP-FRAG. Do gaussian.py -h to get some help about it. An example of the command would be:

gaussian.py 10 32 snic2017-12-59 10:00:00

It means that we launch 10 distinct SLURM jobs, each with 32 processors, with account snic2017-12-59 and time limit 10 hours.

3) When all Gaussian calculations are done (you may need to launch gaussian.py several times in case Check_GProgress tells you that some calculations are still to do), you can calculate the transfer integrals between transport units. In config.txt, set the number of orbitals you want to use (NbrLevels2Read) and the orbital type, i.e. HOMO for hole transport and LUMO for electron transport. Note the number of figures after levelEnergy should match NbrLevels2Read. Intra-molecular transfer integrals can be set to a user-defined fixed value if you want with H_intra (if < 0, results from Gaussian calculations are used), and inter-level tranfer integral is set with H_level. The command to run the transfer integral calculation is:

LOE-CTP-FRAG build morpho.pdb config.txt morpho_homo.grp

The .grp is the output file with transfer integrals. Do not forget to put a different name for this file when you calculate holes transfer integrals and electron transfer integrals.

4) Now you can calculate a mobility !! The command is:

LOE-CTP-FRAG mobility morpho.pdb config.txt morpho_homo.grp mob

The mob is just a folder in which results will be put (Can be useful to distinguish hole and electron transport for instance). In the config.txt file you can change a couple of parameters if you want before running the calculation. For instance, the w0 parameters are the prefactors for the hopping rates for the different kinds of hopping (inter, intra- molecular, inter-levels). The nbrLevels is the number of levels you want to use for each transport unit. Of course, it has to be lower or equal to the number of orbitals calculated previously (i.e. NbrLevels2Read). Also, note you can use either the orbital energies coming from ZINDO calculations (in case DOS broadening < 0) or energies coming from a truncated gaussian random distribution with broadening dosBroadening. In this later case, for each orbitals corresponding to a given level (e.g. HOMO, HOMO-1, HOMO-2...), the gaussian random distribution is centered on the values indicated after levelEnergy. Note the folder where the results are saved will have a specific number; this corresponds to the seed of the gaussian random distribution, and can be used later on to redo the exact same calculation or check how is the density of states from the gaussian random distribution. Finally, you can change, electric field, temperature and charge concentration values, and for each you can indicate several values, e.g.:
///Field values for mobility calculation(V/cm, value in x, y, z directions)///
3
1e4 0 0
0 1e4 0
0 0 1e4
This will calculate mobility in the x, y and z directions

5) A couple of tools to interpret the mobility values:

LOE-CTP-FRAG export_transfer morpho.pdb config.txt morpho_homo.grp transfer.txt -> Output transfer integral distribution in the text file transfer.txt

LOE-CTP-FRAG export_energy morpho.pdb config.txt morpho_homo.grp energy.txt -> Output site energies (i.e. density of states) in the text file energy.txt (if you indicate a seed as a last optionnal argument, you will get DOS from the gaussian random distribution; otherwise, you will get DOS from ZINDO calculation)

LOE-CTP-FRAG perco morpho.pdb config.txt morpho_homo.grp perco.txt -> Output the percolation curve in perco.txt; parameters for this are set in the config.txt file (Hmin,Hmax,dH indicate minimum threshold for percolation, maximum value, and step between to percolation calculation)

LOE-CTP-FRAG export_geometry morpho.pdb config.txt morpho_homo.grp geometry.txt -> Output some link properties between transport units. You get the geometry file geometry.txt with distance between transport units, and relative misorientation.
