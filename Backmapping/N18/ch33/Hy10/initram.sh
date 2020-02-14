#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -A snic2017-12-59
#SBATCH -n 16

export PATH=$PATH:$(pwd)


PROGRAM=initram.sh
VERSION=0.1
VERSTAG=devel-20130709-21-TAW
AUTHOR="Tsjerk A. Wassenaar, PhD"
YEAR="2013"
AFFILIATION="
University of Calgary
2500 University Drive NW
Calgary, Alberta
Canada, T2N 1N4"

CMD="$0 $@"

DESCRIPTION=$(cat << __DESCRIPTION__
This is a convenience script to convert a coarse-grained system
into a united-atom or all-atom representation by projection
and subsequent relaxation. The script is a wrapper around 
backward.py, which is called for the projection of the coarse-
grained system to unit-atom or all-atom, and GROMACS, which is
used for energy minimization and further relaxation using 
molecular dynamics.

The script requires a COARSE-GRAINED input structure and a 
corresponding ATOMISTIC topology file. The topology file will
be used to define the target particle list, allowing setting 
protonation states, virtual sites, and even mutations by 
providing an alternative topology. Atom list, atom names and 
residue names from the topology file take precedence over the
ones in the structure file.

The script has a number of options to control the backmapping
process. The default options are chosen to give a robust protocol,
which can be used for backmapping of virtually any system. This 
protocol is tested for systems containing protein, membrane and 
solvent, with a total target size of up to a million particles.

The default protocol consists of the following steps:

1. Projection
2. Energy minimization (500 steps), excluding non-bonded interactions 
   between atoms in selected groups, e.g. membrane and protein.
3. Energy minimization (500 steps)
4. Position restrained NVT simulation (300K), with time step 0.2 fs
5. Position restrained NVT simulation (300K), with time step 0.5 fs
6. Position restrained NVT simulation (300K), with time step 1.0 fs
7. Position restrained NVT simulation (300K), with time step 2.0 fs

These steps can be modified using the options, as listed below. 

It is possible to extend the protocol using a number of user-supplied
run parameter (.mdp) files. These will be run after step 7. Multiple
mdp files can be supplied, which will be processed in the order given
(e.g., -mdp param1.mdp -mdp param2.mdp).

At the end of the process timing information is given summarizing the 
time per step and the cumulative run time. 

__DESCRIPTION__
)

module load GROMACS/2018.1-nsc2-gcc-2018a-eb

GROMPP="gmx_mpi grompp"
MDRUN="gmx_mpi mdrun"

# These will be looked for before running
DEPENDENCIES=(backward.py gmx_mpi)



# Get script directory
#
# A: Test if program name contains a slash        (1) 
#    If so, change directory to program directory (2)
# B: Echo working directory, which will always be 
#    the program directory
#
#       |---------------A---------------|  |B|
#       |--------1--------|    |---2----|
SDIR=$( [[ $0 != ${0%/*} ]] && cd ${0%/*}; pwd )
SDIR=$(pwd)

# - main options:
INP=
TOP=
OUT=backmapped.gro
OTP=backmapped.top
NDX=backmapped.ndx
RAW=projected.gro
BW=0-backward.gro
CG=martini
AA=gromos54a7
NUM=1
NP=0
KICK=0.0
KEEP=false
TRJ=true
STSTEP=1
EMSTEPS=0
NBSTEPS=2000
#MDSTEPS=5000
###0.000fs
DT=0.0002,0.0005,0.001,0.002
MDSTEPS=1000,1000,1000,800000
MDP=()
POSRE=false
export NSTXTCOUT=1000
export GMX_MAXCONSTRWARN=-1
GMX_MAXCONSTRWARN=-1



# - em/md options

OPTIONS=$(cat << __OPTIONS__

## OPTIONS ##

  INITRAM options:
    -f       Input coarse grained structure                               *FILE:  None
    -p       Input atomistic target topology                              *FILE:  None
    -po      Output processed topology                                    *FILE:  $OTP
    -o       Output atomistic structure                                   *FILE:  $OUT
    -from    Coarse grained force field                                   *STR:   $CG
    -to      Target forcefield                                            *STR:   $AA
    -n       Number of runs                                               *INT:   $NUM
    -np      Number of processors/threads                                 *INT:   $NP
    -fc      Position restraint force constant                            *FLOAT: None
    -kick    Random kick size                                             *FLOAT: $KICK
    -keep    Do not delete intermediate files                             *BOOL:  false
    -trj     Write a trajectory for each stage of relaxation              *BOOL:  false
    -em      Number of steps for EM with bonded interactions only         *INT:   $EMSTEPS
    -nb      Number of steps for EM with nonbonded interactions           *INT:   $NBSTEPS
    -md      Number of steps for MD cycles                                *INT:   $MDSTEPS
    -dt      Time steps for successive MD cycles (comma separated list)   *STR:   $DT                          
    -nopr    Disable position restraints                                  *BOOL:  false
    -mdp     User provided MDP file for post-hoc simulations              *FILE:  None
             Multiple instances possible
__OPTIONS__
)


USAGE ()
{
    cat << __USAGE__

$PROGRAM version $VERSION:

$DESCRIPTION

$OPTIONS

(c)$YEAR $AUTHOR
$AFFILIATION

__USAGE__
}


BAD_OPTION ()
{
  echo
  echo "Unknown option "$1" found on command-line"
  echo "It may be a good idea to read the usage:"
  echo -e $USAGE

  exit 1
}


while [ -n "$1" ]; do
  case $1 in
        -h)                   USAGE      ; exit 0 ;;
    # File options
        -f)	              INP=$2      ; shift 2; continue ;;
        -p)                   TOP=$2      ; shift 2; continue ;;
       -po)                   OTP=$2      ; shift 2; continue ;;
        -o)                   OUT=$2      ; shift 2; continue ;;
        -n)                   NUM=$2      ; shift 2; continue ;;
       -np)                    NP=$2      ; shift 2; continue ;;
     -from)                    CG=$2      ; shift 2; continue ;;
       -to)                    AA=$2      ; shift 2; continue ;;
     -kick)                  KICK=$2      ; shift 2; continue ;;
     -keep)                  KEEP=true    ; shift  ; continue ;;
      -trj)                   TRJ=true    ; shift  ; continue ;;
       -em)               EMSTEPS=$2      ; shift 2; continue ;;
       -nb)               NBSTEPS=$2      ; shift 2; continue ;;
       -md)               MDSTEPS=$2      ; shift 2; continue ;; 
       -dt)                    DT=$2      ; shift 2; continue ;; 
     -nopr)                 POSRE=false   ; shift  ; continue ;;
      -mdp)       MDP[${#MDP[@]}]=$2      ; shift 2; continue ;;
         *)       BAD_OPTION $1;;
  esac
done


cat << __RUNINFO__

$PROGRAM version $VERSION:

(c)$YEAR $AUTHOR
$AFFILIATION

Now executing...

$CMD

__RUNINFO__


# Bookkeeping

# Check input
[[ -z $INP ]] && echo FATAL ERROR: MISSING INPUT STRUCTURE FILE && exit 1
[[ -z $TOP ]] && echo FATAL ERROR: MISSING INPUT TOPOLOGY FILE && exit 1

# Check dependencies
missing=false
echo Checking dependencies:
for i in ${DEPENDENCIES[@]}
do
  echo -n "$i ... "
  which $i || which $SDIR/$i || missing=true
  $missing && echo Missing dependency: $i && exit 1
done

# Array for temporary files
GARBAGE=()

# Function for trashing temporary files
trash()
{
    for item in $@; do GARBAGE[${#GARBAGE[@]}]=$item; done
}


# Define for position restraints
$POSRE && MDPDEF=-DPOSRES || MDPDEF=

# Starting time
START=$(date +%s)

WATER=spc216_one_NB.gro

###########################################
## STEP 1: GENERATING STARTING STRUCTURE ##
###########################################

echo '*********** STEP 1 **********'

GRO=$BW

if [ $STSTEP -le 1 ]
then
	B="$SDIR/backward.py -f $INP -raw $RAW -o $GRO -kick $KICK -sol -p $TOP -po $OTP -n $NDX -from $CG -to $AA"
	echo $B; $B || exit
	
	#Match topology and gromacs atom names
	#python rename.py $OTP $GRO

	#Solvate the box
	#mpprun -n 1 gmx_mpi solvate -cp $GRO -cs $WATER -o $GRO -p $OTP

	#Create index file
	#mpprun -n 1 gmx_mpi make_ndx -f $GRO -o $NDX
fi

BM=$(date +%s)

trash $RAW $GRO 


########################################
## STEP 2: EM WITHOUT NB INTERACTIONS ##
########################################

echo '*********** STEP 2 **********'

# Step counter
i=1

# Basename
BASE=$i-EM

# Run parameter file
mdp=$BASE.mdp


# If TRJ is true then write each frame
$TRJ && NSTXTCOUT=1 || NSTXTCOUT=0



cat << __MDP__ > $mdp

;define=-DFLEXIBLE
integrator=steep
nsteps=$EMSTEPS
emstep=0.001
emtol=10.0
nstcgsteep=1000
pbc=xyz
nstenergy=10

; Table extension is needed initially (or not?)
table-extension=2

; Define energy group and cutoff-scheme
; NB interactions can be excluded within groups PDT and STYR
energygrps=Protein Membrane Solvent

; Freeze atoms
;freezegrps=Protein Membrane
;freezedim=Y Y Y Y Y Y 

; Usually no trajectory is written, 
; but this can be changed (-trj)
nstxout=$NSTXTCOUT

;Nicolas Rolland / modif for PEDOT-TOS GAFF
cutoff-scheme=Verlet
rlist=1.2
coulombtype=P3M-AD
rcoulomb=1.2
vdwtype=PME
rvdw=1.2
pme-order=4
fourierspacing=0.12
ewald-rtol=1e-5

__MDP__


# Maximum global warning allowed for Gromacs run
MAXWARN=10
# Maximum contraints warnings

if [ $STSTEP -le 2 ]
then
	# Set up for running
	G="srun --mpi=pmi2 -n 1 $GROMPP -f $mdp -c $GRO -n $NDX -p $OTP -o $BASE -maxwarn $MAXWARN"
	#G="gmx grompp -f $mdp -c $GRO -n $NDX -p $OTP -o $BASE -maxwarn $MAXWARN"
	echo $G; $G || exit

	# Run
	M="mpprun $MDRUN -deffnm $BASE -v -nt $NP -pin on"
	echo $M; $M 
fi

# Timing
EM1=$(date +%s)

# Mark files for deletion
trash $BASE.*

# Bookkeeping (increment run counter)
GRO=$((i++))-EM.gro


############################################
## STEP 3: EM WITH NONBONDED INTERACTIONS ##
############################################

echo '********** STEP 3 ********'


# Basename for this cycle
BASE=$i-EM

# Turn all nonbonded interactions on and set the table_extension to default
# Change the number of steps
sed -e '/^nsteps/s/=.*$/='$NBSTEPS'/' -e '/^ *table-extension/s/^/;/' -e '/^ *energygrp_excl/s/^/;/' $mdp > $BASE.mdp
mdp=$BASE.mdp

if [ $STSTEP -le 3 ]
then
	# Set up for running
	G="srun --mpi=pmi2 -n 1 $GROMPP -f $mdp -c $GRO -n $NDX -p $OTP -o $BASE -maxwarn $MAXWARN"
	#G="gmx grompp -f $mdp -c $GRO -n $NDX -p $OTP -o $BASE -maxwarn $MAXWARN"
	echo $G; $G || exit

	# Run
	M="mpprun $MDRUN -deffnm $BASE -v -nt $NP -pin on"
	echo $M; $M 
fi

# Timing
EM2=$(date +%s)

# Mark files for deletion
trash $BASE.*

# Bookkeeping (increment run counter)
GRO=$i-EM.gro




##############################################
## STEP 4: NVT MD WITH INCREASING TIME STEP ##
##############################################

echo '********** STEP 4 *********'

NSTXTCOUT=1000

# Array for collecting timing information
PRMD=()

# Set file tag
$POSRE && tag=mdpr || tag=mdnopr

# Unpack the list of time steps
ifs=$IFS
IFS=,
DT=($DT)
MDSTEPS=($MDSTEPS)
IFS=$ifs

# Count for rerun
COUNT=3
# Count for MDSTEPS
COUNT_MD=0

for DELTA_T in ${DT[@]}
do

    $((++COUNT))
    BASE=$((++i))-$tag-$DELTA_T 
    mdp=$BASE.mdp

    # Is the md a continuation?
    CONTINUATION=yes
    GENVEL=no

    if [ $COUNT == 4 ]
    then
	CONTINUATION=no
	GENVEL=yes
    fi

cat << __MDP__ > $mdp
define                   = $MDPDEF
integrator               = md
nsteps                   = ${MDSTEPS[$COUNT_MD]}
dt                       = $DELTA_T
pbc                      = xyz

;Output
nstxout= $NSTXTCOUT
nstvout=$NSTXTCOUT
nstfout=$NSTXTCOUT
nstxout-compressed=$NSTXTCOUT

; Freeze atoms
;energygrps=pedot tos SOL

;Neighbor search
rlist                    = 1.2
coulombtype              = P3M-AD
rcoulomb                 = 1.2
;vdwtype                  = PME
rvdw                     = 1.2
pme-order                = 4
fourierspacing           = 0.12
ewald-rtol               = 1e-5

;Temperature couling
tcoupl                   = v-rescale
ref_t                    = 300
tau_t                    = 0.1
tc_grps                  = System

;Bond parameters
continuation             = $CONTINUATION
constraint_algorithm     = lincs
constraints              = all-bonds
lincs_iter               = 1
lincs_order              = 4

;Pressure coupling
pcoupl                   = berendsen
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = 1.0
compressibility          = 4.5e-5
refcoord_scaling         = com


;Velocity generation
gen_vel                  = $GENVEL
gen_temp                 = 300

;Dispersion correction
DispCorr                 = EnerPres



__MDP__

if [ $STSTEP -le $COUNT ]
then
  # Set up run
  G="srun --mpi=pmi2 -n 1 $GROMPP -f $mdp -c $GRO -r $BW -p $OTP -o $BASE -maxwarn $MAXWARN"
  #G="gmx grompp -f $mdp -c $GRO -r $BW -p $OTP -o $BASE -maxwarn $MAXWARN"
  echo $G; $G || exit

  # Perform run
  M="mpprun $MDRUN -deffnm $BASE -v -nt $NP -pin on"
  echo $M; $M 
fi

  # Collect timing information
  PRMD[${#PRMD[@]}]=$(date +%s)  

  # Collect garbage
  trash $BASE.*

  # Increment counter
  GRO=$BASE.gro

  $((++COUNT_MD))


done



########################################################
## LAST STEP : PRODUCE FINAL FILE, GENERATE RUN STATS ##
########################################################

# Copy last structure for output 
cp $GRO $OUT

# Trash files if not keeping them
if ! $KEEP
then
    trash *mdout.mdp*

    echo Removing files:
    echo ${GARBAGE[@]}
    rm ${GARBAGE[@]}
fi


echo
echo "Timing (seconds):"
echo ===========================
printf "%-15s %5d %5d\n" Backmapping: $((BM-START)) $((BM-START))
printf "%-15s %5d %5d\n" EM1: $((EM1-BM)) $((EM1-START))
printf "%-15s %5d %5d\n" EM2: $((EM2-EM1)) $((EM2-START))
P=$EM2
for ((i=0; i<${#DT[@]}; i++))
do
    T=${PRMD[$i]}
    printf "%-15s %5d %5d\n" PRMD-${DT[$i]}: $((T-P)) $((T-START))
    P=$T
done
for ((i=0; i<${#MDP[@]}; i++))
do
    T=${MD[$i]}
    printf "%-15s %5d %5d\n" MD-${MDP[$i]}: $((T-P)) $((T-START))
    P=$T
done
echo ---------------------------
printf "%-15s %5d\n\n\n" Total: $((T-START))

