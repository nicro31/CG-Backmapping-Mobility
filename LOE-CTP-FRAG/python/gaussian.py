#!/usr/bin/env python3

import glob
import os
import numpy as np
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("nbChunck", help="Number of different jobs to launch",type=int)
parser.add_argument("nbSharedProc", help="Number of shared processors for each Gaussian job",type=int)
parser.add_argument("slurmAccount", help="Slurm account for the jobs, same format than for slurm",type=str)
parser.add_argument("jobLimit", help="Time limit for the jobs, same format than for slurm",type=str)
args = parser.parse_args()

# Parameters for the calculation repartition
nbChunck = args.nbChunck
nbSharedProc = args.nbSharedProc
slurmAccount = args.slurmAccount
jobLimit = args.jobLimit

# Get all directories containing gaussian calculations to perform
sep1 = '/'
sep2 = '.'
siteDir = glob.glob('*/site*')
print(siteDir)
siteDir = [d.split(sep1,1)[0] for d in siteDir if d.split(sep2,1)[1]=='com']
print("Total number of calculations")
print(len(siteDir))

# Get only directories corresponding to calculations that are not done yet
dirToDo = [d for d in siteDir if not any(fname.endswith('.pun') for fname in os.listdir(d))]
print("Calculations not done yet")
print(len(dirToDo))

# Divide the calculation into chunck of smallest calculations that will correspond to different sbatch
chuncks=np.array_split(np.array(dirToDo), nbChunck)
for idx,c in enumerate(chuncks):
	fileName = "run" + str(idx) + ".sh"
	file = open(fileName,"w")
	file.write("%s\n" % "#!/usr/bin/env bash")  
	file.write("%s%s%s%i\n" % ("#SBATCH -U \"",slurmAccount,"\" -N 1 -n ",nbSharedProc))
	file.write("%s%s\n\n" % ("#SBATCH -t  ",jobLimit))
	file.write("""#Import gaussian module depending on the machine it's running on
if [ $NSC_RESOURCE_NAME == "tetralith" ]; then
	module load Gaussian/09.E.01-avx-nsc1-bdist
# 	. $g09root/g09/bsd/g09.profile
    
elif [ $NSC_RESOURCE_NAME == "sigma" ]; then
    	module load Gaussian/09.E.01-avx-nsc1-bdist

elif [ $NSC_RESOURCE_NAME == "kebnekaise" ]; then
    	module load gaussian/09.e.01-AVX
else
    	echo "i'm operating on a not recognized machine, i don't know which gaussian09 module to load"
    	exit -1
fi

""")
	file.write("echo running > status")
	file.write(str(idx))
	file.write(".txt \n")
	file.write("declare -a dirs=(")
	for d in c:
		file.write("%s%s%s " % ("\"",d,"\""))
	file.write(")\n\n")
	file.write("""for f in "${dirs[@]}"
do
	cd $f
	fCom=$f'.com'
	g09 $fCom
	mv fort.7 $f'.pun'
	echo $fCom
	cd ..
done
""")
	file.write("echo done > status")
	file.write(str(idx))
	file.write(".txt \n")
	file.close()
	sbatchCommand = "sbatch " + fileName 
	os.system(sbatchCommand)


