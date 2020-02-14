#!/usr/bin/env python3

# Needed modules
import time
import os
import os.path
import glob
import numpy as np
import argparse
import subprocess

# Time loop duration
dt = 60

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("nbChunck", help="Number of different jobs to launch",type=int)
parser.add_argument("nbSharedProc", help="Number of shared processors for each Gaussian job",type=int)
parser.add_argument("slurmAccount", help="Slurm account for the jobs, same format than for slurm",type=str)
parser.add_argument("jobLimit", help="Time limit for the jobs, same format than for slurm",type=str)
parser.add_argument("nbSharedProcBuild", help="Number of shared processors for building step",type=int)
parser.add_argument("jobLimitBuild", help="Time limit for building step",type=str)
args = parser.parse_args()

# Parameters for the calculation repartition
nbChunck = args.nbChunck
nbSharedProc = args.nbSharedProc
slurmAccount = args.slurmAccount
jobLimit = args.jobLimit
nbSharedProcBuild = args.nbSharedProcBuild
jobLimitBuild = args.jobLimitBuild

# List of calculation folders
calc_folders = ["/proj/morphology/users/x_nicro/T_Dependence/TvsM/Calculation/R2/300/frame9",
  "/proj/morphology/users/x_nicro/T_Dependence/TvsM/Calculation/R2/300/frame10",
  "/proj/morphology/users/x_nicro/T_Dependence/TvsM/Calculation/R2/300/frame11",
  "/proj/morphology/users/x_nicro/T_Dependence/TvsM/Calculation/R2/300/frame12",
  "/proj/morphology/users/x_nicro/T_Dependence/TvsM/Calculation/R2/300/frame13",
  "/proj/morphology/users/x_nicro/T_Dependence/TvsM/Calculation/R2/300/frame14"]


#############################
##### Program time loop #####
#############################
curr_cf = 0
while True and curr_cf < len(calc_folders):

	# Go to current folder to process
	os.chdir(calc_folders[curr_cf])
	# Save total number of calculations remaining
	nb_calc_remain = 0
	
	# Get all directories containing gaussian calculations to perform
	sep1 = '/'
	sep2 = '.'
	siteDir = glob.glob('*/site*')
	siteDir = [d.split(sep1,1)[0] for d in siteDir if d.split(sep2,1)[1]=='com']
	print("Total number of calculations -> ", len(siteDir))

	# Divide the calculation into chunck of smallest calculations that will correspond to different sbatch
	chuncks=np.array_split(np.array(siteDir), nbChunck)
	
	# Get a list of job running
	jobs = subprocess.getoutput("squeue -u x_nicro")
	joblist = jobs.split('\n')
	joblist = joblist[1:]
	joblist = [j.split()[2] for j in joblist]

	# Loop over chunks
	for idx,c in enumerate(chuncks):
		
		# Check if this chunk is already running
		statusName = "status" + str(idx) + ".txt"
		#isRunning = os.path.isfile(statusName)
		jobName = "run" + str(idx)
		isRunning = jobName in joblist
		if isRunning:
			nb_calc_remain += 1
			#f=open(statusName, "r")
			#contents = f.read()
			#if contents.rstrip() == "done":
			#	isRunning = False
			#else:
			#	nb_calc_remain += int(contents.rstrip())
			#f.close()
				

		# Launch calculation if chunk is not running
		if not isRunning:
			# First get only directories corresponding to calculations that are not done yet
			dirToDo = [d for d in c if not any(fname.endswith('.pun') for fname in os.listdir(d))]
			nbToDo = len(dirToDo)
			nb_calc_remain += nbToDo
			print("Chunk number ", idx, " -> ", nbToDo, " calc. left")

			# Launch calculations that remain to be done
			if nbToDo > 0:			
				fileName = "run" + str(idx) + ".sh"
				file = open(fileName,"w")
				file.write("%s\n" % "#!/usr/bin/env bash")  
				file.write("%s%s%s%i\n" % ("#SBATCH -U \"",slurmAccount,"\" -N 1 -n ",nbSharedProc))
				file.write("%s%s\n" % ("#SBATCH -t  ",jobLimit))
				#file.write("%s\n" % ("#SBATCH --signal=B:USR1@60"))
				file.write("%s%s\n" % ("#SBATCH -J ",jobName))
				file.write("%s\n\n" % ("SAVEDPWD=$(pwd)"))
				file.write("""term_handler()
{
        echo "function term_handler called." &
        
""")
				file.write("echo \"done\">$SNIC_TMP/status")
				file.write(str(idx))
				file.write(".txt &\n")
				file.write("cp $SNIC_TMP/status")
				file.write(str(idx))
				file.write(".txt ${SAVEDPWD}/ &\n")
				file.write("echo Exiting &\n wait \n")
				file.write("""
        #exit -1
}

# associate the function "term_handler" with the TERM signal
trap 'term_handler' USR1
""")


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
				file.write("declare -a dirs=(")
				for d in dirToDo:
					file.write("%s%s%s " % ("\"",d,"\""))
				file.write(")\n\n")
				file.write("""for f in "${dirs[@]}"
			do
				cd $f
				fCom=$f'.com'
				g09 $fCom
				mv fort.7 $f'.pun'
				#cp $SNIC_TMP/$f'.log' ${SAVEDPWD}/$f 
				#cp $SNIC_TMP/$f'.pun' ${SAVEDPWD}/$f 
				echo $fCom
				cd ..
			done
			""")
				file.write("echo \"done\">status")
				file.write(str(idx))
				file.write(".txt &\n wait \n")
				file.close()
				sbatchCommand = "sbatch " + fileName
				f=open(statusName, "w")
				f.write(str(nbToDo))
				f.close() 
				os.system(sbatchCommand)

	# Check if build is already running
	print("Total calculations remaining -> ", nb_calc_remain)
	if nb_calc_remain == 0:
		buildStatusName = "buildStatus.txt"
		isBuilt = False
		timeLimit = False
		isBuilding = os.path.isfile(buildStatusName)
		contents=""
		if isBuilding:
			f=open(buildStatusName, "r")
			contents = f.read()
			if contents.rstrip() == "done":
				isBuilt = True
			if contents.rstrip() == "timelimit":
				timeLimit = True;
		
		# If already built, go to next folder to process and free some space by removing sites folders
		if isBuilt:
			print("Build done for ", calc_folders[curr_cf])
			curr_cf += 1
			for s in siteDir:
				rm_command = "rm -rf " + s
				os.system(rm_command)
			
		
		# If not building, run the build calculation		
		elif not isBuilding or (isBuilding and timeLimit):
			print("Start building...")
			fileName = "build.sh"
			fileBuildStatus = open("buildStatus.txt","w")
			fileBuildStatus.write("running")
			fileBuildStatus.close()
			file = open(fileName,"w")
			file.write("%s\n" % "#!/usr/bin/env bash")  
			file.write("%s%s%s%i\n" % ("#SBATCH -U \"",slurmAccount,"\" -N 1 -n ",nbSharedProcBuild))
			file.write("%s%s\n\n" % ("#SBATCH -t  ",jobLimitBuild))
			
			file.write("%s\n" % ("#SBATCH --signal=B:USR1@60"))
			file.write("%s\n\n" % ("SAVEDPWD=$(pwd)"))
			file.write("""term_handler()
{
echo "function term_handler called." &

""")
			file.write("echo \"timelimit\">$SNIC_TMP/buildStatus")
			file.write(".txt &\n")
			file.write("cp $SNIC_TMP/buildStatus")
			file.write(".txt ${SAVEDPWD}/ &\n")
			file.write("echo Exiting &\n wait \n")
			file.write("""
exit -1
}

# associate the function "term_handler" with the TERM signal
trap 'term_handler' USR1
""")
			
			
			pdbFile = siteDir = glob.glob('*.pdb')
			file.write("%s%s%s%s%s\n%s\n" % ("LOE-CTP-FRAG build ",pdbFile[0]," config.txt ",pdbFile[0][:-3], "grp 1 &", "wait"))
			file.close()
			sbatchCommand = "sbatch " + fileName 
			os.system(sbatchCommand)

		else:
			print("Building in progress...")
			

	# Sleep until next monitoring
	print("\n")
	time.sleep(dt)

print("All calculation folders have been processed ! :D")
