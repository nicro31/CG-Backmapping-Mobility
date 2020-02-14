source Load_LOE_CTP_FRAG.sh

declare -a chainLength=("N6" "N12" "N18")
declare -a waterName=("Hy30")


cd /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT

for n in "${chainLength[@]}"
do
    echo  "Going to folder $n"

    for w in "${waterName[@]}"
    do
	echo  "Going to folder $w"

	cp /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/Processing/Jobs/build.sh /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/$n/ch33/$w/Mobility
	
	cd $n
	cd ch33
	cd $w
	cd Mobility

	#Launch build
	sbatch build.sh

	cd ..
	cd ..
	cd .. 
	cd ..

    done
done
