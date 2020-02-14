declare -a chainLength=("N6")
declare -a waterName=("Hy10" "Hy20" "Hy30")


cd /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT

for n in "${chainLength[@]}"
do
    echo  "Going to folder $n"

    for w in "${waterName[@]}"
    do
	echo  "Going to folder $w"

	#Copy initram script
	cp /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/Processing/Jobs/initram.sh $n/ch33/$w
	
	cd $n
	cd ch33
	cd $w

	#Launch backmapping
	sbatch initram.sh -f CG.gro -o AA.gro -p atomistic.top -to amber96 -keep

	cd ..
	cd ..
	cd .. 

    done
done



