declare -a chainLength=("N6")
declare -a waterName=("Hy10" "Hy20" "Hy30")

cd /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT

for n in "${chainLength[@]}"
do
    echo  "Going to folder $n"

    for w in "${waterName[@]}"
    do
	echo  "Going to folder $w"

	#Copy required scripts
	cp /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/Processing/Jobs/Stride_Frames.sh $n/ch33/$w
	cp /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/Processing/Jobs/Extract_Data.sh $n/ch33/$w
	cp /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/Processing/Jobs/Prepare_XRD.sh $n/ch33/$w
	cp /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT/Processing/Jobs/Compute_XRD.sh $n/ch33/$w
	
	cd $n
	cd ch33
	cd $w

	#Extract frames and run XRD calculations
	./Stride_Frames.sh
	./Extract_Data.sh
	./Prepare_XRD.sh
	sbatch Compute_XRD.sh

	cd ..
	cd ..
	cd .. 

    done
done



