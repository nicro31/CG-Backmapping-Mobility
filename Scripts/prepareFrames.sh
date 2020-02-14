source Load_LOE_CTP_FRAG.sh
export OMP_NUM_THREADS=1

declare -a chainLength=("N6" "N12")
declare -a waterName=("Hy10" "Hy20" "Hy30")


cd /proj/morphology/users/x_nicro/Back_Mapping/Final/NPT

for n in "${chainLength[@]}"
do
    echo  "Going to folder $n"

    for w in "${waterName[@]}"
    do
	echo  "Going to folder $w"
	
	cd $n
	cd ch33
	cd $w
	cd Mobility

	#Launch backmapping
	LOE-CTP-FRAG prep AA_mod.gro config.txt AA.pdb

	cd ..
	cd ..
	cd .. 
	cd ..

    done
done

