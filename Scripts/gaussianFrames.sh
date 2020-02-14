source Load_LOE_CTP_FRAG.sh

declare -a chainLength=("N12")
declare -a waterName=("Hy5" "Hy10" "Hy20" "Hy30")


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

	#Launch gaussian
	gaussian.py 100 4 snic2017-12-59 01:00:00

	cd ..
	cd ..
	cd .. 
	cd ..

    done
done
