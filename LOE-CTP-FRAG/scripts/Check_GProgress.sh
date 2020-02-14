#Get all gaussian com to process
declare -a dirs
i=1
for d in */site*
do
    extension="${d##*.}"
    if [ $extension == "com" ]
    then
	dirs[i++]="${d%/*}"
    fi
done
i=$(($i-1))
echo "Number of Gaussian calculations to perform: $i"


countFound=0
countUnfound=0

for f in "${dirs[@]}"
do
	#echo $f
	cd $f
	fCom=$f'.com'

	file=$f'.pun'
	if [ -s "$file" ]
	then
		#echo "$file found."
		countFound=$(( $countFound + 1 ))
	else
		echo "$file not found."
		countUnfound=$(( $countUnfound + 1 ))
		#g09 $fCom
		#mv fort.7 $f'.pun'
	fi
	
	cd ..
done


echo "Number of calculations already performed: $countFound"
echo "Number of calculations not done yet: $countUnfound"
