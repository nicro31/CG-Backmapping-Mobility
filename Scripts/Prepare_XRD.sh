declare -a frameName=("0ps.data" "50ps.data" "100ps.data" "150ps.data" "200ps.data" "250ps.data" "300ps.data" "350ps.data" "400ps.data" "450ps.data" "500ps.data" "550ps.data" "600ps.data" "650ps.data" "700ps.data" "750ps.data" "800ps.data")


cd XRD

	
for f in "${frameName[@]}"
do

scriptName=$f".in"

cat << __SCRIPT__ > $scriptName

# ----------------- Atom Definition Section -----------------
units real
atom_style full
read_data $f
atom_modify sort 0 2.0

# ----------------- Run Section -----------------
compute 2 all xrd 1.540593 H H S C C H O C C C S O H 2Theta 1 33 LP 0 echo
fix 2 all ave/histo/weight 1 1 1 1 33 330 c_2[1] c_2[2] mode vector file $f.xrd

run 0

__SCRIPT__


done

cd ..

