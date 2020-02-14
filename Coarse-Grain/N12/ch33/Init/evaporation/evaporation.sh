# initially bring md.gro from solution phase and change name.
cp ../md.gro md0.gro
# initially bring top.top from solution phase.
cp ../top.top top.top
line=187868
for ((i=1; i<=60; i++)); do

foo=(3820	3396	3038	2735	2474	2249	2054	1882	1732	1599	1480	1374	1280	1194	1117	1048	984	926	873	825	780	739	701	666	634	603	575	549	525	502	481	461	442	424	407	392	377	363	350	337	326	314	304	294	284	275	266	258	250	242	235	228	221	215	209	203	198	192	187	182)

a=${foo[i-1]}
b=$((line-(a*3)))
line=$b
ii=$((i-1))
echo $a
echo $b

grompp_mpi -f LALA.mdp -c md${ii}.gro -p top.top -o LALA.tpr
echo 4 | gmx_mpi genion -s LALA.tpr -o initial${i}.gro -p top.top -pname NAA -np $a -rmin 0.25

sed -i -n '/NAA/!p' initial${i}.gro
sed -i '2 c\  '$b' ' initial${i}.gro
sed -i -n '/NAA/!p' top.top

grompp_mpi -f npt.mdp -c initial${i}.gro -p top.top -o npt${i}.tpr
mpprun mdrun_mpi -rdd 1.4 -deffnm npt${i}

grompp_mpi -f md.mdp -c npt${i}.gro -t npt${i}.cpt -p top.top -o md${i}.tpr
mpprun mdrun_mpi -rdd 1.4 -deffnm md${i}

done
