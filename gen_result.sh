#echo >>show_energy.csv
#echo >>show_energy.csv
#echo >>show_rel.csv
#echo >>show_rel.csv
#echo pevfs >> show_energy.csv
#echo pdvfs >> show_rel.csv
#for r in {9,8,81,82,83}
#do
	 
#	 python get_power.py $1/load${r}/rel_c1_pdvfs 360000 >> show_energy.csv	
#	 python get_rel.py $1/load${r}/rel_c1_pdvfs >> show_rel.csv

#done
path=pwd
echo >>show_energy.csv
echo >>show_energy.csv
echo >>show_rel.csv
echo >>show_rel.csv
echo ours >> show_energy.csv
echo ours >> show_rel.csv
for r in {9,91,92,93,94,8,81,82,83,84,85,86,87,7,71,72,73,74,75,76,77,78,6,61,62,63,64,65,5,51,52,53,54,55}
do
	 
	 python get_power.py $1/load${r}/rel_c1_ours_no_reserve 360000 >> show_energy.csv	
	 python get_rel.py $1/load${r}/rel_c1_ours_no_reserve >> show_rel.csv
	 
done
#9,91,92,93,94,8,81,82,83,84,85,86,87,7,71,72,73,74,75,76,77,78,6,61,62,63,64,65,5,51,52,53,54,55
#for r in {9,91,92,93,94,8,81,82,83,7,71,72,73,74,6,5,4,41,42,43,44,45,46,3,31,32,33,34}
#do
	 
#	 python get_power.py $1/load${r}/rel_c1_ours_expire 360000 >> show_energy.csv	
#	 python get_rel.py $1/load${r}/rel_c1_ours_expire >> show_rel.csv
	 
#done
#for r in {9,8,81,82,83,7,71,72,73,74,6,5,4,41,42,43,44,45,46,3,31,32,33,34,2,21,22,23,24}
#do
 
#	 python get_power.py results/load${r}/rel_c1_ndvfs 360000 >> show_energy.csv	
#	 python get_rel.py results/load${r}/rel_c1_ndvfs >> show_rel.csv
 
#done