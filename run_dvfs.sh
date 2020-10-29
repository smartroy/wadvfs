#for r in {9,91,92,93,94,8,81,82,83,7,71,72,73,74,6,5,4,41,42,43,44,45,46,3,31,32,33,34}
#do
#	mkdir results/load${r}
#	cp wkld${r} results/load${r}/wkld${r}
	
	#./test_core wkld${r} pdvfs
	#mv power_out_test results/load${r}/power_out_test_pdvfs
	#mv rel_c1 results/load${r}/rel_c1_pdvfs
	
#	./test_core wkld${r} wadvfs
#	mv power_out_test results/load${r}/power_out_test_ours_expire
#	mv rel_c1 results/load${r}/rel_c1_ours_expire

#done

#for r in {91,92,93,94,8,81,82,83,84,85,86,87,7,71,72,73,74,75,76,77,78,6,61,62,63,64,65,5,51,52,53,54,55}
for r in {1..100}
do
	#mkdir results/load${r}
	#cp wkld${r} results/load${r}/wkld${r}
	
	#./test_core wkld${r} pdvfs
	#mv power_out_test results/load_tmpr${r}/power_out_test_pdvfs
	#mv rel_c1 results/load_tmpr${r}/rel_c1_pdvfs
	#mv tmpr_c1 results/load_tmpr${r}/tmpr_c1_pdvfs
	for m in {1..5}
	do
		for p in {10..15}
		do
			./test_core wklds_th20/wkld${r}_${m} wadvfs $p
			mv rel_c1 results/random_test_th15/rel_ours_${r}_${m}_${p}
			#./test_core wklds/wkld${r}_${m} pdvfs $p
			#mv rel_c1 results/random_test/rel_pdvfs_${r}_${m}_${p}
		done
	done
	#mv power_out_test results/load${r}/power_out_test_ours_all
	
	#mv tmpr_c1 results/load${r}/tmpr_c1_ours
done
#9,91,92,93,94,8,81,82,83,7,71,72,73,74,6,5,4,41,42,43,44,45,46,3,31,32,33,34,75,76,77,78,61,62,63,64,65,5,51,52,53,54,55
#9,91,92,93,94,8,81,82,83,7,71,72,73,74,6,5,75,76,77,78,61,62,63,64,65,5,51,52,53,54,55