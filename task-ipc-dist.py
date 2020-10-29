#!/usr/bin/python

import sys
import os

if __name__ == '__main__':
	input = sys.argv[1]
	wkld_dir = sys.argv[2]
	threshold = int(sys.argv[3])
	
	raw_dist = open(input,"r")

	raw_dist.readline()
	
	accu_dict={}
	new_dir = wkld_dir+'_th'+str(threshold)
	os.mkdir(new_dir)
	while 1:

		task = raw_dist.readline().rstrip()
		if not task:
			break
		task = task.split(",")
		
		accu_list=[0]
		for i in range(1,len(task)):
			accu_list.append(round(accu_list[i-1]+float(task[i]),1))
		
		accu_dict[task[0]] = accu_list
	
	for root,dirs,files in os.walk(wkld_dir):

		for name in files:
	# name = 'wkld1_1'
	# root = 'wklds'
			wkld = open(os.path.join(root,name),"r")
			new_wkld = open(os.path.join(new_dir,name),"w")
			while 1:
				task = wkld.readline().rstrip()
				if not task:
					break
				task = task.split()
				# print(task[-1])
				task[-1] = str(round(1 - accu_dict[task[0]][threshold+1],1))
				# print(type(threshold))
				# print(accu_dict[task[0]][threshold])
				# print(task)
				new_wkld.write(' '.join(task)+'\n')
			new_wkld.close()



