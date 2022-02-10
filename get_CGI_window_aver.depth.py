#!/bin/env python

'''
input----- all CpG island and franking region window depth
out------- 20/50 window average depth

'''


import sys

out=open(sys.argv[2],'w')
win={}

for i in range(int(sys.argv[3])):
	key=i+1
	win.setdefault(key,[0,0])


#print (win)
time=0
for line in open(sys.argv[1],'r'):
	aa=line[:-1].split('\t')
	num=int(aa[3])
	if num in win:
		win[num][0]=win[num][0]+int(aa[-1])/(int(aa[2])-int(aa[1]))
		win[num][1]=win[num][1]+1
	
	else:
		continue

for key in win:
	out.write('\t'.join([str(key),str(win[key][0]),str(win[key][1]),str(win[key][0]/win[key][1])])+'\n')

#	print(key,win[key]/time)	


out.close()
