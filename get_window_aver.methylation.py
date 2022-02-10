#!/bin/env python

'''
input----- element or hsitone window mc，all c counts
out-------  window average methy ratio

'''
import numpy
import sys

out=open(sys.argv[2],'w')

win={}
for i in range(int(sys.argv[3])):
	key=i+1
	win.setdefault(key,[0,0])

time={}
for line in open(sys.argv[1],'r'):
	aa=line[:-1].split('\t')
	if aa[3] != 'CG':
		continue
	key='_'.join([aa[-4],aa[-3],aa[-2],aa[-1]])
	ratio=int(aa[4])/int(aa[5])
	region.setdefault(key,[]).append(ratio)
	if key not in time:
		time[key]=int(aa[5])
	else:
		time[key] += int(aa[5])

for key in time:
	if time[key] <3 :
		continue
	me_rate=numpy.mean(region[key])
	num=int(key.split('_')[-1])
	if num in win:
		win[num][0]=win[num][0]+me_rate
		win[num][1]=win[num][1]+1
	else:
		continue

for key in win:
	out.write('\t'.join([str(key),str(win[key][0]),str(win[key][1]),str(win[key][0]/win[key][1])])+'\n')


out.close()
