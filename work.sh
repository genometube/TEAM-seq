#!/usr/bash

for i in `cat barcode_name.txt`
do 

echo " cat xa*_out/"$i".fq |gzip >cell56/"$i".fq.gz"  >$i"_gz.sh"
#echo "bash Demultiplexing.sh" $i  >$i".sh"

done


#for i in `cat name.txt`
#do
#echo "bash Demultiplexing.sh" $i  >$i".sh"
#done
