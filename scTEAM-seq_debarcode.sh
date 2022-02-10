#!/usr/bash

###   debarcode and get clean data
bash pipe.sh raw_sc_data


###   split clean fastq to single cell fastq data
for i in  `cat barcode_name.txt`
do
zcat raw_sc_data_out/*.csv.gz | grep $i >$i'.csv'
/ifswh1/BC_RD/RD_COM/USER/chenweitian/software/seqtk-master/seqtk subseq raw_sc_data_out/raw_sc_data_clean.fastq.gz $i'.csv' |gzip > $i'.fq.gz'

done

