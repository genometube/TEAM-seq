#!/bin/bash

data="/zfswh1/BC_RD_P0/fangyitong/LIANTI/s6_scLIANTI/rawdata/$1.fastq.gz"
mkdir "/zfswh1/BC_RD_P0/fangyitong/LIANTI/s6_scLIANTI/s1_barcode/$1_out"
cd "/zfswh1/BC_RD_P0/fangyitong/LIANTI/s6_scLIANTI/s1_barcode/$1_out"
out_dir="/zfswh1/BC_RD_P0/fangyitong/LIANTI/s6_scLIANTI/s1_barcode/$1_out"


##Cut P5/P7 Adapter

/ifswh1/BC_PS/fangyitong/softwares/fastp -i ${data} -o "$1_1.fastq.gz" --adapter_sequence ATCTCGTATG --disable_quality_filtering &&\
/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_1.fastq.gz" -o "$1_11.fastq.gz" --adapter_sequence GTGTAGATCT --disable_quality_filtering &&\
/ifswh1/BC_RD/RD_COM/USER/chenweitian/software/seqtk-master/seqtk seq -r "$1_11.fastq.gz" |gzip > "$1_1r.fastq.gz" &&\
/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_1r.fastq.gz" -o "$1_2r.fastq.gz" --adapter_sequence GTGTAGATCT --disable_quality_filtering &&\
/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_2r.fastq.gz" -o "$1_22r.fastq.gz" --adapter_sequence ATCTCGTATG --disable_quality_filtering &&\
/ifswh1/BC_RD/RD_COM/USER/chenweitian/software/seqtk-master/seqtk seq -r "$1_22r.fastq.gz" |gzip >"$1_cutP5P7.fastq.gz" &&\
rm "fastp.html" "fastp.json" &&\
rm "$1_1.fastq.gz" "$1_11.fastq.gz" "$1_1r.fastq.gz" "$1_2r.fastq.gz" "$1_22r.fastq.gz" &&\
echo ========== Cut P5/P7 Adapter of "$1" is done ==========

## Demultiplexing

zcat "$1_cutP5P7.fastq.gz" |perl -nlE 'if(/runid/){say $_}else{if((length $_)>100){say substr($_,0,50).substr($_,-50)} else{say $_}}' |gzip > "$1_100.fastq.gz" &&\
/ifswh1/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2018a/software/bowtie2/bowtie2 -x /zfswh1/BC_RD_P0/P19Z12200N0095_chenwenfang/chenwenfang/Project/sc_LIANTI/bowtie2_hela_293/hela_293_barcode.fa -U "$1_100.fastq.gz" -S "$1.sam" --local --np 0 --mp 1,0 --rdg 0,1 --rfg 0,1 --score-min L,52,0 2>debarcode.rate &&\
/zfswh1/BC_RD_P0/zhangchen3/software/samtools/bin/samtools view -F4 -f16 "$1.sam" |cut -f1,3 |gzip > "$1_for.csv.gz" &&\
/zfswh1/BC_RD_P0/zhangchen3/software/samtools/bin/samtools view -F4 -F16 "$1.sam" |cut -f1,3 |gzip > "$1_rev.csv.gz" &&\
rm "$1_100.fastq.gz" "$1.sam" &&\
echo ========== Demultiplexing of "$1" is done ==========

## Data cleaning

/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_cutP5P7.fastq.gz" -o "$1_1.fastq.gz" --adapter_sequence CTATCTCTTATACACAT --disable_quality_filtering &&\
/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_1.fastq.gz" -o "$1_2.fastq.gz" --adapter_sequence TTGTTTTTTATATATATT --disable_quality_filtering &&\
/ifswh1/BC_RD/RD_COM/USER/chenweitian/software/seqtk-master/seqtk seq -r "$1_2.fastq.gz" |gzip > "$1_2r.fastq.gz" &&\
/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_2r.fastq.gz" -o "$1_3.fastq.gz" --adapter_sequence TTGTTTTTTATATATATT --disable_quality_filtering &&\
/ifswh1/BC_PS/fangyitong/softwares/fastp -i "$1_3.fastq.gz" -o "$1_4.fastq.gz" --adapter_sequence CTATCTCTTATACACAT --disable_quality_filtering &&\
/ifswh1/BC_RD/RD_COM/USER/chenweitian/software/seqtk-master/seqtk seq -r "$1_4.fastq.gz" |gzip > "$1_clean.fastq.gz" &&\
rm "fastp.html" "fastp.json" &&\
rm "$1_1.fastq.gz" "$1_2.fastq.gz" "$1_2r.fastq.gz" "$1_3.fastq.gz" "$1_4.fastq.gz" &&\
echo ========== Data cleaning of "$1" is done ==========
