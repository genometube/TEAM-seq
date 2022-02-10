## build single cell TEAM-seq barcodes bowtie2 index

awk '{print ">"$1"\n""NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGGAGA"$2 $3"AAATATATATAAAAAACAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"}' hela_293T_barcode.txt >hela_293_barcode.fa

/share/app/bowtie2-2.2.5/bowtie2-build /zfswh1/BC_RD_P0/P19Z12200N0095_chenwenfang/chenwenfang/Project/sc_LIANTI/bowtie2_hela_293/hela_293_barcode.fa /zfswh1/BC_RD_P0/P19Z12200N0095_chenwenfang/chenwenfang/Project/sc_LIANTI/bowtie2_hela_293/hela_293_barcode.fa
