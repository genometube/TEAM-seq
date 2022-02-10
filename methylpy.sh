echo ===`date`===
source /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/zhangchen3/anaconda3/bin/activate
conda activate methylpy
## Forward mapping
methylpy single-end-pipeline --read-files /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_zhangchen3/P18Z12200N0396_zhangchen3/LIANTI/data/A1_B1_B2_B3/LIANTI_A1.reads.fastq.gz --sample A1  --forward-ref /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/index/methylpy_index_minimap2_hg19_f.fasta --reverse-ref /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/index/methylpy_index_minimap2_hg19_f.fasta --ref-fasta /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_zhangchen3/P18Z12200N0396_zhangchen3/DB/methylation/methylpy_index_minimap2/hg19.fa --trim-reads False --sort-mem 2G --num-procs 16 --aligner minimap2 --aligner-options "-a -x map-ont --MD --secondary=no" --remove-clonal True --path-to-picard /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/zhangchen3/anaconda3/envs/methylpy/share/picard-2.25.2-0 --remove-chr-prefix False --num-upstream-bases 1 --keep-temp-files True --pbat True --min-mapq 15
echo ========== Forward mapping is done =========

## Extract unmapped reads
awk  '{if($2 !=4 ) print $1}' A1_libA_forward_strand_hits.sam |awk -F "\!" '{print $1}'| sort|uniq >A1_libA_forward_strand_hits.sam_1
awk  '{if($2 !=4 ) print $1}'  A1_libA_reverse_strand_hits.sam |awk -F "\!" '{print $1}'| sort|uniq >A1_libA__reverse_strand_hits.sam_1
cat A1_libA_forward_strand_hits.sam_1 A1_libA_reverse_strand_hits.sam_1 |sort |uniq >map_ID.txt
/share/app/seqkit/0.14.0-dev/seqkit grep -v -f  map_ID.txt /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_zhangchen3/P18Z12200N0396_zhangchen3/LIANTI/data/A1_B1_B2_B3/LIANTI_A1.reads.fastq.gz |gzip >A1_reads_unmap.fastq.gz
echo ========== Unmapped reads extracted =========

## Reverse mapping
methylpy single-end-pipeline --read-files A1_reads_unmap.fastq.gz  --sample A1_other  --forward-ref /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/index/methylpy_index_minimap2_hg19_r.fasta --reverse-ref /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/index/methylpy_index_minimap2_hg19_r.fasta --ref-fasta /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_zhangchen3/P18Z12200N0396_zhangchen3/DB/methylation/methylpy_index_minimap2/hg19.fa --trim-reads False --sort-mem 2G --num-procs 16 --aligner minimap2 --aligner-options "-a -x map-ont --MD --secondary=no" --remove-clonal True --path-to-picard /hwfssz4/BC_PUB/Software/07.User-defined/02.Research/zhangchen3/anaconda3/envs/methylpy/share/picard-2.25.2-0 --remove-chr-prefix False --num-upstream-bases 1 --keep-temp-files True --pbat True --min-mapq 15
echo ========== Reverse mapping is done =========
## Merge forward and reverse
samtools merge -h A1_processed_reads_no_clonal.bam A1_total.bam A1_processed_reads_no_clonal.bam A1_other_processed_reads_no_clonal.bam
echo ================ Merge is done =============
## call methy
methylpy call-methylation-state --input-file A1_total.bam --paired-end False --sample A1_total --ref-fasta /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/new_reverse_index/new_ref/hg19_up_ref.fa  --num-procs 8 --num-upstream-bases 1 --keep-temp-files True --remove-chr-prefix False  --unmethylated-control Control: --binom-test True

echo ===`date`===



echo ===`date`===
