#!/bin/bash


samtools bedcov  CpGisland_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/A1/alignment/phat/add_control/A1_total.bam >A1_total_CpG.depth 
samtools bedcov  CpGisland_up4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/A1/alignment/phat/add_control/A1_total.bam >A1_total_up4kb.depth
samtools bedcov CpGisland_down4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/A1/alignment/phat/add_control/A1_total.bam >A1_total_down4kb.depth

python /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/coverage_depth/get_window_depth_20.py A1_total_CpG.depth A1_total_CpG.aver.depth 
python /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/coverage_depth/get_window_depth_50.py A1_total_up4kb.depth A1_total_up4kb.aver.depth
python /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/coverage_depth/get_window_depth_50.py A1_total_down4kb.depth A1_total_down4kb.aver.depth







