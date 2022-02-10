#!/bin/bash


#samtools bedcov  CpGisland_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/A1/alignment/phat/add_control/A1_total.bam >A1_total_CpG.depth 

#samtools bedcov  CpGisland_up4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/A1/alignment/phat/add_control/A1_total.bam >A1_total_up4kb.depth

#samtools bedcov CpGisland_down4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/A1/alignment/phat/add_control/A1_total.bam >A1_total_down4kb.depth



samtools bedcov  CpGisland_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/3B/pipeline/alignment/B2/pbat/add_control/B2_total.bam >B2_total_CpG.depth


samtools bedcov  CpGisland_up4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/3B/pipeline/alignment/B2/pbat/add_control/B2_total.bam >B2_total_up4kb.depth



samtools bedcov CpGisland_down4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/3B/pipeline/alignment/B2/pbat/add_control/B2_total.bam  >B2_total_down4kb.depth


samtools bedcov  CpGisland_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/3B/pipeline/alignment/test/alignment/B3/add_control/B3_total.bam >B3_total_CpG.depth


samtools bedcov CpGisland_up4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/3B/pipeline/alignment/test/alignment/B3/add_control/B3_total.bam >B3_total_up4kb.depth


samtools bedcov CpGisland_down4kb_window.txt /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/3B/pipeline/alignment/test/alignment/B3/add_control/B3_total.bam >B3_total_down4kb.depth


