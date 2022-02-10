#!/bin/bash


## Other element or histone average mathylation results were calculated in the same way

/share/app/bedtools/2.29.2/bin/bedtools intersect -a /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0020_chenwenfang/cpg/allc_LIANTI_A1_new_type3.bed -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/CpGisland_window.txt -wa -wb > allc_LIANTI_A1_ctype_CpG.bed
/share/app/bedtools/2.29.2/bin/bedtools intersect -a /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0020_chenwenfang/cpg/allc_LIANTI_A1_new_type3.bed -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/CpGisland_up4kb_window.txt -wa -wb >allc_LIANTI_A1_ctype_up4kb.bed
/share/app/bedtools/2.29.2/bin/bedtools intersect -a /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0020_chenwenfang/cpg/allc_LIANTI_A1_new_type3.bed -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/CpGisland_down4kb_window.txt -wa -wb >allc_LIANTI_A1_ctype_down4kb.bed

python /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/hg19/get_window_aver.methylation.py allc_LIANTI_A1_ctype_CpG.bed allc_LIANTI_A1_ctype_CpG.avermethy.bed 20
python /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/hg19/get_window_aver.methylation.py allc_LIANTI_A1_ctype_up4kb.bed allc_LIANTI_A1_ctype_up4kb.avermethy.bed 50
python /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/hg19/get_window_aver.methylation.py allc_LIANTI_A1_ctype_down4kb.bed allc_LIANTI_A1_ctype_down4kb.avermethy.bed 50


