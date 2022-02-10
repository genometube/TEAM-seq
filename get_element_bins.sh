#!/bin/bash

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K27ac_region.txt -n 50 -i winnum >H3K27ac_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K27ac_up5kb.txt -n 50 -i winnum >H3K27ac_up5kb_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K27ac_down5kb.txt -n 50 -i winnum >H3K27ac_down5kb_window.txt


/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K27me3_region.txt -n 50 -i winnum >H3K27me3_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K27me3_up5kb.txt -n 50 -i winnum >H3K27me3_up5kb_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K27me3_down5kb.txt -n 50 -i winnum >H3K27me3_down5kb_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K4me3_region.txt -n 50 -i winnum >H3K4me3_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K4me3_up5kb.txt -n 50 -i winnum >H3K4me3_up5kb_window.txt

/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/anno/H3K4me3_down5kb.txt -n 50 -i winnum >H3K4me3_down5kb_window.txt


#/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/Ara_CGI_up4kb.bed -n 50 -i winnum > /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/Ara_CGI_up4kb_bin.bed

#/share/app/bedtools/2.29.2/bin/bedtools makewindows -b /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/Ara_CGI_down4kb.bed -n 50 -i winnum > /zfssz4/BC_RD_P1/PROJECT/P19Z12200N0088_chenwenfang/P19Z12200N0088_chenwenfang/LIANTI_Methy/compare/CpG/permethy_rate/Ara_CGI_down4kb_bin.bed

