# fasta mode
docker run --rm fs000:5000/harukao/python4makensr:2.7.9 python nsrmaker.py \
  fasta \
  -f rRNA/Hsapiens_rRNA_45S.fasta \
  -f rRNA/Hsapiens_rRNA_16S.fasta \
  -f rRNA/Hsapiens_rRNA_12S.fasta \
  -f rRNA/Hsapiens_rRNA_5S.fasta \
  -s Hsapiens \
  -n 6 \
  -e 6 \
  -t 1 \
  -o test/test2_Hsapiens_rRNA_45S16S12S5S


# Converntional mode (but have to add something)
docker run --rm fs000:5000/harukao/python4makensr:2.7.9 python nsrmaker.py \
  -r 45S \
  -r 16S \
  -r 12S \
  -r 5S \
  -s Hsapiens \
  -n 6 \
  -e 6 \
  -t 1 \
  -o test/test3_Hsapiens_rRNA_45S16S12S5S


diff test/test2_Hsapiens_rRNA_45S16S12S5S_6mer_6basematchremove.csv test/test3_Hsapiens_rRNA_45S16S12S5S_6mer_6basematchremove.csv

diff test/test2_Hsapiens_rRNA_45S16S12S5S_6mer_6basematchremove.csv results/NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S.csv
diff test/test3_Hsapiens_rRNA_45S16S12S5S_6mer_6basematchremove.csv results/NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S.csv
