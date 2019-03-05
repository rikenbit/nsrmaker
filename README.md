# nsrmaker

A python tool to make not-so-random primers (NSRs).

## How to build an environment to use `nsrmaker`

```
# Build
docker build -t harukao/python4makensr:2.7.9 .

# Test
git clone git@github.com:rikenbit/nsrmaker.git
cd nsrmaker
docker run --rm harukao/python4makensr:2.7.9 python nsrmaker.py -h
```

## How to use
```
docker run --rm fs000:5000/harukao/python4makensr:2.7.9 python PATH_TO_nsrmaker2/nsrmaker2/nsrmaker.py \
  fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_12S.fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_16S.fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_18S.fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_28S.fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_45S.fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_5p8S.fasta \
  -f ../../bin/nsrmaker2/rRNA/Mmusculus_rRNA_5S.fasta \
  -f products/make_mm10_rmsk_rRNA_fasta/mm10_rmsk_rRNA.fa \
  -s Mmusculus \
  -n 6 \
  -e 6 \
  -t 1 \
  -o Mmusculus_6mer_NSR
```

## Maintainers
- Hiroki Danno
- Haruka Ozaki
