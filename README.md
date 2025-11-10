# nsrmaker

A python tool to make not-so-random primers (NSRs).

## How to pull docker image to use `nsrmaker`

```bash
docker pull [username]/python4makensr:2.7.9

cd nsrmaker*

docker run --rm -v $PWD:/nsrmaker [username]/python4makensr:2.7.9 python /nsrmaker/nsrmaker.py -h
```

## Help

```bash
Usage:
    nsrmaker.py (-s species) (-n length) (-e exclude) (-r rRNA...) (-t trimsd) (-o outprefix)
    nsrmaker.py fasta (-f fasta...) (-s species) (-n length) (-e exclude) (-t trimsd) (-o outprefix)
    nsrmaker.py -h | --help
    nsrmaker.py -v | --version

Options:
    -f fasta...     List of fasta files to be removed
    -s species      Species
    -n length       Length of NSR
    -e exclude      x base matching to be removed
    -r rRNA...      List of rRNA to be removed
    -t trimsd       SD value to trim NSR by its Tm
    -o outprefix    Prefix of output files
    -h --help       Show this screen
    -v --version    Show version
```

## How to use

```bash
# On `nsrmaker` or `nsrmaker-master`

mkdir nsrmaker_out

docker run --rm -v $PWD:/nsrmaker -v $PWD/nsrmaker_out:/nsrmaker_out \
  [username]/python4makensr:2.7.9 python /nsrmaker/nsrmaker.py \
  fasta \ # FASTA mode
  -f /nsrmaker/rRNA/Mmusculus_rRNA_12S.fasta \ # You can add a FASTA file
  -f /nsrmaker/rRNA/Mmusculus_rRNA_16S.fasta \
  -f /nsrmaker/rRNA/Mmusculus_rRNA_18S.fasta \
  -f /nsrmaker/rRNA/Mmusculus_rRNA_28S.fasta \
  -f /nsrmaker/rRNA/Mmusculus_rRNA_45S.fasta \
  -f /nsrmaker/rRNA/Mmusculus_rRNA_5p8S.fasta \
  -f /nsrmaker/rRNA/Mmusculus_rRNA_5S.fasta \
  -s Mmusculus \ # Species
  -n 6 \ # Length of NSR
  -e 6 \ # x base matching to be removed
  -t 1 \ # trimsd
  -o /nsrmaker_out/Mmusculus_6mer_NSR # `Mmusculus_6mer_NSR` is arbitrary
```


## How to build an environment to use `nsrmaker`

```bash
# Build
docker build -t harukao/python4makensr:2.7.9 .

# Test
git clone git@github.com:rikenbit/nsrmaker.git
cd nsrmaker
docker run -v $PWD:/nsrmaker --rm harukao/python4makensr:2.7.9 python /nsrmaker/nsrmaker.py -h
```

## How to use (on eloood)

```bash
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
## License
Copyright (c) 2025 Hiroki Danno and Haruka Ozaki, Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach Released under the Artistic License 2.0.

## Maintainers
- Hiroki Danno
- Haruka Ozaki
