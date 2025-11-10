# nsrmaker

`nsrmaker` is a Python tool for designing high-quality not-so-random primers (NSRs), mainly used to minimize rRNA bias in RNA-seq and full-length cDNA library preparation. By systematically excluding oligonucleotide sequences complementary to specific rRNA regions (from FASTA or built-in database), it enables efficient and non-biased cDNA synthesis targeting mRNA molecules.

Key Features
* Automatic generation of all possible oligo sequences of designated length (e.g., N=6-mers)
* Rapid exclusion of candidates matching rRNA sequences (contiguous base match threshold adjustable)
* Integrated calculation and histogram visualization of melting temperatures (Tm) for quality control
* Selection of optimal primer set using standard deviation-based filtering (trimSD parameter)
* Flexible input options (fasta mode or rRNA subunit names; custom output prefixes)
* Full reproducibility via Docker containers

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
docker build -t [username]/python4makensr:2.7.9 .

# Test
git clone git@github.com:rikenbit/nsrmaker.git
cd nsrmaker
docker run -v $PWD:/nsrmaker --rm [username]/python4makensr:2.7.9 python /nsrmaker/nsrmaker.py -h
```

## Example: How to use (on our environment)

```bash
docker run --rm fs000:5000/[username]/python4makensr:2.7.9 python PATH_TO_nsrmaker2/nsrmaker2/nsrmaker.py \
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

## Citation
If you use nsrmaker in your research, please cite the following paper:

[Hayashi T, Ozaki H, Sasagawa Y, Umeda M, Danno H, Nikaido I. Single-cell full-length total RNA sequencing uncovers dynamics of recursive splicing and enhancer RNAs. Nat Commun. 2018 Feb 12;9(1):619. doi: 10.1038/s41467-018-02866-0.](https://www.nature.com/articles/s41467-018-02866-0)

## License
Copyright (c) 2025 Hiroki Danno and Haruka Ozaki, Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach Released under the Artistic License 2.0.

## Maintainers
- Hiroki Danno
- Haruka Ozaki
