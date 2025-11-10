# Pre-designed NSR Primer Sets

This directory contains **ready-to-use not-so-random (NSR) primer sets** designed for various species and experimental conditions. These pre-computed datasets allow you to immediately use optimized primers without running the design pipeline yourself.

## Table of Contents

- [Available Species](#available-species)
- [Design Parameters](#design-parameters)
- [File Naming Convention](#file-naming-convention)
- [Species-Specific NSR Sets](#species-specific-nsr-sets)
  - [Homo sapiens (Human)](#homo-sapiens-human)
  - [Mus musculus (Mouse)](#mus-musculus-mouse)
  - [Escherichia coli](#escherichia-coli)
  - [Xenopus laevis](#xenopus-laevis)

---

## Available Species

| Species | Common Name | Available Files |
|---------|-------------|-----------------|
| Homo sapiens | Human | 28+ variants |
| Mus musculus | Mouse | 12+ variants |
| Escherichia coli | E. coli | 8+ variants |
| Xenopus laevis | African clawed frog | 2 variants |

---

## Design Parameters

NSR sets are characterized by the following parameters:

### 1. Primer Length (N-mer)

The length of the oligonucleotide sequence:
- **5-mer**: 5 nucleotides (1,024 possible sequences)
- **6-mer**: 6 nucleotides (4,096 possible sequences) - **Most commonly used**
- **7-mer**: 7 nucleotides (16,384 possible sequences)

**Recommendation**: 6-mer provides optimal balance between specificity and coverage.

### 2. rRNA Subunits Excluded

The rRNA sequences used for exclusion filtering. Common subunits include:
- **5S rRNA**: Small ribosomal subunit
- **5.8S rRNA**: Eukaryotic-specific small RNA
- **12S rRNA**: Mitochondrial small subunit
- **16S rRNA**: Prokaryotic small subunit / Mitochondrial large subunit
- **18S rRNA**: Eukaryotic small subunit
- **23S rRNA**: Prokaryotic large subunit
- **28S rRNA**: Eukaryotic large subunit
- **45S rRNA**: Pre-rRNA precursor (contains 18S, 5.8S, 28S)

**Common combinations:**
- **Eukaryotes**: `45S16S12S5S` (comprehensive coverage)
- **Prokaryotes**: `23S16S5S`

### 3. Base Match Removal Threshold

Minimum number of consecutive complementary bases to exclude a primer:
- **0 bases**: No exclusion (all possible N-mers)
- **3 bases**: Lenient filtering
- **5 bases**: Moderate filtering
- **6 bases**: Stringent filtering - **Recommended**
- **7 bases**: Very stringent filtering

**Recommendation**: 6-base match removal provides strong rRNA depletion while maintaining good mRNA coverage.

### 4. Tm Filtering (trimSD)

Standard deviation-based filtering of melting temperature (Tm) to reduce bias:
- **No suffix**: All primers included
- **trimSD_0.5**: Keep primers within 0.5 SD of mean Tm
- **trimSD_1.0**: Keep primers within 1.0 SD of mean Tm (moderate filtering)
- **trimSD_1.5**: Keep primers within 1.5 SD of mean Tm
- **trimSD_2.5**: Keep primers within 2.5 SD of mean Tm (lenient filtering)

Lower trimSD values result in more uniform Tm distribution but fewer primers.

---

## File Naming Convention

Files follow a systematic naming pattern:

**Pattern:**  
`NSR_[Species]_[Length]mer_[ExclusionBases]basematchremove_[rRNASubunits][_trimSD_[SD]].csv`

**Components:**
- `[Species]` = Species name (e.g., Hsapiens, Mmusculus, Ecoli)
- `[Length]` = Primer length (5, 6, or 7)
- `[ExclusionBases]` = Number of consecutive bases for exclusion (0, 3, 5, 6, 7)
- `[rRNASubunits]` = Combination of rRNA targets (e.g., 45S16S12S5S)
- `[SD]` = Standard deviation for Tm filtering (optional: 0.5, 1.0, 1.5, 2.5)

### Examples:

**Human 6-mer with 6-base exclusion:**

`NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S.csv`

- Species: Homo sapiens
- Length: 6-mer
- Exclusion: 6 consecutive bases
- rRNA targets: 45S, 16S, 12S, 5S

**Mouse 6-mer:**

`NSR_Mmusculus_6mer_28S18S5S.csv`

- Species: Mus musculus
- Length: 6-mer
- rRNA targets: 28S, 18S, 5S

**With Tm filtering:**

`NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_1.0.csv`

- Additional Tm filtering with SD=1.0

### Tm Distribution Plots

Some entries include `.png` files showing Tm distribution histograms:

- `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_tm.png`
- `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_1.0_tm.png`

---

## Species-Specific NSR Sets

### Homo sapiens (Human)

#### Standard 6-mer Sets

| File | rRNA Targets | Base Exclusion | Tm Filter | Notes |
|------|--------------|----------------|-----------|-------|
| `NSR_Hsapiens.csv` | - | - | - | Basic set |
| `NSR_Hsapiens_6mer_0S.csv` | None | 0 | - | All possible 6-mers |
| `NSR_Hsapiens_6mer_28S18S16S12S.csv` | 28S, 18S, 16S, 12S | - | - | Nuclear + mitochondrial |
| `NSR_Hsapiens_6mer_28S18S16S12S5S.csv` | 28S, 18S, 16S, 12S, 5S | - | - | Comprehensive set |
| `NSR_Hsapiens_6mer_45S16S12S5S.csv` | 45S, 16S, 12S, 5S | - | - | Using pre-rRNA precursor |

#### Stringent Filtering (6-base match removal)

| File | Tm Filter | Description |
|------|-----------|-------------|
| `NSR_Hsapiens_6mer_0basematchremove_45S16S12S5S.csv` | - | No base exclusion baseline |
| `NSR_Hsapiens_6mer_5basematchremove_45S16S12S5S.csv` | - | 5-base match removal |
| `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S.csv` | - | **Recommended** 6-base removal |
| `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_0.5.csv` | 0.5 SD | Highly uniform Tm |
| `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_1.0.csv` | 1.0 SD | **Recommended** balanced set |
| `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_1.5.csv` | 1.5 SD | More primers, less uniform |
| `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_2.5.csv` | 2.5 SD | Lenient Tm filtering |

#### 5-mer and 7-mer Sets

| File | Description |
|------|-------------|
| `NSR_Hsapiens_5mer_28S18S16S12S5S.csv` | 5-mer comprehensive set |
| `NSR_Hsapiens_7mer_5basematchremove_45S16S12S5S.csv` | 7-mer with 5-base removal |
| `NSR_Hsapiens_7mer_6basematchremove_45S16S12S5S.csv` | 7-mer with 6-base removal |
| `NSR_Hsapiens_7mer_7basematchremove_45S16S12S5S.csv` | 7-mer with 7-base removal |

#### Tm Distribution Visualizations

- `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_tm.png`
- `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_0.5_tm.png`
- `NSR_Hsapiens_6mer_6basematchremove_45S16S12S5S_trimSD_1.0_tm.png`

---

### Mus musculus (Mouse)

#### Standard 6-mer Sets

| File | rRNA Targets | Description |
|------|--------------|-------------|
| `NSR_Mmusculus.csv` | - | Basic set |
| `NSR_Mmusculus_6mer_28S.csv` | 28S only | Large subunit only |
| `NSR_Mmusculus_6mer_28S18S.csv` | 28S, 18S | Nuclear rRNA |
| `NSR_Mmusculus_6mer_28S18S5S.csv` | 28S, 18S, 5S | Nuclear comprehensive |
| `NSR_Mmusculus_6mer_45S16S12S5S.csv` | 45S, 16S, 12S, 5S | **Recommended** comprehensive |
| `NSR_Mmusculusn6_28S18S.csv` | 28S, 18S | Alternative naming |

#### 5-mer and 7-mer Sets

| File | Description |
|------|-------------|
| `NSR_Mmusculus_5mer_28S18S5S.csv` | 5-mer with 28S, 18S, 5S |
| `NSR_Mmusculus_5mer_28S18S16S12S5S.csv` | 5-mer comprehensive |
| `NSR_Mmusculus_7mer_28S18S.csv` | 7-mer with 28S, 18S |
| `NSR_Mmusculus_7mer_28S18S5S.csv` | 7-mer comprehensive |

---

### Escherichia coli

#### Standard 6-mer Sets

| File | Base Exclusion | Description |
|------|----------------|-------------|
| `NSR_Ecoli.csv` | - | Basic set |
| `NSR_Ecoli_6mer_23S16S5S.csv` | - | Prokaryotic rRNA set |
| `NSR_Ecoli_6mer_23S16S5S_remove3p5match.csv` | 3 or 5 | Legacy naming |
| `NSR_Ecoli_6mer_5basematchremove_23S16S5S.csv` | 5 | Moderate filtering |
| `NSR_Ecoli_6mer_6basematchremove_23S16S5S.csv` | 6 | **Recommended** stringent |
| `NSR_Ecoli_remove5match.csv` | 5 | Alternative format |

#### 7-mer Sets

| File | Base Exclusion | Description |
|------|----------------|-------------|
| `NSR_Ecoli_7mer_5basematchremove_23S16S5S.csv` | 5 | 7-mer with moderate filtering |
| `NSR_Ecoli_7mer_6basematchremove_23S16S5S.csv` | 6 | 7-mer with stringent filtering |

---

### Xenopus laevis

| File | Description |
|------|-------------|
| `NSR_Xlaevis.csv` | Basic set |
| `NSR_Xlaevis_6mer_45S16S12S5So5Ss.csv` | Comprehensive 6-mer (includes both 5S variants) |

**Note**: X. laevis has two 5S rRNA genes (5So and 5Ss), both are included in the comprehensive set.


## File Formats
