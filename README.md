# Palbociclib + IO for Patient Cohort A (p11044)

## Tumor and Immune Effects of CDK4/6 Inhibition Combined with PD-1 Blockade in Dedifferentiated Liposarcoma (DDLPS)

This repository contains the analysis code and data processing pipelines accompanying the manuscript:

Rosenbaum E, Gularte-Mérida R, Seffar E, Lee J, Adamow M, Bradic M, Dickson MA, Avutu V, Banks LB, Chan JE, Chi P, Gounder MM, Kelly CM, Keohan ML, Maki RG, Movva S, Reed DR, Desir R, Biniakewitz M, Erinjeri JP, Lefkowitz RA, Wong P, Antonescu CR, Qin L-X, Panageas KS, Shen R, Singer S, Koff A, Tap WD, D'Angelo SP. *Tumor and Immune Effects of CDK4/6 + retifanlimab in DDLPS.* Cancer Research Communications. 2025;6(2):437–446.

**DOI:** [10.1158/2767-9764.CRC-25-0334](https://doi.org/10.1158/2767-9764.CRC-25-0334)

**ClinicalTrials.gov:** [NCT04438824](https://clinicaltrials.gov/study/NCT04438824)

---

## Study Overview

This phase IIa clinical trial investigated the combination of palbociclib (CDK4/6 inhibitor) and retifanlimab (anti-PD-1) in patients with advanced dedifferentiated liposarcoma (DDLPS). Palbociclib had a 2-week lead-in followed by concomitant treatment with retifanlimab. Correlative analyses included single-cell RNA sequencing of tumor biopsies and high-parameter flow cytometry of peripheral blood to characterize changes in both the tumor microenvironment and circulating immune populations. 

---

## Data Availability

Raw and final gene count matrices are publicly available via **Gene Expression Omnibus**:

🗂️ **GEO Accession:** [GSE320212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320212)

---

## Analysis Overview

### scRNA-seq 

1. **Library preparation:** 10X Genomics Chromium platform; sequenced on 10X Genomics Chromium instrument to ~25,000 reads per cell
2. **Alignment:** Cell Ranger v7.1 aligned to GRCh38.p5 reference genome
3. **Quality control:**
   - Doublet removal via [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
   - Mitochondrial expression filtering via [miQC](https://github.com/greenelab/miQC)
   - Cell filtering: 256–5,792 expressed genes, ≥362 total reads, <50% ribosomal expression, <3% hemoglobin genes, <0.7% platelet markers
4. **Cell type classification:** [Azimuth](https://azimuth.hubmapconsortium.org/) using Adipose and PBMC references
5. **Clustering & analysis:** [Seurat v4](https://satijalab.org/seurat/) with SCTransform (20,000 cell subsampling, 35 PCs, regressing on percent mitochondrial expression)
6. **Cancer cell classification:** Six DDLPS-specific signatures (adipocyte differentiation, stemness, ECM remodeling, hypoxia, angiogenesis, invasion) per [Gruel et al., Nat Commun 2024](https://doi.org/10.1038/s41467-024-52428-6)
7. **Longitudinal clone analysis:** [Harmony](https://github.com/immunogenomics/harmony) integration for pre/post-retifanlimab biopsy alignment; differential expression via Seurat `FindMarkers` (Bonferroni-corrected P ≤ 0.01)

### Flow Cytometry

1. **Marker Panel:** 29-color panel on BD FACSymphony
2. **Cell Stratification:** T-cell gating and analysis using [staRgate](https://github.com/leejasme/staRgate)
3. **Cell Composition:** Topic modeling approach for T-cell composition analysis per [Peng et al., Cell Rep Methods 2023](https://doi.org/10.1016/j.crmeth.2023.100546)

---

## Software & Dependencies

| Tool | Version | RRID | Purpose |
|------|---------|------|---------|
| [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) | v7.1 | SCR_023221 | scRNA-seq alignment & quantification |
| [Seurat](https://satijalab.org/seurat/) | v4 | SCR_007322 | Single-cell analysis framework |
| [SCTransform](https://github.com/satijalab/sctransform) | — | SCR_022146 | Normalization |
| [Azimuth](https://azimuth.hubmapconsortium.org/) | — | SCR_021084 | Automated cell type annotation |
| [Harmony](https://github.com/immunogenomics/harmony) | — | SCR_022206 | Batch integration |
| [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) | — | — | Doublet detection |
| [miQC](https://github.com/greenelab/miQC) | — | — | Mitochondrial QC |
| [R](https://www.r-project.org/) | v4.0.5 | SCR_001905 | Statistical computing |
| [staRgate](https://github.com/leejasme/staRgate) | — | — | Flow cytometry gating |

---

## Repository Structure

```
palbo_io_p11044/
├── data-raw
│   └── DATASET.R
├── DESCRIPTION
├── LICENSE
├── NAMESPACE
├── R
│   ├── dubFinder.R
│   ├── figure.1.fx.R
│   ├── figure.3.fx.R
│   ├── helper_fx.R
│   └── prepare_seu.R
├── README.md
└── vignettes
    ├── 00_pre_proc_p11044.R
    ├── 01_split.cells.R
    ├── Figure_1.R
    ├── Figure_2.R
    ├── Figure_3_and_4.R
    ├── Figure_4_T_Cells.R
    ├── Figure3_Panel C_D_E_F.R
    ├── main.R
    ├── manuscript_flow_figs.R
    └── Supplementary_Figures.R

4 directories, 20 files
```

---

## Funding

This study was supported by:

- Geoffrey Beene Foundation
- Cycle for Survival
- FDA R01 FD007528
- NIH/NCI Cancer Center Support Grant P30 CA008748 (MSKCC)
- Marie-Josée and Henry R. Kravis Center for Molecular Oncology
- Additional funding and retifanlimab was provided by Incyte Corporation

---

## Citation

If you use code or data from this repository, please cite:

```bibtex
@article{Rosenbaum2025,
  title={Tumor and Immune Effects of CDK4/6 + retifanlimab in DDLPS},
  author={Rosenbaum, Eli and Gularte-M{\'e}rida, Ricardo and Seffar, Elias and Lee, Jasme and Adamow, Margaret and Bradic, Martina and Dickson, Mark A and Avutu, Vamsi and Banks, Liam B and Chan, John E and Chi, Ping and Gounder, Mrinal M and Kelly, Ciara M and Keohan, Mary Louise and Maki, Robert G and Movva, Sujana and Reed, Damon R and Desir, Reginald and Biniakewitz, Meredith and Erinjeri, Joseph P and Lefkowitz, Robert A and Wong, Phillip and Antonescu, Cristina R and Qin, Li-Xuan and Panageas, Katherine S and Shen, Ronglai and Singer, Samuel and Koff, Andrew and Tap, William D and D'Angelo, Sandra P},
  journal={Cancer Research Communications},
  volume={6},
  number={2},
  pages={437--446},
  year={2025},
  doi={10.1158/2767-9764.CRC-25-0334}
}
