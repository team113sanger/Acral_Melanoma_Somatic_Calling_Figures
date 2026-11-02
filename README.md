# Somatic variant calling analysis — acral melanoma PDX study

Scripts that generate the somatic-mutation figures of the acral melanoma
patient-derived xenograft (PDX) study: the mutational landscape of the
Brazilian cohort, and the comparison between patient tumours and their matched
PDX models.

## Repository layout

```
Somatic_variant_calling_analysis/
├── LICENSE
├── README.md
└── analysis/
    ├── renv.lock          package versions
    ├── .Rprofile          activates renv
    ├── data/              inputs 
    ├── results/           outputs 
    └── scripts/
        ├── Figure_2A.R
        └── Figure_3A_and_C.R
```

## Running the analyses

All scripts expect to be run from `analysis/`, the directory that holds
`renv.lock`. They use relative paths, so nothing needs editing.

```bash
git clone https://github.com/<user>/Somatic_variant_calling_analysis.git
cd Somatic_variant_calling_analysis/analysis

# place the input files in analysis/data/ — see Data availability below

R -e 'renv::restore()'          # once, installs the recorded package versions
Rscript scripts/Figure_2A.R
Rscript scripts/Figure_3A_and_C.R
```

If a script stops with `'data/' not found`, the working directory is wrong:
`cd` into `analysis/` and run it again.

## Scripts

### `Figure_2A.R`

Oncoplot of the Brazilian acral melanoma cohort, one tumour per patient,
labelled with the AM identifier.

Produces two figures: the top 20 COSMIC Cancer Gene Census genes, and a
pre-specified 20-gene list restricted to genes mutated in at least two samples.
The top bar shows mutations per Mb.

### `Figure_3A_C.R`

Patient tumour versus matched PDX at passage 1.

On its first run the script writes `data/mapk_pathway_genes.txt` from MSigDB and
reads it from disk thereafter, so the gene set is versioned with the analysis.

## Environment

R 4.4.1. Package versions are recorded in `analysis/renv.lock` and restored with
[`renv`](https://rstudio.github.io/renv/):

```R
if (!require("renv")) install.packages("renv")
renv::restore()
```

Opening R inside `analysis/` activates the project library automatically through
`.Rprofile`. Do not set `.libPaths()` or `RENV_PROJECT` manually — either will
override the recorded environment.

## Data availability

Download them into `analysis/data/` before running. Each script lists the files
it needs in its header and stops with a named error if one is missing.

The COSMIC Cancer Gene Census (v97) must be obtained separately from
<https://cancer.sanger.ac.uk/census>; it requires registration and cannot be
redistributed.

## Contact

Annie Cristhine Moraes Sousa-Squiavinato
annie.sousa@ensino.inca.gov.br · as81@sanger.ac.uk
