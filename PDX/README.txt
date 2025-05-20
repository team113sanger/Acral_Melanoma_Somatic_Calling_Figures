Fri 14 Mar 18:57:33 GMT 2025

Directories and files:

files.log
- a log file of files copied and coverted to xlsx format

NOTE: The DOMAINS column in the .xlsx files may be trunctated to 1000 characters.
(Excel has a ~32,000 character limit, but here we limit it to 1000. Use the original MAF
as to view all DOMAINS informtion)

# Main directories

one_tumour_per_patient
- variants from tumours when considering only one tumour per patient 

independent_tumours (if applicable)
- variants from all tumours, including multiple independent tumours from the same patient.
- does NOT include samples taken from the SAME tumour unless specified by Ingrid/Louise.

all_tumours
- variants from all tumours, including duplicates (samples from the same tumours) 

# Sub-directories

all_samples
- variants from all matched t/n pairs and unmatched tumours
- *keep.maf *keep.xlsx - all variants in coding and splice sites regions (including synonymous variants)
- *keepPA.maf *keepPA.xslx - protein-altering and splice site altering variants only

matched_samples
- variants from matched t/n pairs only
- *keep.maf *keep.xlsx - all variants in coding and splice sites regions (including synonymous variants)
- *keepPA.maf *keepPA.xslx - protein-altering and splice site altering variants only

QC_keep
QC_keepPA
- Quality control plots showing VAF vs depth for variant calls, variant call types,
top 12 recurrent sites and genes. FOR QC PURPOSES ONLY!

# Other files

sample_list.tsv 
- tumour/normal pairs used for analysis, including unmatched and insilico normal PDv38is_wes_v2

sample_list_matched.tsv
- matched tumour/normal pairs used for analysis

mutations_per_megabase.tsv
- mutation rate by sample, calculated for 'keep' MAFs only


For information on MAF files used for DERMATLAS:

https://drive.google.com/drive/u/1/folders/1CjKEPyBFIm3hwbKi1ja4xybMyF6uOMOC

------------------------------------------------------------
release_v4
├── all_tumours
│   ├── all_samples
│   │   ├── 6633_3248-filtered_mutations_all_allTum_keep.maf
│   │   ├── 6633_3248-filtered_mutations_all_allTum_keepPA.maf
│   │   ├── 6633_3248-filtered_mutations_all_allTum_keepPA.xlsx
│   │   ├── 6633_3248-filtered_mutations_all_allTum_keep.xlsx
│   │   ├── mutations_per_Mb.tsv
│   │   ├── QC_keep
│   │   │   ├── AF_vs_depth_indel_samples.pdf
│   │   │   ├── AF_vs_depth_recurrent_genes.pdf
│   │   │   ├── AF_vs_depth_recurrent_sites.pdf
│   │   │   ├── AF_vs_depth_snv_samples.pdf
│   │   │   ├── gene_tileplot_gene_order.txt
│   │   │   ├── gene_tileplot.pdf
│   │   │   ├── gene_tileplot_sample_order.txt
│   │   │   ├── gene_tileplot.tsv
│   │   │   ├── mutation_types_barplot_proportion_samples.pdf
│   │   │   ├── mutation_types_barplot_samples.pdf
│   │   │   ├── top_genes.list
│   │   │   ├── top_recurrently_mutated_genes.tsv
│   │   │   └── top_recurrently_mutated_sites.tsv
│   │   ├── QC_keepPA
│   │   │   ├── AF_vs_depth_indel_samples.pdf
│   │   │   ├── AF_vs_depth_recurrent_genes.pdf
│   │   │   ├── AF_vs_depth_recurrent_sites.pdf
│   │   │   ├── AF_vs_depth_snv_samples.pdf
│   │   │   ├── gene_tileplot_gene_order.txt
│   │   │   ├── gene_tileplot.pdf
│   │   │   ├── gene_tileplot_sample_order.txt
│   │   │   ├── gene_tileplot.tsv
│   │   │   ├── mutation_types_barplot_proportion_samples.pdf
│   │   │   ├── mutation_types_barplot_samples.pdf
│   │   │   ├── top_genes.list
│   │   │   ├── top_recurrently_mutated_genes.tsv
│   │   │   └── top_recurrently_mutated_sites.tsv
│   │   └── sample_list.tsv
│   └── matched_samples
│       ├── 6633_3248-filtered_mutations_matched_allTum_keep.maf
│       ├── 6633_3248-filtered_mutations_matched_allTum_keepPA.maf
│       ├── 6633_3248-filtered_mutations_matched_allTum_keepPA.xlsx
│       ├── 6633_3248-filtered_mutations_matched_allTum_keep.xlsx
│       ├── mutations_per_Mb.tsv
│       ├── QC_keep
│       │   ├── AF_vs_depth_indel_samples.pdf
│       │   ├── AF_vs_depth_recurrent_genes.pdf
│       │   ├── AF_vs_depth_recurrent_sites.pdf
│       │   ├── AF_vs_depth_snv_samples.pdf
│       │   ├── gene_tileplot_gene_order.txt
│       │   ├── gene_tileplot.pdf
│       │   ├── gene_tileplot_sample_order.txt
│       │   ├── gene_tileplot.tsv
│       │   ├── mutation_types_barplot_proportion_samples.pdf
│       │   ├── mutation_types_barplot_samples.pdf
│       │   ├── top_genes.list
│       │   ├── top_recurrently_mutated_genes.tsv
│       │   └── top_recurrently_mutated_sites.tsv
│       ├── QC_keepPA
│       │   ├── AF_vs_depth_indel_samples.pdf
│       │   ├── AF_vs_depth_recurrent_genes.pdf
│       │   ├── AF_vs_depth_recurrent_sites.pdf
│       │   ├── AF_vs_depth_snv_samples.pdf
│       │   ├── gene_tileplot_gene_order.txt
│       │   ├── gene_tileplot.pdf
│       │   ├── gene_tileplot_sample_order.txt
│       │   ├── gene_tileplot.tsv
│       │   ├── mutation_types_barplot_proportion_samples.pdf
│       │   ├── mutation_types_barplot_samples.pdf
│       │   ├── top_genes.list
│       │   ├── top_recurrently_mutated_genes.tsv
│       │   └── top_recurrently_mutated_sites.tsv
│       └── sample_list_matched.tsv
├── files.log
├── independent_tumours
│   ├── all_samples
│   │   ├── 6633_3248-filtered_mutations_all_indepTum_keep.maf
│   │   ├── 6633_3248-filtered_mutations_all_indepTum_keepPA.maf
│   │   ├── 6633_3248-filtered_mutations_all_indepTum_keepPA.xlsx
│   │   ├── 6633_3248-filtered_mutations_all_indepTum_keep.xlsx
│   │   ├── mutations_per_Mb.tsv
│   │   ├── QC_keep
│   │   │   ├── AF_vs_depth_indel_samples.pdf
│   │   │   ├── AF_vs_depth_recurrent_genes.pdf
│   │   │   ├── AF_vs_depth_recurrent_sites.pdf
│   │   │   ├── AF_vs_depth_snv_samples.pdf
│   │   │   ├── gene_tileplot_gene_order.txt
│   │   │   ├── gene_tileplot.pdf
│   │   │   ├── gene_tileplot_sample_order.txt
│   │   │   ├── gene_tileplot.tsv
│   │   │   ├── mutation_types_barplot_proportion_samples.pdf
│   │   │   ├── mutation_types_barplot_samples.pdf
│   │   │   ├── top_genes.list
│   │   │   ├── top_recurrently_mutated_genes.tsv
│   │   │   └── top_recurrently_mutated_sites.tsv
│   │   ├── QC_keepPA
│   │   │   ├── AF_vs_depth_indel_samples.pdf
│   │   │   ├── AF_vs_depth_recurrent_genes.pdf
│   │   │   ├── AF_vs_depth_recurrent_sites.pdf
│   │   │   ├── AF_vs_depth_snv_samples.pdf
│   │   │   ├── gene_tileplot_gene_order.txt
│   │   │   ├── gene_tileplot.pdf
│   │   │   ├── gene_tileplot_sample_order.txt
│   │   │   ├── gene_tileplot.tsv
│   │   │   ├── mutation_types_barplot_proportion_samples.pdf
│   │   │   ├── mutation_types_barplot_samples.pdf
│   │   │   ├── top_genes.list
│   │   │   ├── top_recurrently_mutated_genes.tsv
│   │   │   └── top_recurrently_mutated_sites.tsv
│   │   └── sample_list.tsv
│   └── matched_samples
│       ├── 6633_3248-filtered_mutations_matched_indepTum_keep.maf
│       ├── 6633_3248-filtered_mutations_matched_indepTum_keepPA.maf
│       ├── 6633_3248-filtered_mutations_matched_indepTum_keepPA.xlsx
│       ├── 6633_3248-filtered_mutations_matched_indepTum_keep.xlsx
│       ├── mutations_per_Mb.tsv
│       ├── QC_keep
│       │   ├── AF_vs_depth_indel_samples.pdf
│       │   ├── AF_vs_depth_recurrent_genes.pdf
│       │   ├── AF_vs_depth_recurrent_sites.pdf
│       │   ├── AF_vs_depth_snv_samples.pdf
│       │   ├── gene_tileplot_gene_order.txt
│       │   ├── gene_tileplot.pdf
│       │   ├── gene_tileplot_sample_order.txt
│       │   ├── gene_tileplot.tsv
│       │   ├── mutation_types_barplot_proportion_samples.pdf
│       │   ├── mutation_types_barplot_samples.pdf
│       │   ├── top_genes.list
│       │   ├── top_recurrently_mutated_genes.tsv
│       │   └── top_recurrently_mutated_sites.tsv
│       ├── QC_keepPA
│       │   ├── AF_vs_depth_indel_samples.pdf
│       │   ├── AF_vs_depth_recurrent_genes.pdf
│       │   ├── AF_vs_depth_recurrent_sites.pdf
│       │   ├── AF_vs_depth_snv_samples.pdf
│       │   ├── gene_tileplot_gene_order.txt
│       │   ├── gene_tileplot.pdf
│       │   ├── gene_tileplot_sample_order.txt
│       │   ├── gene_tileplot.tsv
│       │   ├── mutation_types_barplot_proportion_samples.pdf
│       │   ├── mutation_types_barplot_samples.pdf
│       │   ├── top_genes.list
│       │   ├── top_recurrently_mutated_genes.tsv
│       │   └── top_recurrently_mutated_sites.tsv
│       └── sample_list_matched.tsv
├── one_tumour_per_patient
│   ├── all_samples
│   │   ├── 6633_3248-filtered_mutations_all_oneTum_keep.maf
│   │   ├── 6633_3248-filtered_mutations_all_oneTum_keepPA.maf
│   │   ├── 6633_3248-filtered_mutations_all_oneTum_keepPA.xlsx
│   │   ├── 6633_3248-filtered_mutations_all_oneTum_keep.xlsx
│   │   ├── mutations_per_Mb.tsv
│   │   ├── QC_keep
│   │   │   ├── AF_vs_depth_indel_samples.pdf
│   │   │   ├── AF_vs_depth_recurrent_genes.pdf
│   │   │   ├── AF_vs_depth_recurrent_sites.pdf
│   │   │   ├── AF_vs_depth_snv_samples.pdf
│   │   │   ├── gene_tileplot_gene_order.txt
│   │   │   ├── gene_tileplot.pdf
│   │   │   ├── gene_tileplot_sample_order.txt
│   │   │   ├── gene_tileplot.tsv
│   │   │   ├── mutation_types_barplot_proportion_samples.pdf
│   │   │   ├── mutation_types_barplot_samples.pdf
│   │   │   ├── top_genes.list
│   │   │   ├── top_recurrently_mutated_genes.tsv
│   │   │   └── top_recurrently_mutated_sites.tsv
│   │   ├── QC_keepPA
│   │   │   ├── AF_vs_depth_indel_samples.pdf
│   │   │   ├── AF_vs_depth_recurrent_genes.pdf
│   │   │   ├── AF_vs_depth_recurrent_sites.pdf
│   │   │   ├── AF_vs_depth_snv_samples.pdf
│   │   │   ├── gene_tileplot_gene_order.txt
│   │   │   ├── gene_tileplot.pdf
│   │   │   ├── gene_tileplot_sample_order.txt
│   │   │   ├── gene_tileplot.tsv
│   │   │   ├── mutation_types_barplot_proportion_samples.pdf
│   │   │   ├── mutation_types_barplot_samples.pdf
│   │   │   ├── top_genes.list
│   │   │   ├── top_recurrently_mutated_genes.tsv
│   │   │   └── top_recurrently_mutated_sites.tsv
│   │   └── sample_list.tsv
│   └── matched_samples
│       ├── 6633_3248-filtered_mutations_matched_oneTum_keep.maf
│       ├── 6633_3248-filtered_mutations_matched_oneTum_keepPA.maf
│       ├── 6633_3248-filtered_mutations_matched_oneTum_keepPA.xlsx
│       ├── 6633_3248-filtered_mutations_matched_oneTum_keep.xlsx
│       ├── mutations_per_Mb.tsv
│       ├── QC_keep
│       │   ├── AF_vs_depth_indel_samples.pdf
│       │   ├── AF_vs_depth_recurrent_genes.pdf
│       │   ├── AF_vs_depth_recurrent_sites.pdf
│       │   ├── AF_vs_depth_snv_samples.pdf
│       │   ├── gene_tileplot_gene_order.txt
│       │   ├── gene_tileplot.pdf
│       │   ├── gene_tileplot_sample_order.txt
│       │   ├── gene_tileplot.tsv
│       │   ├── mutation_types_barplot_proportion_samples.pdf
│       │   ├── mutation_types_barplot_samples.pdf
│       │   ├── top_genes.list
│       │   ├── top_recurrently_mutated_genes.tsv
│       │   └── top_recurrently_mutated_sites.tsv
│       ├── QC_keepPA
│       │   ├── AF_vs_depth_indel_samples.pdf
│       │   ├── AF_vs_depth_recurrent_genes.pdf
│       │   ├── AF_vs_depth_recurrent_sites.pdf
│       │   ├── AF_vs_depth_snv_samples.pdf
│       │   ├── gene_tileplot_gene_order.txt
│       │   ├── gene_tileplot.pdf
│       │   ├── gene_tileplot_sample_order.txt
│       │   ├── gene_tileplot.tsv
│       │   ├── mutation_types_barplot_proportion_samples.pdf
│       │   ├── mutation_types_barplot_samples.pdf
│       │   ├── top_genes.list
│       │   ├── top_recurrently_mutated_genes.tsv
│       │   └── top_recurrently_mutated_sites.tsv
│       └── sample_list_matched.tsv
└── README.txt

21 directories, 194 files
