#Author: Annie Squiavinato
#OCT-2026
#______________________________________________________________________________
# Oncoplot of the Brazilian acral melanoma cohort, one tumour per patient.
# Two figures: the top COSMIC Cancer Gene Census genes, and a pre-specified
# gene list restricted to genes mutated in at least MIN_CASES samples.
# Samples are labelled with the AM identifier.
#
# INPUTS  
#   6633_2729_3248-filtered_mutations_matched_allTum_keepPA_KEEPONLY.maf
#       somatic mutations, matched tumour-normal, protein-altering only
#   clinical_info_one_sp_april_26.csv
#       one row per patient: ID_Sanger, ID_lab, Anatomic_site, Tumor_type,
#       Clinical_Stage
#   cancer_gene_census.v97.genes.tsv
#       COSMIC Cancer Gene Census v97, downloaded from
#       https://cancer.sanger.ac.uk/census (registration required)
#   6633_2729_3248-filtered_mutations_matched_allTum_keep_KEEPONLY_mutations_per_MB.tsv
#       mutations per Mb per sample, used for the top bar
#
# OUTPUTS 
#   oncoplot_one_per_patient_top20.pdf
#   oncoplot_genes_of_interest_min2.pdf
#   COSMIC_matched_samples_BRA_AM.txt   CGC genes mutated in the cohort
#   gene_frequency_all.csv              mutated samples per gene
#
# Post-processing in Affinity for the submitted figures.

# Dependencies are managed with renv: renv::restore() rebuilds the
# environment recorded in renv.lock.
#_______________________________________________________________________________

library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(maftools)

# Run from analysis/ — the directory that holds renv.lock, data/, scripts/
# and results/. Run renv::restore() once before the first use.
IN  <- "analysis/data"
OUT <- "analysis/results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the analysis/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_N     <- 20
MIN_CASES <- 2

F_MAF  <- file.path(IN, "6633_2729_3248-filtered_mutations_matched_allTum_keepPA_KEEPONLY.maf")
F_CLIN <- file.path(IN, "clinical_info_one_sp_april_26.csv")
F_CGC  <- file.path(IN, "cancer_gene_census.v97.genes.tsv")
F_TMB  <- file.path(IN, "6633_2729_3248-filtered_mutations_matched_allTum_keep_KEEPONLY_mutations_per_MB.tsv")
for (f in c(F_MAF, F_CLIN, F_CGC, F_TMB)) if (!file.exists(f)) stop("missing: ", f)

TOP_N <- 20

# ---------------------------------------------------------------- colours
# Named vectors: a category keeps its colour whatever order the levels
# come in. A level with no entry here is drawn white and reads as missing.
VC_COL <- c(
  Missense_Mutation      = "#1B7F7C",
  Nonsense_Mutation      = "#F39B7F",
  In_Frame_Del           = "#7E6148",
  Frame_Shift_Ins        = "#3C5488",
  Frame_Shift_Del        = "#91D1C2",
  Multi_Hit              = "#8491B4",
  In_Frame_Ins           = "#B09C85",
  Splice_Site            = "#DC0000",
  Nonstop_Mutation       = "#4DBBD5",
  Translation_Start_Site = "#E64B35"
)

TUMOUR_COL <- c(
  Primary               = "#2E3192",
  Lymph_node_metastasis = "#6A6BAF",
  Recurrence            = "#D9D9EC",
  Not_applicable        = "grey80"
)

SITE_COL <- c(
  Plantar                   = "#E8452B",
  Subungual                 = "#8B1A1A",
  Lymph_node                = "#F9DBD5",
  `in-transit_metastasis`   = "#F4A582",
  Not_applicable            = "grey80"
)

STAGE_COL <- c(
  Ib  = "#40E0C0", IIa = "#2E8B57", IIb = "#A8E6DC",
  IIc = "#A8D5A2", III = "#2A9D8F", IV  = "#1D5E5A",
  Not_applicable = "grey80"
)

# ------------------------------------------------------- cohort + clinical
clin_raw <- read_csv(F_CLIN, show_col_types = FALSE) |>
  select(where(~ !all(is.na(.))))

NEEDED <- c("ID_Sanger", "ID_lab", "Anatomic_site", "Tumor_type", "Clinical_Stage")
if (!all(NEEDED %in% names(clin_raw)))
  stop("missing columns: ", paste(setdiff(NEEDED, names(clin_raw)), collapse = ", "),
       "\navailable: ", paste(names(clin_raw), collapse = ", "))

clin <- clin_raw |>
  transmute(
    PD     = str_squish(ID_Sanger),
    AM     = str_squish(ID_lab),
    Tumor_type     = str_squish(Tumor_type)     |> str_replace_all(" ", "_"),
    Anatomic_site  = str_squish(Anatomic_site)  |> str_replace_all(" ", "_"),
    Clinical_Stage = str_squish(Clinical_Stage)
  ) |>
  mutate(across(c(Tumor_type, Anatomic_site, Clinical_Stage),
                \(x) coalesce(na_if(x, ""), "Not_applicable")))

# PD70976a is labelled AM057a in the source file, but AM057a is PD70966a,
# a different donor. The correct label is AM066a, consistent with
# PD70976c = AM066b in the copy-number metadata.
FIX <- c(PD70976a = "AM066a")
stopifnot(all(names(FIX) %in% clin$PD))
clin <- clin |> mutate(AM = coalesce(unname(FIX[PD]), AM))
message("relabelled: ", paste(names(FIX), "->", FIX, collapse = ", "))

stopifnot(!anyDuplicated(clin$PD), !anyDuplicated(clin$AM))
message(nrow(clin), " samples, ", n_distinct(clin$AM), " AM ids")

# One Sanger donor must map to one AM patient number.
donor <- sub("[a-z]+$", "", clin$PD)
chk <- tapply(sub("[a-z]*$", "", clin$AM), donor, \(x) length(unique(x)))
if (any(chk > 1))
  stop("Sanger donor mapping to more than one AM number: ",
       paste(names(chk)[chk > 1], collapse = ", "))

message("Tumor_type:");     print(table(clin$Tumor_type))
message("Anatomic_site:");  print(table(clin$Anatomic_site))
message("Clinical_Stage:"); print(table(clin$Clinical_Stage))

# All three variables are reported at once so the palettes can be
# completed in one pass.
PALS <- list(Tumor_type = TUMOUR_COL, Anatomic_site = SITE_COL,
             Clinical_Stage = STAGE_COL)
miss <- lapply(names(PALS), \(v) setdiff(unique(clin[[v]]), names(PALS[[v]])))
names(miss) <- names(PALS)
miss <- miss[lengths(miss) > 0]
if (length(miss))
  stop("levels with no colour:\n",
       paste0("  ", names(miss), ": ", sapply(miss, paste, collapse = ", "),
              collapse = "\n"),
       "\nadd them to the palettes at the top")

# ------------------------------------------------------------------- MAF
maf_df <- as.data.frame(fread(F_MAF, sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE))

sem_mut <- setdiff(clin$PD, unique(maf_df$Tumor_Sample_Barcode))
if (length(sem_mut))
  message("no mutation in the MAF (kept, drawn empty): ",
          paste(sem_mut, collapse = ", "))

# Both tables must carry the same identifiers: renaming only the MAF makes
# read.maf drop the clinical annotation for every sample.
am_of <- setNames(clin$AM, clin$PD)

maf_df <- maf_df |>
  filter(Tumor_Sample_Barcode %in% clin$PD) |>
  mutate(Tumor_Sample_Barcode = unname(am_of[Tumor_Sample_Barcode]))

clin_maf <- clin |>
  transmute(Tumor_Sample_Barcode = AM, Tumor_type, Anatomic_site, Clinical_Stage)

maf <- read.maf(maf = maf_df, clinicalData = clin_maf, verbose = FALSE)
message("samples in the MAF object: ",
        as.numeric(maf@summary[ID == "Samples", summary]), " / ", nrow(clin))

# -------------------------------------------------------------- TMB bar
# CAVEAT: this file is derived from the "keep" MAF, while the oncoplot
# cells come from "keepPA". The filters differ (keep retains synonymous and
# non-coding calls), so the top bar and the cells do not share a
# denominator. State this in the figure legend, or regenerate the TMB from
# the same MAF with maftools::tmb(maf, captureSize = <kit size in Mb>).
tmb_raw <- fread(F_TMB, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (!grepl("^PD", tmb_raw[[1]][1])) tmb_raw <- tmb_raw[-1, ]

i_s <- which(vapply(tmb_raw, \(x) mean(grepl("^PD", trimws(as.character(x)))) > .8,
                    logical(1)))[1]
i_v <- which(vapply(tmb_raw, \(x) {
  v <- suppressWarnings(as.numeric(as.character(x))); mean(!is.na(v)) > .8
}, logical(1)))[1]
if (is.na(i_s) || is.na(i_v))
  stop("could not identify the sample and value columns in ", basename(F_TMB))

tmb <- data.frame(
  PD  = str_squish(as.character(tmb_raw[[i_s]])),
  TMB = suppressWarnings(as.numeric(as.character(tmb_raw[[i_v]])))
) |>
  filter(PD %in% clin$PD) |>
  transmute(Tumor_Sample_Barcode = unname(am_of[PD]), TMB)

# A sample with no TMB value would get no bar, which reads as zero
# mutations rather than missing data.
falta_tmb <- setdiff(clin$AM, tmb$Tumor_Sample_Barcode)
if (length(falta_tmb))
  stop("no TMB value for: ", paste(falta_tmb, collapse = ", "))
stopifnot(!anyDuplicated(tmb$Tumor_Sample_Barcode), !any(is.na(tmb$TMB)))
message("TMB range: ", paste(round(range(tmb$TMB), 2), collapse = " - "),
        " mut/Mb")

# ------------------------------------------------------------------ genes
cgc <- read.csv(F_CGC, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
i   <- grep("^gene", tolower(names(cgc)))[1]
stopifnot(!is.na(i))
CGC <- unique(str_squish(cgc[[i]]))

gsum   <- getGeneSummary(maf)
in_cgc <- gsum$Hugo_Symbol[gsum$Hugo_Symbol %in% CGC]
genes  <- head(in_cgc, TOP_N)
message("genes: ", paste(genes, collapse = ", "))
writeLines(in_cgc, file.path(OUT, "COSMIC_matched_samples_BRA_AM.txt"))

# --------------------------------------------------------------- figure 1
fig <- file.path(OUT, paste0("oncoplot_one_per_patient_top", TOP_N, ".pdf"))
pdf(fig, width = 18, height = 16)
oncoplot(
  maf                     = maf,
  genes                   = genes,
  colors                  = VC_COL,
  topBarData              = tmb,     # mutations/Mb instead of raw counts
  clinicalFeatures        = c("Tumor_type", "Clinical_Stage", "Anatomic_site"),
  annotationColor         = list(Tumor_type     = TUMOUR_COL,
                                 Clinical_Stage = STAGE_COL,
                                 Anatomic_site  = SITE_COL),
  sortByMutation          = TRUE,
  removeNonMutated        = FALSE,   # keep samples with none of these genes mutated
  showTumorSampleBarcodes = TRUE,
  SampleNamefontSize      = 1.2,
  annotationFontSize      = 1.6,
  legendFontSize          = 1.6,
  fontSize                = 1,
  barcode_mar             = 7,
  gene_mar                = 8,
  draw_titv               = FALSE,
  bgCol                   = "#E5E5E5"
)
dev.off()
message("figure: ", fig)

# --------------------------------------------------------------- figure 2
# Same cohort, pre-specified gene list, genes mutated in >= MIN_CASES
# samples.
MIN_CASES <- 2

# MUC16 is deliberately absent: it is the second entry of the FLAGS list
# (Shyr et al., BMC Med Genomics 2014), 14,508 aa, and accumulates rare
# non-synonymous variants in healthy population exomes in proportion to
# coding length. Excluded a priori.
GOI <- unique(c(
  "KIT", "NRAS", "BRAF", "CDKN2A", "CSMD3", "NF1", "ANK1", "ASXL1", "CBL",
  "CNTNAP2", "CUX1", "FAM47C", "FAT1", "HRAS", "IL7R", "KIF5B", "KRAS",
  "LIFR", "LRP1B", "LZTR1"
))
stopifnot(length(GOI) == 20)

freq_all <- gsum |> as.data.frame() |> select(Hugo_Symbol, MutatedSamples)

keep <- freq_all |>
  filter(Hugo_Symbol %in% GOI, MutatedSamples >= MIN_CASES) |>
  arrange(desc(MutatedSamples), Hugo_Symbol)
goi_keep <- keep$Hugo_Symbol

message("\nMIN_CASES = ", MIN_CASES, " of ", nrow(clin), " samples")
message("  genes shown: ", length(goi_keep), "/", length(GOI), " -> ",
        paste(goi_keep, collapse = ", "))

# Pre-specified genes that do not make the figure: absence is a result,
# and a reader cannot tell "not mutated" from "not looked at".
few  <- freq_all |> filter(Hugo_Symbol %in% GOI, MutatedSamples < MIN_CASES)
none <- setdiff(GOI, freq_all$Hugo_Symbol)
if (nrow(few))
  message("  from the list, only 1 case: ",
          paste(few$Hugo_Symbol, collapse = ", "))
if (length(none))
  message("  from the list, no mutation at all: ",
          paste(none, collapse = ", "))

# Quantify what the a priori FLAGS exclusion costs, so it is not silent.
muc <- freq_all |> filter(Hugo_Symbol == "MUC16")
if (nrow(muc))
  message("  MUC16 excluded a priori (FLAGS), had ", muc$MutatedSamples,
          " cases")

stopifnot(length(goi_keep) > 0)

write_csv(freq_all |> arrange(desc(MutatedSamples)),
          file.path(OUT, "gene_frequency_all.csv"))

fig2 <- file.path(OUT, paste0("oncoplot_genes_of_interest_min", MIN_CASES, ".pdf"))
pdf(fig2, width = 18, height = 16)
oncoplot(
  maf                     = maf,
  genes                   = goi_keep,
  colors                  = VC_COL,
  topBarData              = tmb,
  clinicalFeatures        = c("Tumor_type", "Clinical_Stage", "Anatomic_site"),
  annotationColor         = list(Tumor_type     = TUMOUR_COL,
                                 Clinical_Stage = STAGE_COL,
                                 Anatomic_site  = SITE_COL),
  sortByMutation          = TRUE,
  removeNonMutated        = FALSE,
  showTumorSampleBarcodes = TRUE,
  SampleNamefontSize      = 1.2,
  annotationFontSize      = 1.6,
  legendFontSize          = 1.6,
  fontSize                = 1,
  barcode_mar             = 7,
  gene_mar                = 8,
  draw_titv               = FALSE,
  bgCol                   = "#E5E5E5"
)
dev.off()
message("figure: ", fig2)
