#Author: Annie Squiavinato
#Oct-2026
#_______________________________________________________________________________
# Patient tumour vs matched PDX-X1.
#   Fig 3 panel a  coBarplot, driver groups (KIT / NF1 / BRAF / RAS collapsed)
#   Fig 3 panel c  coBarplot, the remaining genes of interest
#   plus     SNP/del/ins burden correlation patient vs PDX-X1
# The CNApp score panels live in their own script; the shared lesion palette
# is built here and exported for it.
#
# INPUTS  
#   6633_2729_3248-filtered_mutations_matched_allTum_keepPA_KEEPONLY.maf
#       somatic mutations, matched tumour-normal, protein-altering only
#   6633_METADATA_table_2729_3248_all_data.tsv
#       WES metadata: Tumor_Sample_Barcode, Supplier_Sample_ID, sample_type,
#       PDX_passage, SAMPLE_USE
#   cancer_gene_census.v97.genes.tsv
#       COSMIC Cancer Gene Census v97, https://cancer.sanger.ac.uk/census
# or use COSMIC_matched_samples_BRA_AM.txt file (only genes used on this figure)
#   mapk_pathway_genes.txt
#       written on first run from MSigDB (MAPK_SET), then read from disk
#  /data/metadata_copy_number_data_filtered_by_cellularity_all_samples.csv
#       optional; only widens the lesion palette so it matches the CNApp script
#
# OUTPUTS (under SVC/results)
#   coBarplot_driver_groups.pdf
#   coBarplot_GOI_remaining.pdf
#   Correlation_mutation_burden_tumour_vs_pdx.pdf
#   lesion_colours.tsv                  palette shared with the CNApp script
#   tumour_pdx_pairs_MAF.tsv            the tumour-PDX pairs analysed
#   gene_frequency_CGC.tsv              all CGC genes, both cohorts
#   gene_frequency_MAPK_CGC.tsv         CGC genes in the MAPK set
#   gene_frequency_driver_groups.tsv    panel a
#   gene_frequency_GOI_remaining.tsv    panel c
#   mutation_burden_tumour_vs_pdx.tsv   counts per variant type, per pair
#
# Post-processing in Affinity for the submitted figures.

# Package versions and the library path are managed with renv: run
# renv::restore() once to rebuild the environment recorded in renv.lock.
#_______________________________________________________________________________

library(maftools)
library(dplyr)      # after maftools, so dplyr verbs win
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggsci)
library(msigdbr)

# Run from analysis/ — the directory that holds renv.lock,
# data/, scripts/ and results/. Run renv::restore() once before the first use.
IN  <- "analysis/data"
OUT <- "analysis/results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the analysis/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MIN_CASES <- 2
MAPK_SET  <- "REACTOME_MAPK1_MAPK3_SIGNALING"

F_MAF  <- file.path(IN, "6633_2729_3248-filtered_mutations_matched_allTum_keepPA_KEEPONLY.maf")
F_WES  <- file.path(IN, "6633_METADATA_table_2729_3248_all_data.tsv")
F_CGC  <- file.path(IN, "cancer_gene_census.v97.genes.tsv")
F_MAPK <- file.path(IN, "mapk_pathway_genes.txt")
F_CNA  <- file.path(IN, "metadata_copy_number_data_filtered_by_cellularity_all_samples.csv")

# F_MAPK is written on the first run from MSigDB; F_CNA is optional and only
# widens the lesion palette, so neither is required to exist here.
for (f in c(F_MAF, F_WES, F_CGC))
  if (!file.exists(f)) stop("missing input: ", f)

if (!file.exists(F_CNA))
  message("optional input not found, palette will be narrower: ", F_CNA)

# The 20 pre-specified acral melanoma genes, the same list as the oncoplot.
# Panel b is this set minus the driver groups of panel a, derived rather than
# typed out, so the two figures cannot drift apart.
GOI <- c("KIT", "NRAS", "BRAF", "CDKN2A", "CSMD3", "NF1", "ANK1", "ASXL1",
         "CBL", "CNTNAP2", "CUX1", "FAM47C", "FAT1", "HRAS", "IL7R",
         "KIF5B", "KRAS", "LIFR", "LRP1B", "LZTR1")
stopifnot(length(GOI) == 20, !anyDuplicated(GOI))

RAS_GENES <- c("NRAS", "KRAS", "HRAS")
GROUPS    <- c("KIT", "NF1", "BRAF", "RAS")

# Lesion labels never legitimately contain whitespace, and the source has it
# in four positions ("AM031-X1 _", "AM049_ X1", "AM0104a- X1",
# "AM0104a-X1- Cell"). The "_X1" / "-X1" separator is also inconsistent;
# without normalising it those pairs silently fail to join.
am_label  <- function(x) x |>
  str_remove_all("\\s") |>
  str_remove("_(gDNA|RNA|gRNA)_tumour$") |>
  str_replace_all("_X", "-X")
am_lesion <- function(x) str_remove(am_label(x), "-X\\d+.*$")

# ---------------------------------------------------------------- palette
# Union of both metadata files in numeric order, so a lesion keeps its colour
# across panels even though the cohorts differ.
lev <- c(
  if (file.exists(F_CNA)) am_lesion(read_csv(F_CNA, show_col_types = FALSE)$Lab_id),
  am_lesion(read_tsv(F_WES, show_col_types = FALSE)$Supplier_Sample_ID)
)
# [a-z]* matters: labels carry a letter suffix (AM015a, AM015b) and
# "^AM\\d+$" would drop every one of them, leaving those points uncoloured.
lev <- unique(lev[grepl("^AM\\d+[a-z]*$", lev)])
lev <- lev[order(as.integer(str_extract(lev, "\\d+")), str_extract(lev, "[a-z]*$"))]

# pal_igv("default")(51) returns only 49 distinct values (#CC9900FF and
# #99CC00FF each appear twice), so several palettes are pooled.
pool <- unique(c(pal_igv("default")(51), pal_ucscgb("default")(26),
                 pal_d3("category20")(20), pal_simpsons("springfield")(16),
                 pal_futurama("planetexpress")(12)))
if (length(pool) < length(lev))
  stop(length(lev), " lesions but only ", length(pool), " distinct colours")

PAL <- setNames(pool[seq_along(lev)], lev)
stopifnot(!anyDuplicated(PAL), !any(is.na(PAL)))
write_tsv(tibble(lesion = lev, colour = unname(PAL)),
          file.path(OUT, "lesion_colours.tsv"))
message(length(PAL), " lesions in the palette")

# ------------------------------------------------------------------ pairs
wes <- read_tsv(F_WES, show_col_types = FALSE)
stopifnot(all(c("Tumor_Sample_Barcode", "Supplier_Sample_ID", "sample_type",
                "PDX_passage", "SAMPLE_USE") %in% names(wes)))

maf <- read.maf(F_MAF, verbose = FALSE)
tsb <- as.character(getSampleSummary(maf)$Tumor_Sample_Barcode)
message("samples in MAF: ", length(tsb))

# SAMPLE_USE is the curated exclusion list and the single source of truth;
# a hard-coded vector drifts away from it.
rm_use <- wes |> filter(str_detect(coalesce(SAMPLE_USE, ""), "^Remove"),
                        Tumor_Sample_Barcode %in% tsb)
message("dropped by SAMPLE_USE: ", nrow(rm_use), " -> ",
        paste(rm_use$Tumor_Sample_Barcode, collapse = ", "))

m <- wes |>
  filter(!str_detect(coalesce(SAMPLE_USE, ""), "^Remove")) |>
  transmute(Sample = str_squish(as.character(Tumor_Sample_Barcode)),
            AM     = am_label(Supplier_Sample_ID),
            type   = str_squish(sample_type),
            pass   = str_squish(PDX_passage)) |>
  filter(nzchar(Sample))

# PD53368a is the patient tumour AM049; its label carries a spurious "_ X1".
# Its PDX (PD53368c) was never sequenced, so AM049 correctly has no pair.
FIX <- tibble::tribble(~Sample, ~AM, "PD53368a", "AM049")
stopifnot(all(FIX$Sample %in% m$Sample))
m <- rows_update(m, FIX, by = "Sample")

# Pairing is on the lesion label, not the patient number: AM022a and AM022b
# are different lesions of the same donor, each with its own PDX.
tum <- m |> filter(type == "Tumour") |> transmute(AM_tum = AM, PD_tum = Sample)
pdx <- m |> filter(type == "PDX", pass == "P1") |>
  transmute(AM_tum = str_remove(AM, "-X1$"), AM_pdx = AM, PD_pdx = Sample)

# A duplicated label inflates the join when it matches on the other side and
# drops the pair silently when it does not.
dup <- c(tum$AM_tum[duplicated(tum$AM_tum)], pdx$AM_tum[duplicated(pdx$AM_tum)])
if (length(dup)) stop("duplicated lesion labels: ",
                      paste(unique(dup), collapse = ", "))

pairs <- inner_join(tum, pdx, by = "AM_tum") |>
  filter(PD_tum %in% tsb, PD_pdx %in% tsb) |>
  mutate(lesion = AM_tum, patient = str_extract(AM_tum, "^AM\\d+"))

stopifnot(nrow(pairs) > 0, all(pairs$lesion %in% names(PAL)),
          !anyDuplicated(pairs$PD_tum), !anyDuplicated(pairs$PD_pdx))
N_PAIRS <- nrow(pairs)
message(N_PAIRS, " pairs, ", n_distinct(pairs$patient), " patients")
rep2 <- count(pairs, patient) |> filter(n > 1)
if (nrow(rep2)) message("two lesions: ", paste(rep2$patient, collapse = ", "))
write_tsv(pairs, file.path(OUT, "tumour_pdx_pairs_MAF.tsv"))

m_tum <- subsetMaf(maf, tsb = pairs$PD_tum, mafObj = TRUE)
m_pdx <- subsetMaf(maf, tsb = pairs$PD_pdx, mafObj = TRUE)
n_tum <- as.numeric(m_tum@summary[ID == "Samples", summary])
n_pdx <- as.numeric(m_pdx@summary[ID == "Samples", summary])
if (is.na(n_tum) || n_tum != n_pdx || n_tum != N_PAIRS)
  warning("cohort sizes do not match the number of pairs")

# ------------------------------------------------------------- frequency
# Ranking over both cohorts combined; ranking on patients alone would hide a
# gene mutated only in the PDXs.
gs <- \(x, nm) getGeneSummary(x) |> as.data.frame() |>
  select(Hugo_Symbol, !!nm := MutatedSamples)

freq <- full_join(gs(m_tum, "n_tum"), gs(m_pdx, "n_pdx"), by = "Hugo_Symbol") |>
  mutate(across(c(n_tum, n_pdx), \(x) replace_na(x, 0)), total = n_tum + n_pdx) |>
  arrange(desc(total), desc(n_tum), Hugo_Symbol)

# `total` only orders the genes: the same lesion is counted once per cohort,
# so it is not a number of cases. The percentages are per cohort.
add_pct <- \(d) d |>
  mutate(pct_tum = round(100 * n_tum / N_PAIRS, 1),
         pct_pdx = round(100 * n_pdx / N_PAIRS, 1),
         n_pairs = N_PAIRS) |>
  arrange(desc(total), desc(n_tum), Hugo_Symbol)

# ------------------------------------------------------------------- CGC
# Full Census, not the BRA-filtered list: that one was derived from these
# same samples, so ranking within it would be circular.
cgc  <- read_tsv(F_CGC, show_col_types = FALSE)
norm <- gsub("[^a-z]", "", tolower(names(cgc)))
i    <- na.omit(match(c("genesymbol","hugosymbol","genename","symbol","gene"), norm))[1]
stopifnot(!is.na(i))
CGC  <- unique(str_squish(as.character(cgc[[i]])))
CGC  <- CGC[nzchar(CGC)]
stopifnot(length(CGC) > 100)     # v97 has ~740; far fewer means the wrong column

freq_cgc <- freq |> filter(Hugo_Symbol %in% CGC) |> add_pct()
write_tsv(freq_cgc, file.path(OUT, "gene_frequency_CGC.tsv"))

fora <- setdiff(GOI, CGC)
if (length(fora)) message("genes of interest NOT in the CGC: ",
                          paste(fora, collapse = ", "))

# ----------------------------------------------------------- MAPK pathway
# Supplementary table only, no figure. Reactome ERK1/2 arm rather than the
# whole MAPK family, which also carries p38 and JNK. Written to file on
# first run so the gene set is versioned with the analysis: MSigDB content
# changes between releases.
if (!file.exists(F_MAPK)) {
  g <- msigdbr(species = "Homo sapiens", collection = "C2",
               subcollection = "CP:REACTOME") |>
    filter(gs_name == MAPK_SET) |>
    pull(gene_symbol) |> unique() |> sort()
  stopifnot(length(g) > 0)
  writeLines(g, F_MAPK)
  message("wrote ", length(g), " genes from ", MAPK_SET, " to ", F_MAPK)
}
MAPK <- readLines(F_MAPK, warn = FALSE) |> str_squish() |> unique()
MAPK <- MAPK[nzchar(MAPK)]
message("MAPK genes: ", length(MAPK), " from ", basename(F_MAPK))

freq_mapk <- freq_cgc |> filter(Hugo_Symbol %in% MAPK)
write_tsv(freq_mapk, file.path(OUT, "gene_frequency_MAPK_CGC.tsv"))

# --------------------------------------------------------------- colours
# Multi_Hit is created by maftools at plot time for samples with more than
# one mutation in a gene; leaving it out of the colours draws it white.
vc <- sort(unique(c(as.character(m_tum@data$Variant_Classification),
                    as.character(m_pdx@data$Variant_Classification))))
vc <- unique(c(vc, "Multi_Hit"))
stopifnot(length(vc) <= 10)
VC_COL <- setNames(pal_npg("nrc")(10)[seq_along(vc)], vc)
VC_COL["Multi_Hit"] <- "grey40"

# Bar thickness is (plot area) / (n genes), and coBarplot has a fixed
# overhead (title, axis, legend) that does not scale with the number of
# rows. Height must be base + k * n for the bars to come out the same
# thickness in both panels, which is what joining them requires.
H_BASE <- 2.6    # title + axis + legend
H_ROW  <- 0.40   # per gene

co <- function(g, file, m1 = m_tum, m2 = m_pdx) {
  h <- H_BASE + H_ROW * length(g)
  pdf(file, width = 6, height = h)
  par(oma = c(0, 0, 1.5, 0))        # coBarplot sets its own mar via layout
  coBarplot(m1, m2, m1Name = "patients", m2Name = "AM-PDXs",
            genes = rev(g),         # coBarplot draws the first gene at the bottom
            colors = VC_COL, geneSize = 0.6, showPct = TRUE)
  dev.off()
  message("figure: ", file, "  (", length(g), " genes, height ",
          round(h, 2), ")")
}

# ------------------------------------------------- panel a: driver groups
# Collapsing on Hugo_Symbol before read.maf makes maftools count the group
# as one gene. A sample carrying both NRAS and KRAS then has two mutations
# in RAS and is drawn as Multi_Hit.
d <- as.data.frame(maf@data) |>
  mutate(Hugo_Symbol = ifelse(Hugo_Symbol %in% RAS_GENES, "RAS", Hugo_Symbol))

maf_g   <- read.maf(d, verbose = FALSE)
m_tum_g <- subsetMaf(maf_g, tsb = pairs$PD_tum, mafObj = TRUE)
m_pdx_g <- subsetMaf(maf_g, tsb = pairs$PD_pdx, mafObj = TRUE)

# Collapsing merges genes, not samples: the counts must not move.
stopifnot(as.numeric(m_tum_g@summary[ID == "Samples", summary]) == N_PAIRS,
          as.numeric(m_pdx_g@summary[ID == "Samples", summary]) == N_PAIRS)

freq_g <- full_join(gs(m_tum_g, "n_tum"), gs(m_pdx_g, "n_pdx"),
                    by = "Hugo_Symbol") |>
  mutate(across(c(n_tum, n_pdx), \(x) replace_na(x, 0)), total = n_tum + n_pdx) |>
  filter(Hugo_Symbol %in% GROUPS) |>
  add_pct()

print(as.data.frame(freq_g), row.names = FALSE)
write_tsv(freq_g, file.path(OUT, "gene_frequency_driver_groups.tsv"))

# RAS counts samples, not mutations: a sample mutated in two RAS genes counts
# once, so the collapsed row is smaller than the sum of the three original
# rows. Report which samples cause that gap.
multi_ras <- bind_rows(
  as.data.frame(m_tum@data)[, c("Hugo_Symbol", "Tumor_Sample_Barcode")],
  as.data.frame(m_pdx@data)[, c("Hugo_Symbol", "Tumor_Sample_Barcode")]) |>
  filter(Hugo_Symbol %in% RAS_GENES) |>
  distinct(Tumor_Sample_Barcode, Hugo_Symbol) |>
  count(Tumor_Sample_Barcode, name = "n_ras_genes") |>
  filter(n_ras_genes > 1)
if (nrow(multi_ras)) {
  message("samples mutated in more than one RAS gene:")
  print(as.data.frame(multi_ras), row.names = FALSE)
} else message("no sample carries more than one RAS gene")

co(freq_g$Hugo_Symbol, file.path(OUT, "coBarplot_driver_groups.pdf"),
   m1 = m_tum_g, m2 = m_pdx_g)

# ------------------------------------------- panel b: remaining GOI genes
# The 20 genes of interest minus the driver groups of panel a, derived so
# the two panels are complementary by construction.
REST <- setdiff(GOI, c("KIT", "NF1", "BRAF", RAS_GENES))
stopifnot(length(REST) == 14)

freq_rest <- freq |>
  filter(Hugo_Symbol %in% REST) |>
  add_pct()

# A gene of interest with no mutation is a result, but it would be drawn as
# an empty row, so it is reported and left out of the figure.
sem_mut <- setdiff(REST, freq_rest$Hugo_Symbol[freq_rest$total > 0])
if (length(sem_mut)) message("GOI with no mutation (not drawn): ",
                             paste(sem_mut, collapse = ", "))
baixo <- filter(freq_rest, total > 0, total < MIN_CASES)
if (nrow(baixo)) message("GOI with only 1 case (still drawn): ",
                         paste(baixo$Hugo_Symbol, collapse = ", "))

genes_rest <- freq_rest |> filter(total > 0) |> pull(Hugo_Symbol)
message("panel b: ", length(genes_rest), " genes -> ",
        paste(genes_rest, collapse = ", "))
stopifnot(length(genes_rest) > 0)
write_tsv(freq_rest, file.path(OUT, "gene_frequency_GOI_remaining.tsv"))
co(genes_rest, file.path(OUT, "coBarplot_GOI_remaining.pdf"))

# ---------------------------------------------------------- burden
# Breakdown by Variant_Type (SNP / INS / DEL), not by Variant_Classification:
# the axis is the kind of change, not its consequence on the protein.
cnt <- as.data.frame(maf@data) |>
  transmute(Sample = as.character(Tumor_Sample_Barcode),
            Variant_Type = as.character(Variant_Type)) |>
  count(Sample, Variant_Type) |>
  pivot_wider(names_from = Variant_Type, values_from = n, values_fill = 0) |>
  mutate(total = rowSums(across(-Sample)))

# getSampleSummary sums the same calls by a different route, so the two
# totals must agree; if they do not, one is counting something else.
chk <- getSampleSummary(maf) |> as.data.frame() |>
  transmute(Sample = as.character(Tumor_Sample_Barcode), ref = total) |>
  inner_join(select(cnt, Sample, total), by = "Sample")
stopifnot(nrow(chk) == nrow(cnt), all(chk$ref == chk$total))
message("variant types: ",
        paste(setdiff(names(cnt), c("Sample", "total")), collapse = ", "))

burden <- pairs |>
  left_join(rename_with(cnt, \(x) paste0(x, "_tum"), -Sample),
            by = c("PD_tum" = "Sample")) |>
  left_join(rename_with(cnt, \(x) paste0(x, "_pdx"), -Sample),
            by = c("PD_pdx" = "Sample")) |>
  rename(mut_tum = total_tum, mut_pdx = total_pdx) |>
  filter(!is.na(mut_tum), !is.na(mut_pdx))

fit <- lm(mut_pdx ~ mut_tum, data = burden)
ci  <- confint(fit)["mut_tum", ]
lim <- c(0, max(burden$mut_tum, burden$mut_pdx) * 1.05)
message("R2 = ", round(summary(fit)$r.squared, 3),
        "  slope = ", round(coef(fit)[2], 3),
        " (95% CI ", round(ci[1], 3), "-", round(ci[2], 3), ")",
        "  n = ", nrow(burden))

p <- ggplot(burden, aes(mut_tum, mut_pdx)) +
  # y = x reference: R2 alone cannot tell slope 1.0 from slope 1.3
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey35", fill = "grey80", alpha = 0.4) +
  geom_point(aes(colour = lesion), shape = 17, size = 4, alpha = 0.8) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("R² = %.2f\nslope = %.2f (%.2f-%.2f)\nn = %d",
                           summary(fit)$r.squared, coef(fit)[2], ci[1], ci[2],
                           nrow(burden)),
           hjust = 1.05, vjust = -0.25, size = 3.6, lineheight = 1.05) +
  scale_colour_manual(values = PAL) +
  coord_cartesian(xlim = lim, ylim = lim) +
  labs(x = "SNP/del/ins patient tumours", y = "SNP/del/ins AM-PDXs",
       colour = "Lesion") +
  theme_classic(base_size = 13) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3)))

print(p)
ggsave(file.path(OUT, "Correlation_mutation_burden_tumour_vs_pdx.pdf"),
       p, device = "pdf", width = 5.52, height = 5.5)

# One row per pair, counts per variant type on both sides plus the total.
burden_out <- burden |>
  select(patient, lesion, AM_tum, AM_pdx, PD_tum, PD_pdx,
         ends_with("_tum"), ends_with("_pdx")) |>
  rename(total_tum = mut_tum, total_pdx = mut_pdx)
write_tsv(burden_out, file.path(OUT, "mutation_burden_tumour_vs_pdx.tsv"))
message("\ndone")
