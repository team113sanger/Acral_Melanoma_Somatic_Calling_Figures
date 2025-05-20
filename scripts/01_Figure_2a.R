### ou diretório acessível em todos os nós
mylib <- "/data04/projects04/PatriciaPossik/pdx_la_sanger/lib/Rpackages/4.4.1"
libdir <- "/data04/tools/R/Rpackages/4.4.1"

### comando para efetivar a seleção de libPaths
.libPaths(c(mylib,libdir))
.libPaths()
### set CRAN Mirror https://cran-r.c3sl.ufpr.br/
options(repos = c("CRAN" = "https://cran-r.c3sl.ufpr.br/"))

#setting directory
# Start by reading the libraries that are required
# 0. Set the variables for the Project -------
# Check for requirements
#PREPARE the environment 

##Dependencies for this project are recorded and managed with renv
## To setup your Rstudio environment you might need to run:
#renv::init() # To initialise revn
renv::restore()
#renv::init() # To initialise revn
#renv::install("glue") # To add a package from CRAN
#renv::install("BiocManager") # To add a package from Bioconductor
#renv::install("bioC::ggplot2") # To add a package from Bioconductor
#renv::snapshot() # To record any packages you've used and their versions and sources into the renv.lock file
#renv::restore() # To rebuild an environment from the renv.lockfile

library(renv) # To retain package versions
library(here)
library(RColorBrewer)
library(writexl)
library(data.table)
library(maftools)
library(ggplot2)
renv::snapshot()
# setting directory  

# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/01_Figure_2a.R")

# 0 Set the directroies to work ----
wdir<-here()
pdx_dir<-file.path(wdir, "PDX")
tum_dir<-file.path(wdir, "Tumour")
cosmic_genesv97<-read.csv(file =file.path(wdir, "COSMIC","cancer_gene_census.v97.genes.tsv" ), header = T, sep = "\t", stringsAsFactors = F )


# 1 Merge and plot the data for Combined all_samples  matched & unmatched ----
prefix<-"all_samples"
suffix<-"all_samples"
comball_dir<-file.path(wdir, prefix)

#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, "all_samples", "6633_3248-filtered_mutations_all_allTum_keep.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, "all_samples", "6633_2729-filtered_mutations_all_allTum_keep.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)

#Read the PDX and Tumour TMBs mafs to merge them
pdx_tmb<- as.data.frame(data.table::fread(input = file.path(pdx_dir, "all_samples", "mutations_per_Mb.tsv"), 
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_tmb<- as.data.frame(data.table::fread(input = file.path(tum_dir, "all_samples", "mutations_per_Mb.tsv"),
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_tmb<- rbind(pdx_tmb, tum_tmb)
colnames(comb_tmb)<-c("PDID", "mutations_per_MB")

################### oncoplot to paper using only one sample per patient

wdir<-here()
tum_dir<-file.path(wdir, "patient_sample")
# You have to downlaod the COSMIC v97 (see COSMIC folder)
cosmic_genesv97<-read.csv(file =file.path(wdir, "COSMIC","cancer_gene_census.v97.genes.tsv" ), header = T, sep = "\t", stringsAsFactors = F )


# 1 Merge and plot the data for all_samples matched & unmatched - Patient only ----
prefix<-"patient_samples"
suffix<-"patient_samples"
tumour_only_dir<-file.path(wdir, prefix)

#Read the Tumour TMBs mafs to merge them
#Create the directory for the files 
dir.create(tumour_only_dir, recursive = T)
#Save the MAF files
patient_maf_fname<-file.path(tumour_only_dir, "6633_2729_3248-filtered_mutations_all_allTum_keep")
data.table::fwrite(tum_maf, paste0(patient_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(tum_maf, paste0(patient_maf_fname, ".xlsx"),col_names = T, format_headers = T )
#Save the TMB 
writexl::write_xlsx(tum_tmb, paste0(patient_maf_fname, "mutations_per_MB.xlsx"),col_names = T, format_headers = T )

# 2 Read the Combined patient_samples matched & unmatched maf file for oncoplots ----
#Try reading the maf
patient_maf<-read.maf(maf=paste0(patient_maf_fname, ".maf"), clinicalData = file.path(wdir,"data","METADATA_all_clinical_reviewed_dec_2024.tsv") )


#get rid of samples of interest in MAF files, Tumor_Sample_Barcode (tirando as unmacthed)
filterd_maf_new <- filterMaf(maf = patient_maf, tsb = c("PD53335e", "PD53336f", "PD53340a", "PD53344f", "PD53356c", "PD53360a"))
#get rid of samples of interest in MAF files, Tumor_Sample_Barcode (tirando as amostras repetidas dos pacientes,
#para ficar one sample per patient)
filterd_maf_new <- filterMaf(maf = filterd_maf_new, tsb = c("PD53332e", "PD53337c", "PD53343i", "PD53349e", 
                                                            "PD53350e", "PD53352e", "PD53357d", "PD53357e"))

tcgaCompare(maf = filterd_maf_new, cohortName = c("Brazilian_AM_cohort"),
            logscale= TRUE, capture_size = 48.225157, rm_zero = FALSE)

laml.titv = titv(maf = filterd_maf_new, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

data.table::fwrite(tum_maf, paste0(filterd_maf_new.maf, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
# Salvar o MAF
maftools::write.maf(maf = filterd_maf_new, file = output_file)
#Set the colours for the Tumour_site
pcol<-c("#295888", "#FCCDE5", "#FFCC99")
#pcol<-c("#295884", "#D20C60")
names(pcol)<- unique(filterd_maf_new@clinical.data$Tumor_Type)
#Set the colours for the  clinical_stage
clini_stage<- c("#CCEBC5", "#FFED6F", "#8DD3C7","#BC80BD","#BEBADA","#A6761D")  
names(clini_stage)<- unique(filterd_maf_new@clinical.data$Clinical_Stage )
#Set the colour for anatomic site
anatomic<- c("#1B9E77", "#B3DE69",  "#D20C60") 
names(anatomic)<- unique(filterd_maf_new@clinical.data$Anatomic_Site )

clin_col<-list(Tumor_Type=pcol,
               Clinical_Stage=clini_stage,
               Anatomic_Site=anatomic)


# The results is in patient_samples folder
# Oncoplot showing ALLS COSMIC genes from the CANCER gene census only
topn<-20
# THIS IS how the overlaps should be taken 
observed_mut_canc_genes <-filterd_maf_new@gene.summary$Hugo_Symbol[filterd_maf_new@gene.summary$Hugo_Symbol %in% cosmic_genesv97$Gene.Symbol]
writeLines(observed_mut_canc_genes, "COSMIC_matched_samples_BRA_AM.txt")

fig_fnamepdf<-file.path(tumour_only_dir, paste0("FINAL3_oncoplot_tumors", topn, "_COSMICv97cgenes_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)

oncoplot(filterd_maf_new,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not
         genes = observed_mut_canc_genes[1:topn], # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         SampleNamefontSize = 0.9, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         annotationFontSize = 1.6,
         fontSize = 1,
         legendFontSize = 1.6, 
         barcode_mar = 7,
         gene_mar=8,
         clinicalFeatures = c('Tumor_Type','Clinical_Stage','Anatomic_Site'),
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = TRUE #Draw the SNVs by substitution type 
)
dev.off() #plot saved in Patient_Samples folder

#to see th type of mutations
library(dplyr)
genes_interesse <- c("BRAF", "KRAS", "HRAS", "NRAS", "NF1", "KIT")
variantes <- filterd_maf_new@data %>%
  filter(Hugo_Symbol %in% genes_interesse) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short)

write.csv(variantes, file = "variants_drivers_per_sample_matched_one_sample_per_patient.csv", row.names = FALSE)




#######
fig_fnamepdf<-file.path(tumour_only_dir, paste0("FINAL2_oncoplot_tumors", topn, "_COSMICv97cgenes_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 20, height = 8)
oncoplot(filterd_maf_new,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not
         genes = c("NRAS","HRAS","KRAS","KIT","NF1"), # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 0.9, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         annotationFontSize = 1.6,
         fontSize = 1,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=8,
         clinicalFeatures = c('Tumor_Type','Clinical_Stage','Anatomic_Site'),
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()

oncoplot(maf = filterd_maf_new, pathways = "smgbp", gene_mar = 8, fontSize = 0.8, topPathways = 5)
plotmafSummary(maf = filterd_maf_new, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#exclusive/co-occurance event analysis on top 25 mutated genes. 
somaticInteractions(maf = filterd_maf_new, top = 25, colNC = 9,
                    nShiftSymbols = 5,
                    sigSymbolsSize = 1,
                    sigSymbolsFontSize = 0.6, 
                    countsFontSize = 0.4, pvalue = c(0.05, 0.1))

#Detecting cancer driver genes based on positional clustering
laml.sig = oncodrive(maf = filterd_maf_new, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)
#We can plot the results using plotOncodrive.

plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

dgi = drugInteractions(maf = filterd_maf_new, fontSize = 0.8, top = 10)

# Defina o nome do arquivo e as dimensões desejadas
png("drug_interactions.png", width = 800, height = 500) # largura ajustada para 800 pixels e altura para 1200 pixels
# Gere o gráfico
dgi = drugInteractions(maf = filterd_maf_new, fontSize = 1.5, top = 15)
# Salve o gráfico
dev.off()


