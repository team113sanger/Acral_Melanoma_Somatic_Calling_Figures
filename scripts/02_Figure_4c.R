
### Author: Annie Squiavinato  and Antonio
### date: "2025-04-25"  

### Directory for the new library 24-06-2024  
### or directory accessible on all nodes  
mylib <- "/data04/projects04/PatriciaPossik/pdx_la_sanger/lib/Rpackages/4.4.1"  
libdir <- "/data04/tools/R/Rpackages/4.4.1"  

### Command to apply the libPaths selection  
.libPaths(c(mylib, libdir))  
.libPaths()  

# setting directory  
setwd("~/pdx_la_sanger/results/INCA_res")  
### set CRAN Mirror https://cran-r.c3sl.ufpr.br/  
options(repos = c("CRAN" = "https://cran-r.c3sl.ufpr.br/"))  

### Directory for the new library 24-06-2024  
libdir <- "/scr/R/Rpackages/4.4.1"  

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

# setting directory  
setwd("~/pdx_la_sanger/results/INCA_res/Somatic_variant_calling_analysis/")  
# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/02_Figure_4c.R")

library(renv) # To retain package versions
library(here)
library(RColorBrewer)
library(writexl)
library(data.table)
library(maftools)
library(BiocManager)
renv::snapshot()

#To record any packages you've used and their versions and sources into the renv.lock file
renv::snapshot()
# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/05_Figure_4c.R")


# 0 Set the directroies to work ----
wdir<-here()
data_dir <- file.path(wdir, "data")
pdx_dir<-file.path(wdir, "PDX")
tum_dir<-file.path(wdir, "Tumour")
# You have to downlaod the COSMIC v97 (see COSMIC folder)
cosmic_genesv97<-read.csv(file =file.path(wdir, "COSMIC","cancer_gene_census.v97.genes.tsv" ), header = T, sep = "\t", stringsAsFactors = F )

# 3 Merge and plot the data for Combined all_samples  matched & unmatched ----
prefix<-"matched_samples"
suffix<-"matched_samples"
comball_dir<-file.path(wdir, prefix)

#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, prefix, "6633_3248-filtered_mutations_matched_allTum_keep.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, prefix, "6633_2729-filtered_mutations_matched_allTum_keep.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)

#Read the PDX and Tumour TMBs mafs to merge them
pdx_tmb<- as.data.frame(data.table::fread(input = file.path(pdx_dir, prefix, "mutations_per_Mb.tsv"), 
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_tmb<- as.data.frame(data.table::fread(input = file.path(tum_dir, prefix, "mutations_per_Mb.tsv"),
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_tmb<- rbind(pdx_tmb, tum_tmb)
colnames(comb_tmb)<-c("PDID", "mutations_per_MB")
#Create the directory for the files 
dir.create(comball_dir, recursive = T)
#Save the MAF files
comb_maf_fname<-file.path(comball_dir, "6633_2729_3248-filtered_mutations_matched_allTum_keep")
data.table::fwrite(comb_maf, paste0(comb_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(comb_maf, paste0(comb_maf_fname, ".xlsx"),col_names = T, format_headers = T )
#Save the TMB 
writexl::write_xlsx(comb_tmb, paste0(comb_maf_fname, "mutations_per_MB.xlsx"),col_names = T, format_headers = T )
#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, prefix, "6633_3248-filtered_mutations_matched_allTum_keepPA.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, prefix, "6633_2729-filtered_mutations_matched_allTum_keepPA.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)
#Save the MAF files
comb_maf_fname<-file.path(comball_dir, "6633_2729_3248-filtered_mutations_matched_allTum_keepPA")
data.table::fwrite(comb_maf, paste0(comb_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(comb_maf, paste0(comb_maf_fname, ".xlsx"),col_names = T, format_headers = T )

# 4 Read the Combined all_samples matched & unmatched maf file for oncoplots ----
#Try reading the maf
comb_maf<-read.maf(maf=paste0(comb_maf_fname, ".maf"), clinicalData = file.path(wdir, "data", "METADATA_all_clinical_reviewed_dec_2024.tsv") )


#get rid of samples of interest in MAF files, Tumor_Sample_Barcode (tirando as unmacthed)
tumor_new_maf <- subsetMaf(maf = comb_maf, tsb = c("PD53330c", "PD53332d", "PD53333f", "PD53343h", "PD53343i", 
                                                   "PD53347c", "PD53349d","PD53349e", "PD53350d", "PD53350e", "PD53352d", 
                                                   "PD53352e", "PD53355c", "PD53357d", "PD53359c", "PD53364c"
))
                                                   
                                                   

#get rid of samples of interest in MAF files, Tumor_Sample_Barcode (tirando as amostras de PDX)
#tumor_new_maf <- filterMaf(maf = tumor_new_maf, tsb = c("PD53333a_hum", "PD53333c_hum", "PD53333d_hum", "PD53333e_hum", "PD53335a_hum",
                                                   #"PD53335c_hum", "PD53335d_hum", "PD53336a_hum", "PD53336c_hum", "PD53336d_hum",
                                                  # "PD53336e_hum", "PD53343a_hum", "PD53343d_hum", "PD53343e_hum", "PD53343c_hum",
                                                  # "PD53343f_hum", "PD53343g_hum", "PD53344a_hum", "PD53344c_hum", "PD53344d_hum",
                                                  # "PD53344e_hum"))
#removing tumor samples
pdx_new_maf <- subsetMaf(maf = comb_maf, tsb = c("PD53330a_hum", "PD53332a_hum", "PD53333a_hum", "PD53343a_hum", "PD53343c_hum", 
                                                 "PD53347a_hum", "PD53349a_hum","PD53349c_hum", "PD53350a_hum", "PD53350c_hum", "PD53352a_hum", 
                                                 "PD53352c_hum", "PD53355a_hum", "PD53357a_hum", "PD53359a_hum", "PD53364a_hum"))
                                                 
topn<-20

#To make a figure 4c:
# p <- coBarplot(m1 = tumor_new_maf, m2 = pdx_new_maf, genes = c("KIT", "LRP1B", "NRAS", "ABL2", "ARHGEF12", "AXIN2", "CBL", "CTNNB1", "FAT4", "KIF5B"), m1Name = "Patients", m2Name = "PDXs", orderBy ='m1', normalize = TRUE)3
#ggsave("top10_genes_pdx_tumor_COSMIC_related_those_tumor_samples.pdf", plot = p$plot, width = 10, height = 8, units = "in", dpi = 300) 



#Saved in PDX_Tumor_samples folder
pdf(file = here("PDX_Tumor_samples", "Figure_4c.pdf"), width = 10, height = 8)
coBarplot(
  m1 = tumor_new_maf, 
  m2 = pdx_new_maf, 
  genes = c("NRAS", "KIT", "CUX1", "FAM47C", "FAT1", "HRAS", "KRAS", "NF1", "ABL2", "CBL", "LRP1B"),
  m1Name = "Patients", 
  m2Name = "PDXs", 
  orderBy = 'm1', 
  normalize = TRUE
)
dev.off()  # Fecha e salva o PDF

tumor_variant <- tumor_new_maf@variant.type.summary
pdx_variant <- pdx_new_maf@variant.type.summary

#salvando as tabelas
write.csv(tumor_variant , file = "tumor_variant .csv", row.names = TRUE)
write.csv(pdx_variant , file = "pdx_variant .csv", row.names = TRUE)


plotmafSummary(maf = tumor_new_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)



#Set the colours for the Tumour_site
pcol<-c("#295884", "#FCCDE5", "#FFCC99")
#pcol<-c("#295884", "#D20C60")
names(pcol)<- unique(tumor_new_maf@clinical.data$Tumor_Type)
#Set the colours for the  clinical_stage
clini_stage<- c("#CCEBC5", "#FFED6F", "#8DD3C7","#BC80BD","#BEBADA","#A6761D")  
names(clini_stage)<- unique(tumor_new_maf@clinical.data$Clinical_Stage )
#Set the colour for anatomic site
anatomic<- c("#1B9E77", "#B3DE69",  "#D20C60") 
names(anatomic)<- unique(tumor_new_maf@clinical.data$Anatomic_Site )

clin_col<-list(Tumor_Type=pcol,
               Clinical_Stage=clini_stage,
               Anatomic_Site=anatomic)

#to make the figure 2a
# Oncoplot showing ALLS COSMIC genes from the CANCER gene census only
topn<-20
# THIS IS how the overlaps should be taken 
observed_mut_canc_genes <-tumor_new_maf@gene.summary$Hugo_Symbol[tumor_new_maf@gene.summary$Hugo_Symbol %in% cosmic_genesv97$Gene.Symbol]  

oncoplot(tumor_new_maf,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not
         genes = observed_mut_canc_genes[1:topn], # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 0.9, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         annotationFontSize = 1.6,
         fontSize = 1,
         legendFontSize = 1.6, 
         barcode_mar = 4,
         gene_mar=8,
         #clinicalFeatures = c('Tumor_Type','Clinical_Stage','Anatomic_Site'),
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         #annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off() # The plot was saved in patient_samples_folder

#To make a figure 4c:
tumor_variant <- tumor_new_maf@variant.type.summary
pdx_variant <- pdx_new_maf@variant.type.summary

#to extract the number of SNP, DEL, INS for tumour and PDX-X1:
variant_number_table <- data.frame(
  Patients = tumor_variant$Tumor_Sample_Barcode,
  Number_SNP_DEL_INS_tumor = tumor_variant$SNP + tumor_variant$DEL + tumor_variant$INS,
  Number_SNP_DEL_INS_PDX = pdx_variant$SNP + pdx_variant$DEL + pdx_variant$INS
)

#ploting 
variant_number_table %>% ggplot(aes(x = Number_SNP_DEL_INS_tumor, y = Number_SNP_DEL_INS_PDX)) +
  geom_jitter(aes(color = Patients), width = 0.2, height = 0.2, alpha = 0.5, size = 3) +
  labs(x="Number of SNP/DEL/INS Tumor", y= "Number of SNP/DEL/INS PDX") +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "Black") +
  geom_text(x = Inf, y = -Inf, 
            label = paste("R²: ", round(cor_result$estimate, 2), "\n",
                          "p-value: ", format(cor_result$p.value, scientific = TRUE, 3), sep = ""),
            hjust = 1.1, vjust = -1, size = 4) +
  labs(x = "Number_SNP_DEL_INS_tumor", y = "Number_SNP_DEL_INS_PDX", title = "Correlation Plot") +
  theme_classic()


combined_plot <- grid.arrange(plot1, plot2, plot3, ncol = 3)
