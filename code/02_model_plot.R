## lung cancer reply
require("cancereffectsizeR")
require("ces.refset.hg19")
library(data.table)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(stringr)
library(ggplot2)

## load data
cesa <- load_cesa("luad_cesa_object_pan.rds")
cesa_smoking <- load_cesa("luad_cesa_object_smoker.rds")
cesa_nonsmoking <- load_cesa("luad_cesa_object_ns.rds")

sample_pan_forCsea <- readRDS("sample_pan_forCESA.rds")
## CES values of KRAS variants ####
target_effect_KRAS <- rbind(
  cesa$selection$recurrent[gene %in% c("KRAS")][, group := "All"][],
  cesa_smoking$selection$recurrent[gene %in% c("KRAS")][,
    group := "Smoker"
  ][],
  cesa_nonsmoking$selection$recurrents[gene %in% c("KRAS")][,
    group := "Never-smoker"
  ][]
)
target_effect_KRAS[,
  group := factor(
    group,
    levels = c("All", "Never-smoker", "Smoker"),
    ordered = T
  )
]
variant_order <- target_effect_KRAS["All", on = "group"][order(
  selection_intensity,
  decreasing = T
)][, variant_name]
target_effect_KRAS[,
  variant_name := factor(variant_name, levels = variant_order, ordered = T)
]
target_effect_KRAS[,
  prevalence := round(included_with_variant / (included_total + held_out), 3)
]
target_effect_KRAS[,
  ci_low_95_forPlot := ifelse(is.na(ci_low_95), 0, ci_low_95)
]

variant_order2 <- c(
  "G10V",
  "G12A",
  "G12D",
  "G12C",
  "G12V",
  "G12R",
  "G12S",
  "G13D",
  "G13C",
  "G13R",
  "L19F",
  "Q22K",
  "D33E",
  "A59G",
  "A59T",
  "Q61H",
  "Q61R",
  "A146T",
  "A146V"
)
target_effect_KRAS[,
  variant_name_s := stringr::str_split(variant_name, "_", simplify = T)[, 2]
]
target_effect_KRAS_sub <- target_effect_KRAS[variant_name_s %in% variant_order2]
target_effect_KRAS_sub[,
  variant_name_s := factor(variant_name_s, levels = variant_order2, ordered = T)
]

### Fig1 ####
## only show all group
p_ces_KRAS <- ggplot(
  target_effect_KRAS_sub[group %in% c("All") & !is.na(variant_name_s)],
  aes(x = variant_name_s, y = selection_intensity)
) +
  geom_errorbar(
    aes(ymin = ci_low_95_forPlot, ymax = ci_high_95),
    color = "azure4",
    na.rm = TRUE,
    width = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    shape = 21,
    color = "gray20",
    aes(fill = selection_intensity, size = prevalence)
  ) +
  scale_fill_gradient(
    low = "black",
    high = "red",
    breaks = c(25000, 50000, 75000, 100000),
    labels = c("25k", "50k", "75k", "100k")
  ) +
  labs(
    x = "KRAS variant",
    y = "Cancer effect size",
    fill = "Cancer effect",
    size = "Prevalence"
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5e5, 1e5),
    labels = c("0", expression(0.5 %*% 10^5), expression(1 %*% 10^5)) ## adjust the y-axis to display correctly, with tick labels like "0, 0.5 × 10^5, 1 × 10^5", using a superscript "5" instead of the caret notation (^5)
  ) +
  scale_size_continuous(
    breaks = c(0.01, 0.05, 0.10, 0.15),
    labels = c("1%", "5%", "10%", "15%"),
    limits = c(0, max(target_effect_KRAS_sub$prevalence)),
    range = c(2, 6)
  ) +
  theme_classic() +
  theme(
    legend.direction = "vertical",
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.2),
    legend.background = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, face = "bold", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
  ) +
  guides(
    fill = guide_colourbar(order = 1), # color legend first
    size = guide_legend(order = 2) # size legend second
  )

svg(file = "p_CES_KRAS.svg", height = 7, width = 13.5)
p_ces_KRAS
dev.off()


### check the frequency of KRAS G12C, 12D, G12A, G12C in smoker and ns-smoker####
head(target_effect_KRAS)
ggplot(target_effect_KRAS_sub) +
  geom_bar(
    aes(x = variant_name_s, y = prevalence),
    stat = "identity",
    color = "black",
    position = position_dodge(width = 0.9)
  ) +
  facet_wrap(~group, nrow = 3) +
  theme_classic()


## calculate TMB ####
maf_file <- as.data.table(maf_file)

## When calculating Tumor Mutational Burden (TMB), doublet base substitutions (DBS) are typically treated as two separate mutations rather than one. so we use maf_file instead of integrated_lungCancer_maf_sample here
maf_file_sample <- left_join(
  maf_file,
  cesa$samples,
  by = c("Tumor_Sample_Barcode" = "Unique_Patient_Identifier")
)

maf_file_sample_wes <- maf_file_sample[coverage %in% "exome"]
### way1 considering only non-synonymous mutations: #####
# "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site”

# Define non-synonymous mutations based on the priority list
non_synonymous <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Translation_Start_Site"
)

# Function to choose the most impactful classification from a composite classification
## how to deal with composite annotation:
# Separate Composite Annotations, Assign the mutation type based on the priority list, Exclude any classifications that are purely silent or non-coding.
# Keep the one with a higher impact, and count it as 1 mutation in TMB calculation
select_priority <- function(variant) {
  # Split the composite variant classification
  mut_types <- unlist(strsplit(variant, ","))
  # Retain only those that match the non-synonymous list
  impactful <- mut_types[mut_types %in% non_synonymous]
  # If multiple types, take the first one from the prioritized list
  if (length(impactful) > 0) {
    return(impactful[1])
  } else {
    return(NA)
  }
}

#### calculate TMB for WES data ####
# Apply the function to filter and clean up variant classifications
maf_file_sample_wes$Filtered_Classification <- sapply(
  maf_file_sample_wes$Variant_Classification,
  select_priority
)

# Remove rows with NA (non-impactful or silent mutations)
maf_file_sample_wes_filtered <- maf_file_sample_wes[
  !is.na(Filtered_Classification)
]
table(maf_file_sample_wes_filtered$Filtered_Classification)

exome_size_mb <- 38 # Adjust this number as needed

# Count non-synonymous mutations per sample
tmb_table_wes <- maf_file_sample_wes_filtered[, .N, by = Tumor_Sample_Barcode]
tmb_table_wes[, TMB := N / exome_size_mb]
tmb_table_wes[, coverage := "exome"]
head(tmb_table_wes)
head(setorder(tmb_table_wes, -TMB)) ## the highest TMB value is 70
head(setorder(tmb_table_wes, TMB))
dim(tmb_table_wes) ##TMB infor of 1325 WES sample

### way2 considering all somatic mutations #####
#### calculate TMB for WES data ####
tmb_table_wes_v2 <- maf_file_sample_wes[, .N, by = Tumor_Sample_Barcode]
exome_size_mb <- 38 # Adjust this number as needed
tmb_table_wes_v2[, TMB := N / exome_size_mb]
tmb_table_wes_v2[, coverage := "exome"]
head(setorder(tmb_table_wes_v2, -TMB))
summary(tmb_table_wes_v2$TMB) ## the max TMB is 92 here

## TMB- the distribution of TMB of samples with four KRAS variants (KRAS G12D mutation is correlated with lower TMB)####
KRAS_mut_tar <- c("KRAS_G12D", "KRAS_G12C", "KRAS_G12V", "KRAS_G12A")


integrated_lungCancer_maf_sample <- left_join(
  cesa$maf,
  cesa$samples,
  by = "Unique_Patient_Identifier"
) ## need to add sample info

## two or more of the KRAS G12 mutations (G12D, G12C, G12V, G12A) cannot coexist in the same tumor sample because they occur at the same genomic position (codon 12 of the KRAS gene).
## KRAS mutations at codon 12 are oncogenic driver mutations.
## Tumors typically acquire a single driver mutation per oncogene because additional mutations in the same codon do not provide additional selective advantage.\
KRAS_mut_tar_sample <- integrated_lungCancer_maf_sample[
  top_consequence %in% KRAS_mut_tar
]
dim(KRAS_mut_tar_sample) ## 2659, 18
head(KRAS_mut_tar_sample)
table(KRAS_mut_tar_sample$coverage, exclude = NULL) ## only 256 WES samples contain these four mutations
table(KRAS_mut_tar_sample$top_consequence, exclude = NULL) ## there are 1279 sample with G12C mut


### V1: TMB was calculated with only non-synonymous mutations####
KRAS_mut_tar_sample_TMB_wes <- left_join(
  KRAS_mut_tar_sample[coverage %in% c("exome")],
  tmb_table_wes,
  by = c("Unique_Patient_Identifier" = "Tumor_Sample_Barcode")
)
dim(KRAS_mut_tar_sample_TMB_wes) ## 256  21
head(KRAS_mut_tar_sample_TMB_wes)
KRAS_mut_tar_sample_TMB_wes[,
  top_consequence := factor(top_consequence, levels = KRAS_mut_tar, ordered = T)
]
ggplot(data = KRAS_mut_tar_sample_TMB_wes, aes(x = top_consequence, y = TMB)) +
  geom_boxplot(aes(fill = top_consequence)) +
  labs(x = "KRAS mutation", y = "Tumor mutation burden") +
  scale_x_discrete(labels = c("G12D", "G12C", "G12V", "G12A")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold")
  ) +
  guides(fill = "none") +
  ggpubr::stat_compare_means(
    aes(label = ..p.signif..),
    comparisons = list(
      c("KRAS_G12D", "KRAS_G12C"),
      c("KRAS_G12D", "KRAS_G12V"),
      c("KRAS_G12D", "KRAS_G12A")
    ),
    method = "wilcox.test"
  )


### V2: TMB was calculated with all somatic mutations ## Fig2C  ####
KRAS_mut_tar_sample_TMB_wes_v2 <- left_join(
  KRAS_mut_tar_sample[coverage %in% c("exome")],
  tmb_table_wes_v2,
  by = c("Unique_Patient_Identifier" = "Tumor_Sample_Barcode")
)
dim(KRAS_mut_tar_sample_TMB_wes_v2) ## 256  21
head(KRAS_mut_tar_sample_TMB_wes_v2)
KRAS_mut_tar_sample_TMB_wes_v2[,
  top_consequence := factor(top_consequence, levels = KRAS_mut_tar, ordered = T)
]
p_G12_TMB <- ggplot(
  data = KRAS_mut_tar_sample_TMB_wes_v2,
  aes(x = top_consequence, y = TMB)
) +
  geom_boxplot(aes(fill = top_consequence)) + ## , outlier.shape = NA
  labs(x = "KRAS mutation", y = "Tumor mutation burden") +
  scale_x_discrete(labels = c("G12D", "G12C", "G12V", "G12A")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold")
  ) +
  guides(fill = "none") +
  ggpubr::stat_compare_means(
    aes(label = ..p.signif..),
    comparisons = list(
      c("KRAS_G12D", "KRAS_G12C"),
      c("KRAS_G12D", "KRAS_G12V"),
      c("KRAS_G12D", "KRAS_G12A")
    ),
    method = "wilcox.test"
  )

ggsave(file = "p_G12_TMB.svg", plot = p_G12_TMB)
## Fig: the distribution of SBS4 attribution in tumors across different KRAS common mutation subtypes (G12D mut exhibited significantly lower SBS4) ####

### attribution analysis ####
biological_wei <- cesa$mutational_signatures$biological_weights
ncol(biological_wei)
SBS_weight <- biological_wei[, (ncol(biological_wei) - 77):ncol(biological_wei)]
SBS_weight_sum <- rowSums(SBS_weight) ##
unique(SBS_weight_sum)

head(biological_wei)
KRAS_mut_tar_sample_biologWei <- inner_join(
  KRAS_mut_tar_sample,
  biological_wei,
  by = c("Unique_Patient_Identifier")
)
dim(KRAS_mut_tar_sample_biologWei) ## 2604 99 ## there are 54 samples that do not contain biological weights information
length(unique(KRAS_mut_tar_sample_biologWei$Unique_Patient_Identifier)) ## 2604
head(setorder(KRAS_mut_tar_sample_biologWei, -SBS4))
KRAS_mut_tar_sample_biologWei[,
  top_consequence := factor(top_consequence, levels = KRAS_mut_tar, ordered = T)
]
summary(KRAS_mut_tar_sample_biologWei$SBS4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.3419  0.3419  0.3521  0.3419  1.0000
summary(KRAS_mut_tar_sample_biologWei[coverage %in% c("targeted")]$SBS4) ## targeted sample should be removed as they do no really participated in signature attribution analysis
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3419  0.3419  0.3419  0.3431  0.3419  0.3864
summary(KRAS_mut_tar_sample_biologWei[!coverage %in% c("targeted")]$SBS4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.2799  0.4690  0.4324  0.5882  1.0000
length(KRAS_mut_tar_sample_biologWei$Unique_Patient_Identifier)
length(unique(KRAS_mut_tar_sample_biologWei$Unique_Patient_Identifier))
### Fig 2A #####
p_G12C_SBS4 <- ggplot(
  data = KRAS_mut_tar_sample_biologWei[!coverage %in% c("targeted")],
  aes(x = top_consequence, y = SBS4)
) +
  geom_boxplot(aes(fill = top_consequence), outlier.shape = NA) + ## ,outlier.shape = NA
  scale_y_continuous(
    limits = c(0, 1.4),
    expand = c(0, 0),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  labs(
    x = "KRAS mutation",
    y = "       Proportion of mutations
       attributed to SBS4 signature"
  ) +
  scale_x_discrete(labels = c("G12D", "G12C", "G12V", "G12A")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(
      size = 18,
      face = "bold",
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.title.x = element_text(size = 18, face = "bold")
  ) +
  guides(fill = "none") +
  ggpubr::stat_compare_means(
    aes(label = ..p.signif..),
    comparisons = list(
      c("KRAS_G12D", "KRAS_G12C"),
      c("KRAS_G12D", "KRAS_G12V"),
      c("KRAS_G12D", "KRAS_G12A")
    ),
    method = "wilcox.test",
    y.position = c(0.95, 0.98, 1.01)
  )

ggsave(file = "p_G12_SBS4.svg", plot = p_G12C_SBS4)


## linear regression of TMB   and SBS4 attribution (strong positive association, suggesting  that smoking is a major contribution to TMB in lung cancer#####
### V1: TMB was calculated with only non-synonymous mutations ####
tmb_table_wes_biologWei <- inner_join(
  tmb_table_wes,
  biological_wei,
  by = c("Tumor_Sample_Barcode" = "Unique_Patient_Identifier")
) ## "Unique_Patient_Identifier" = "Tumor_Sample_Barcode"
dim(tmb_table_wes_biologWei) ## 1325-1256; 69 WES sample without biological weight information
tmb_table_wes_biologWei[, log10TMB := log10(TMB), by = .I]
head(tmb_table_wes_biologWei)
lm_result <- lm(log10TMB ~ SBS4, data = tmb_table_wes_biologWei)
summary(lm_result)
ggplot(data = tmb_table_wes_biologWei, aes(x = SBS4, y = log10TMB)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(
    x = "Proportion of mutations 
       attributed to SBS4 signature",
    y = expression(Log[10] ~ "of tumor mutation burden")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold")
  )
### V2: TMB was calculated with all somatic mutations ####
tmb_table_wes_biologWei_v2 <- inner_join(
  tmb_table_wes_v2,
  biological_wei,
  by = c("Tumor_Sample_Barcode" = "Unique_Patient_Identifier")
) ## "Unique_Patient_Identifier" = "Tumor_Sample_Barcode"
dim(tmb_table_wes_biologWei_v2) ## 1325-1256; 69 WES sample without biological weight information
tmb_table_wes_biologWei_v2[, log10TMB := log10(TMB), by = .I]
head(tmb_table_wes_biologWei_v2)
lm_result_v2 <- lm(log10TMB ~ SBS4, data = tmb_table_wes_biologWei_v2)
summary(lm_result_v2) ## slop is 1.040, p-value is <2.2e-16

ggplot(data = tmb_table_wes_biologWei_v2, aes(x = SBS4, y = log10TMB)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(
    x = "Proportion of mutations 
       attributed to SBS4 signature",
    y = expression(Log[10] ~ "of tumor mutation burden")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold")
  )

#### Fig2B use 50 SNV as the cutoff, as samples with < 50 SNV are not qualified for signature analysis ####
lm_result_v2_filter <- lm(
  log10TMB ~ SBS4,
  data = tmb_table_wes_biologWei_v2[TMB > 50 / exome_size_mb]
)
summary(lm_result_v2_filter) ## slope is 0.9462, the P-value is <2.2e-16


p_SBS4_TMB <- ggplot(
  data = tmb_table_wes_biologWei_v2[TMB > 50 / exome_size_mb],
  aes(x = SBS4, y = log10TMB)
) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(
    x = "       Proportion of mutations 
       attributed to SBS4 signature",
    y = bquote(bold(Log[10] ~ "of tumor mutation burden"))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(
      size = 18,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    )
  ) +
  scale_y_continuous(
    limits = c(0, 2.0),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.05, 0)))


ggsave(file = "p_SBS4_TMB.svg", plot = p_SBS4_TMB)

## tobacco exposure result in G12C: Significant enrichment of KRAS G12C mutations in smokers compared to non-smokers. ####
### Smoking status here was determined by clinical record for WES/WGS/TGS data####
clin_df = fread(paste0(location_output, 'merged_luad_clinical.txt'))
clin_df_sampleInfo <- inner_join(
  clin_df,
  cesa$samples,
  by = c("Sample ID" = "Unique_Patient_Identifier")
)
dim(clin_df_sampleInfo) ## 9230, 15
table(clin_df_sampleInfo$Smoker, clin_df_sampleInfo$coverage, exclude = NULL)

integrated_lungCancer_maf_sample[,
  G12C_sta := ifelse(top_consequence %in% "KRAS_G12C", "G12C", "No_G12C")
]
sample_G12C_sta <- unique(integrated_lungCancer_maf_sample[, .(
  Unique_Patient_Identifier,
  G12C_sta,
  coverage
)])
dim(sample_G12C_sta) ## 9457, 3
## add G12_status
clin_df_sampleInfo <- inner_join(
  clin_df_sampleInfo,
  sample_G12C_sta,
  by = c("Sample ID" = "Unique_Patient_Identifier")
)
head(clin_df_sampleInfo)

setnames(
  clin_df_sampleInfo,
  old = "Smoker",
  new = "smoking_sta_byCliniDataOnly"
)

contingency_table_byCliniData <- table(
  clin_df_sampleInfo$smoking_sta_byCliniDataOnly,
  clin_df_sampleInfo$G12C_sta
)
# G12C No_G12C
# FALSE   20     760
# TRUE   161     858
print(contingency_table_byCliniData)
contingency_table_byCliniData
dimnames(contingency_table_byCliniData) = list(
  c("Never-Smoker", "Smoker"),
  c("G12C", "No_G12C")
) ## rename the contigenct_table
chisq_result_byCliniData <- chisq.test(contingency_table_byCliniData)
print(chisq_result_byCliniData)
# Pearson's Chi-squared test with Yates' continuity correction
#
# data:  contingency_table_byCliniData
# X-squared = 84.077, df = 1, p-value < 2.2e-16

## SBS4 -> TMB -> IO response ####
### Litch，logistic model, RR-TMB####
# Define the functions
# 1. Define the regression coefficients from SBS4 to log10(TMB)
b0 <- coef(lm_result_v2_filter)[1] ## intercept_from_your_regression  # e.g., -1
b1 <- coef(lm_result_v2_filter)[2] ## slope_from_your_regression      # e.g., 2
predict_log10TMB <- function(SBS4) {
  b0 + b1 * SBS4
}
# 2. Define parameters for the logistic model
OR <- 1.7
p0_TMB <- 0.03 ## 0.03 or 0.05 is more reasonable
a1 <- log(OR) ## 0.5306
a0_TMB <- log(p0_TMB / (1 - p0_TMB)) ## -2.19 for p0 <- 0.1
# Define a function to predict IO response probability from log10(TMB), use a simple logistic model for now

predict_IO_response_TMB <- function(log10TMB) {
  1 / (1 + exp(-(a0_TMB + a1 * (10^log10TMB)))) ## probability form of logistic model is : 1/(1+exp(-(a0+a1*x))) ## x for Litchfield et al is Z_TMB
}
# use real SBS4 data, should not include panel data
nrow(KRAS_mut_tar_sample_biologWei[(!coverage %in% c("targeted"))]) ## 262 sample
nrow(KRAS_mut_tar_sample_biologWei[(coverage %in% c("targeted"))]) ## 2342 panel sample
table(
  KRAS_mut_tar_sample_biologWei[
    (!coverage %in% c("targeted")),
    .(top_consequence, coverage)
  ],
  exclude = NULL
)

SBS4_G12C <- KRAS_mut_tar_sample_biologWei[
  (!coverage %in% c("targeted")) & (top_consequence %in% "KRAS_G12C"),
  SBS4
]
SBS4_G12D <- KRAS_mut_tar_sample_biologWei[
  (!coverage %in% c("targeted")) & (top_consequence %in% "KRAS_G12D"),
  SBS4
]

# Create data frames
df_G12C <- data.frame(
  SBS4 = SBS4_G12C,
  KRAS_variant = "G12C"
)
df_G12D <- data.frame(
  SBS4 = SBS4_G12D,
  KRAS_variant = "G12D"
)

# Combine them
df_all <- rbind(df_G12C, df_G12D)

# Predict log10TMB and IO response
df_all$log10TMB <- predict_log10TMB(df_all$SBS4)
df_all$TMB <- 10^(df_all$log10TMB)

## check the real TMB of G12C/G12D samples
head(KRAS_mut_tar_sample_biologWei)
head(tmb_table_wes_v2)
KRAS_mut_tar_sample_biologWei_TMB_wes_v2 <- inner_join(
  KRAS_mut_tar_sample_biologWei,
  tmb_table_wes_v2,
  by = c("Unique_Patient_Identifier" = "Tumor_Sample_Barcode")
)

df_all$IO_response_prob_TMB_plateau <- 0.6 *
  predict_IO_response_TMB(df_all$log10TMB) ## add 0.6 here

### Samstein et al, four-parameter logistic (4PL) model, HR-TMB####
TMB_hazardRatio <- fread(file = "TMB_HazardRatio_Data.csv")
colnames(TMB_hazardRatio) <- c("TMB", "HR")
#### fit the curve with the data, here use a four-parameter logistic (4PL) curve, common in dose-response:HR(TMB) = d+ (a-d)/(1+(TMB/c)^b)
# install.packages("drc")
library(drc)
model <- drm(HR ~ TMB, fct = LL.4(), data = TMB_hazardRatio)
summary(model)
# Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#
# Parameter estimates:
#
#                Estimate Std. Error t-value   p-value
# b:(Intercept) 2.3466487  0.1992005  11.780 < 2.2e-16 ***
#   c:(Intercept) 0.3561694  0.0059638  59.722 < 2.2e-16 ***
#   d:(Intercept) 1.0111100  0.0498118  20.299 < 2.2e-16 ***
#   e:(Intercept) 7.8421479  0.5458201  14.368 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error:
#
#   0.02466389 (99 degrees of freedom)
predict_HR_TMB <- function(TMB) {
  c + (d - c) / (1 + (TMB / e)^b) ## probability form of logistic model is : 1/(1+exp(-(a0+a1*x))) ## x for Litchfield et al is Z_TMB
}
params <- coef(model)

b <- params["b:(Intercept)"]
c <- params["c:(Intercept)"]
d <- params["d:(Intercept)"]
e <- params["e:(Intercept)"]
df_all$HR_TMB <- predict_HR_TMB(df_all$TMB)


### add two short vertical lines at the median SBS4 for KRAS G12C and KRAS G12D to highlight them
df_all <- as.data.table(df_all)
median_G12C <- median(df_all$SBS4[df_all$KRAS_variant == "G12C"], na.rm = TRUE)
median_G12D <- median(df_all$SBS4[df_all$KRAS_variant == "G12D"], na.rm = TRUE)
df_all$KRAS_variant <- factor(df_all$KRAS_variant, levels = c("G12C", "G12D"))
setorder(df_all, KRAS_variant) ## In ggplot2, when you have geom_point(aes(fill = KRAS_variant)), ggplot plots the points in the order that they appear in your data. as G12C is infront of G12D, so G12D points which were plotted later will be on the top of G12C

#### Fig 2 D ####
RR_median_G12C <- median(df_all[
  KRAS_variant == "G12C",
  IO_response_prob_TMB_plateau
])
RR_median_G12D <- median(df_all[
  KRAS_variant == "G12D",
  IO_response_prob_TMB_plateau
])
# Define the RR(SBS4) function
RR_from_SBS4 <- function(SBS4) {
  TMB <- 10^(0.513 + 0.946 * SBS4)
  RR <- 0.6 * (1 / (1 + exp(-(-3.476 + 0.531 * TMB))))
  return(RR)
}

# Create a data frame with a sequence of SBS4 values
sbs4_seq <- seq(0, 1, length.out = 200)
df_RR <- data.frame(
  SBS4 = sbs4_seq,
  RR = RR_from_SBS4(sbs4_seq)
)
p_RR_SBS4_line_point_blackline <- ggplot() +
  geom_line(data = df_RR, aes(x = SBS4, y = RR), color = "black", size = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(0.25, 0.5, 0.75, 1)
  ) +
  geom_segment(
    data = df_all,
    aes(x = median_G12C, xend = median_G12C, y = 0, yend = RR_median_G12C),
    linetype = "dashed",
    color = "#F8766D"
  ) +
  geom_segment(
    data = df_all,
    aes(x = 0, xend = median_G12C, y = RR_median_G12C, yend = RR_median_G12C),
    linetype = "dashed",
    color = "#F8766D"
  ) +
  geom_segment(
    data = df_all,
    aes(x = median_G12D, xend = median_G12D, y = 0, yend = RR_median_G12D),
    linetype = "dashed",
    color = "#00BFC4"
  ) +
  geom_segment(
    data = df_all,
    aes(x = 0, xend = median_G12D, y = RR_median_G12D, yend = RR_median_G12D),
    linetype = "dashed",
    color = "#00BFC4"
  ) +
  geom_segment(
    data = df_all,
    aes(x = SBS4, y = 0.03, yend = 0.03 + 0.04, color = KRAS_variant),
    size = 0.5,
    alpha = 0.7
  ) +
  theme_classic() +
  labs(
    x = "       Proportion of mutations 
       attributed to SBS4 signature",
    y = "Objective response rate",
    color = "KRAS Variant",
  ) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = c(0.75, 0.2),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold")
  )


ggsave(file = "p_RR_SBS4.svg", plot = p_RR_SBS4_line_point_blackline)

#### Fig 2 E ####
HR_median_G12C <- median(df_all[KRAS_variant == "G12C", HR_TMB])
HR_median_G12D <- median(df_all[KRAS_variant == "G12D", HR_TMB])
# Define the HR(SBS4) function
HR_from_SBS4 <- function(SBS4) {
  TMB <- 10^(0.513 + 0.946 * SBS4)
  HR <- 0.356 + (1.011 - 0.356) / (1 + (TMB / 7.84)^2.35)
  return(HR)
}

# Create a data frame with a sequence of SBS4 values
sbs4_seq <- seq(0, 1, length.out = 200)
df_HR <- data.frame(
  SBS4 = sbs4_seq,
  HR = HR_from_SBS4(sbs4_seq)
)
p_HR_SBS4_line_point_blackline <- ggplot() +
  geom_line(data = df_HR, aes(x = SBS4, y = HR), color = "black", size = 2) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(0.25, 0.5, 0.75, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  geom_segment(
    data = df_all,
    aes(x = median_G12C, xend = median_G12C, y = 0, yend = HR_median_G12C),
    linetype = "dashed",
    color = "#F8766D"
  ) +
  geom_segment(
    data = df_all,
    aes(x = 0, xend = median_G12C, y = HR_median_G12C, yend = HR_median_G12C),
    linetype = "dashed",
    color = "#F8766D"
  ) +
  geom_segment(
    data = df_all,
    aes(x = median_G12D, xend = median_G12D, y = 0, yend = HR_median_G12D),
    linetype = "dashed",
    color = "#00BFC4"
  ) +
  geom_segment(
    data = df_all,
    aes(x = 0, xend = median_G12D, y = HR_median_G12D, yend = HR_median_G12D),
    linetype = "dashed",
    color = "#00BFC4"
  ) +
  geom_segment(
    data = df_all,
    aes(x = SBS4, y = 0.03, yend = 0.03 + 0.04, color = KRAS_variant),
    size = 0.5,
    alpha = 0.7
  ) +
  theme_classic() +
  labs(
    x = "       Proportion of mutations 
       attributed to SBS4 signature",
    y = "Hazard ratio",
    color = "KRAS Variant" #,
  ) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = c(0.75, 0.2),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold")
  )

ggsave(file = "p_HR_SBS4.svg", plot = p_HR_SBS4_line_point_blackline)
