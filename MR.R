library(VariantAnnotation)
library(gwasglue)
library(data.table)
library(dplyr)
library(tidyr)
library(TwoSampleMR)


setwd("D:/MR")

#read exposure
his <- fread("./exposure/34594039-GCST90019001-EFO_0009943.h.tsv.gz", data.table = F)
#filter exposure
his_filter <- his[his$p_value<5e-7,]


his_filter <- format_data(dat = his_filter,
                          type = "exposure",
                          snp_col = "hm_rsid",
                          beta_col = "beta",
                          se_col = "standard_error",
                          eaf_col = "effect_allele_frequency",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          pval_col = "p_value")
colnames(his_filter)



#local clumping
his_clumped <- ld_clump(
  dat = dplyr::tibble(rsid=his_filter$SNP,
                      pval=his_filter$pval.exposure,
                      id=his_filter$exposure),
  clump_kb=500, clump_r2=0.05,clump_p = 1,
  pop = "EAS",
  bfile = "./1kg.v3/EAS", 
  plink_bin = "./plink")


#filter SNM
library(dplyr)
his_clumped <- his_filter %>% filter(SNP %in% his_clumped$rsid)

duplicated(his_clumped$rsid,his_clumped$rsid)
is_duplicate <- his_clumped$SNP %in% his_clumped$rsid
table(is_duplicate)


#read outcome
outcome_data <- fread("./outcome/finngen_R7_ST19_SEQUEL_INJURI_LOWER_LIMB.gz", data.table = F)


#format finn
outcome_data <- format_data(dat = outcome_data,
                            type = "outcome", # 类型是outcome
                            snp_col = "rsids",
                            phenotype_col = "phenotype",
                            beta_col = "beta",
                            se_col = "sebeta",
                            eaf_col="af_alt",
                            effect_allele_col = "alt",
                            other_allele_col = "ref",
                            pval_col = "pval") 


#merge
common_snp <- merge(his_clumped,outcome_data,by="SNP")$SNP
length(common_snp)

outcome_data1 <- outcome_data[outcome_data$SNP %in% common_snp,]

#harmonise
dat_har <- harmonise_data(his_clumped, outcome_data1)



#analysis
res <- mr(dat_har)
res

write.csv(res, "mrresult.csv", row.names=F)

#OR calculation
or <- generate_odds_ratios(res)


#plots
library(ggplot2)
library(ggsci)
p1 <- mr_scatter_plot(res, dat = dat_har)
p1
ggsave(p1[[1]],filename = "k.pdf",dpi = 300)


#sensitive check

#heterogeneity
mr_heterogeneity(dat_har)


#pleiotropy test
mr_pleiotropy_test(dat_har)


#single SNP
res_single <- mr_singlesnp(dat_har)
head(res_single)

p2 <- mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]] , filename = "forest.pdf", dpi = 300)


p3 <- mr_funnel_plot(res_single)
p3[[1]]
ggsave(p3[[1]],filename = "funnel.pdf",dpi = 300)



#leaveoneout
res_loo <- mr_leaveoneout(dat_har)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
ggsave(p4[[1]],filename = "leaveoneout.pdf",dpi = 300)