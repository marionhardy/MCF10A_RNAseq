
# Merge .csv count files

## Assumes you have run STAR on Linux
# And you have an output like : "ReadsPerGene.out.tab"

library(DESeq2)
library(tidyverse)
library(stringr)
library(biomaRt)
library(Hmisc)

# Im missing some data right now so some conditions dont have biological replicates

names = c('GM_none_1',
          'IM_none_1','IM_none_2',
          'IM_none_Glut_1','IM_none_Glut_2',
          'IM_none_EGF_1','IM_none_EGF_2',
          'IM_none_noIns_1','IM_none_noIns_2',
          'IM_none_noGluc_1','IM_none_noGluc_2',
          'IM_none_noHC_1','IM_none_noHC_2',
          'IM_none_noCT_1','IM_none_noCT_2',
          'IM_none_AMPKact_1','IM_none_AMPKact_2',
          'IM_none_AKTinhib_1','IM_none_AKTinhib_2',
          'IM_none_ERKinhib_1','IM_none_ERKinhib_2',
          'IM_none_mTORC1inhib_1','IM_none_mTORC1inhib_2',
          'IM_none_Oligo_1','IM_none_Oligo_2',
          'IM_none_MPCinhib_1','IM_none_MPCinhib_2',
          'IM_none_LDHinhib_1','IM_none_LDHinhib_2',
          'IM_none_IL6_1','IM_none_IL6_2')

files = list.files("./data/", pattern = "ReadsPerGene.out.tab")
files # check that the order of names and files is correct or you will 
# assign wrong names to samples
# itll go 10-11-12-etc-16-1-2-3-...

names = c('IM_none_AKTinhib_1','IM_none_AKTinhib_2',
          'IM_none_ERKinhib_1','IM_none_ERKinhib_2',
          'IM_none_mTORC1inhib_1','IM_none_mTORC1inhib_2',
          'IM_none_Oligo_1','IM_none_Oligo_2',
          'IM_none_MPCinhib_1','IM_none_MPCinhib_2',
          'IM_none_LDHinhib_1','IM_none_LDHinhib_2',
          'IM_none_IL6_1','IM_none_IL6_2',
          'GM_none_1',
          'IM_none_1','IM_none_2',
          'IM_none_Glut_1','IM_none_Glut_2',
          'IM_none_EGF_1','IM_none_EGF_2',
          'IM_none_noIns_1','IM_none_noIns_2',
          'IM_none_noGluc_1','IM_none_noGluc_2',
          'IM_none_noHC_1','IM_none_noHC_2',
          'IM_none_noCT_1','IM_none_noCT_2',
          'IM_none_AMPKact_1','IM_none_AMPKact_2')


for(i in 1:length(files)){
  x = read_tsv(paste0("./data/",files[i]), col_names = F)
  x = x[-c(1:5),-c(3,4)]
  colnames(x) = c("ensembl",paste0(names[i]))
  assign(names[i],x)
}

noquote(names) # Use that to quickly copy paste the names without the quotes

counts <- list(GM_none_1 ,       
               IM_none_1    ,      IM_none_2  ,
               IM_none_Glut_1  ,   IM_none_Glut_2,     IM_none_EGF_1   ,  
               IM_none_EGF_2 ,     IM_none_noIns_1  ,  IM_none_noIns_2  ,
               IM_none_noGluc_1  , IM_none_noGluc_2  , IM_none_noHC_1,    
               IM_none_noHC_2  ,   IM_none_noCT_1 ,    IM_none_noCT_2 ,
               IM_none_AMPKact_1,  IM_none_AMPKact_2 , IM_none_AKTinhib_1,
               IM_none_AKTinhib_2, IM_none_ERKinhib_1, IM_none_ERKinhib_2 ,
               IM_none_mTORC1inhib_1, IM_none_mTORC1inhib_2,
               IM_none_Oligo_1 ,   IM_none_Oligo_2 ,   IM_none_MPCinhib_1,
               IM_none_MPCinhib_2, IM_none_LDHinhib_1, IM_none_LDHinhib_2,
               IM_none_IL6_1, IM_none_IL6_2) %>% 
  purrr::reduce(full_join, by = "ensembl")


strrep =
  sub(pattern = "\\.(.*)","",counts$ensembl)

counts$ensembl = strrep
# How to deal with gene version duplicates?
table(duplicated(counts$ensembl))
table(is.na(counts$ensembl))

counts = counts %>% distinct(ensembl, .keep_all = TRUE)

rownames(counts) = counts$ensembl

write.csv(counts,"./data/MCF10A_counts.csv", row.names = T)

# Create the coldata for the summarized experiment

coldata <- data.frame(
  celltype =c(rep('MCF10A',each = 31)),
  condition = as.factor(c('Growth_Medium',rep(c('IM','IM_Glutamine_pos','IM_EGF_pos','IM_Ins_neg',
                        'IM_Glucose_neg','IM_HC_neg','IM_CT_neg','AMPK_activator',
                        'AKT_inhibitor','ERK_inhibitor','mTORC1_inhibitor',
                        'Oligomycin','MPC_inhibitor','LDH_inhibitor','IL6'), each = 2, times = 1))),
  condition_detail = as.factor(c('Growth_Medium',rep(c('IM','IM_Glutamine_pos','IM_EGF_pos','IM_Ins_neg',
                                    'IM_Glucose_neg','IM_HC_neg','IM_Cholera_toxin_neg','MK8722',
                                    'Ipasertip','PD0325901','Rapamycin',
                                    'Oligomycin','UK5099','Galloflavin','IL6'), each = 2, times = 1))),
  replicate=as.factor(c(1,rep(c(1:2),15))))

rownames(coldata) <- colnames(counts)[-1]

write.csv(coldata,"./data/coldata.csv")

