########################################################################################
######## MAIN SCRIPT FOR CANDIDATE GENE ANALYSES #######################################
########################################################################################

## This script contains the data cleaning/manipulation and regression analyses sections
### for the candidate gene analyses. The population structure analyses require the outputs
### of the 'phylogenetic clustering' scripts.

## This script contains the following sections. Each regression section contains
## univariate and multivariate regressions, with and without adjustment for population
## structure.
  ## Data loading, cleaning and manipulation
  ## Logistic regressions - presence of variants in candidate genes 
  ## Logistic regressions - number of variants per candidate gene
  ## Logistic regressions - by variant
  ## Data manipulation - continuous phenotype
  ## Linear regressions - presence of variants in candidate genes 
  ## Linear regressions - number of variants per candidate gene
  ## Linear regressions - by variant


library(tidyverse)
library(tidyr)
library(stringr)
library(dplyr)


#Read in file (known genes variants)
known_variants <- read_csv("efm_dap_dataset_v1.known_dap_genes.genotype_table.csv")
head(known_variants)

#Data wrangling --> by gene
known_variants_names <- colnames(known_variants)
known_variants3 <- rbind(known_variants, known_variants_names)
colnames(known_variants3)<-c(1)

ncols <- max(str_count(known_variants3[671,], "\t"))+1
ncols #1701
colmn <- paste(1:ncols)
colmn <- as.character(colmn)
ncol(known_variants3)
known_variants3 <- separate(known_variants3,
                            col=1,
                            into=colmn,
                            sep="\t")

known_variants4 <- data.frame(c(known_variants_names))
ncols2 <- max(str_count(known_variants4[1,], "\t"))+1
ncols2 #1701
colmn2 <- paste(1:ncols2)
colmn2 <- as.character(colmn2)

known_variants4 <- separate(known_variants4,
                            col=1,
                            into=colmn2,
                            sep = "\t")
colnames(known_variants3) <- c(known_variants4)
View(known_variants3)

#Importing and data wrangling of annotation file

known_annotation <- read_csv("efm_dap_dataset_v1.known_dap_genes.ann_table.csv")
View(known_annotation) #need to change column header to row1 

known_annotation_row1 <- colnames(known_annotation)
known_annotation <- rbind(known_annotation_row1, known_annotation)
colnames(known_annotation) <- c(1)

known_annotation <- separate(known_annotation,
                             col=1,
                             into=c("variant", "gene", "amino_acid_change", "variant_type", "severity"),
                             sep = "\t")
View(known_annotation)

#Merge dataframes

known_annotation <- t(known_annotation) #transpose annotation part
colnames(known_annotation) <- known_annotation[1,] #make variable names 'variant'
known_annotation <- known_annotation[-1,] #remove duplicate row
known_variants3 <- column_to_rownames(known_variants3, var="sample") #first make them have the same number of columns(variant)

known_full_table <- rbind(known_annotation, known_variants3) #concatenate by rows - this should match column names?
View(known_full_table)

View(known_full_table$CP003351.C.903984.A)

#Remove synonymous mutations (no amino acid change)

transp_known_full_table <- t(known_full_table) #transpose again so severity is a column
View(transp_known_full_table)

transp_known_full_table <- as.data.frame(transp_known_full_table)
transp_known_full_table <- rownames_to_column(transp_known_full_table)
transp_known_full_table <- as_tibble(transp_known_full_table) #convert back to tibble (lose rownames but we don't need these for next step) to do filter
View(transp_known_full_table)
known_table_nonsyn <- filter(transp_known_full_table, severity!="LOW") #NB no rownames. 672 samples x548variants

known_table_nonsyn <-known_table_nonsyn[-167,] #remove 2872 variant
dim(known_table_nonsyn) # [1] 547 676
View(known_table_nonsyn)

#Merge columns by gene

known_table_nonsyn1 <- known_table_nonsyn %>%
  select(!c(rowname, amino_acid_change,variant_type, severity))#first get rid of non-gene extra columns 
View(known_table_nonsyn1)

known_table_nonsyn1 <- select(known_table_nonsyn1, -sample)
View(known_table_nonsyn1) #547x671

known_table_nonsyn1 <- pivot_longer(known_table_nonsyn1, cols=2:671, names_to = "sample_id", values_to="count")

known_table_nonsyn1$count <- as.numeric(known_table_nonsyn1$count)
known_table_nonsyn2 <- known_table_nonsyn1 %>%
  aggregate(by=list(known_table_nonsyn1$sample_id, known_table_nonsyn1$gene), FUN=mean)
head(known_table_nonsyn2)

View(known_table_nonsyn2) #think this works for the moment, though not tidy

mutated_samples_per_gene <- aggregate(known_table_nonsyn1$count, by=list(known_table_nonsyn1$gene), FUN=sum)
mutated_samples_per_gene$Count <- mutated_samples_per_gene$x/2
mutated_samples_per_gene <- select(mutated_samples_per_gene, -x)
View(mutated_samples_per_gene)
write.table(mutated_samples_per_gene, file="mutated_samples_per_gene.txt", sep="\t", row.names = F)

#phenotypic resistance metadata - binary

phenotype <- read_csv("efm_dap_dataset_v1.daptomycin_metadata.csv")
phenotype <- separate(phenotype,
                      col=1,
                      into=c("sample_id", "daptomycin_mic", "daptomycin_MIC_cutoff", "daptomycin_RIS", "daptomycin_DST_method", "study"),
                      sep = "\t")

check1 <- as.data.frame(colnames(known_table_nonsyn))
check2 <- as.data.frame(phenotype$sample_id)
install.packages("sqldf")
library(sqldf)
res <- sqldf('SELECT * FROM check1 EXCEPT SELECT * FROM check2')
print(res)

phenotype <- filter(phenotype, daptomycin_RIS!="NA") ##remove 1 sample with no RIS or MIC data
phenotype <- filter(phenotype, daptomycin_mic!="NA") ##remove samples with no MIC data, unsure of MIC so cannot recode adequately
View(phenotype)

#recode MIC breakpoints. Recode MICs first, as per below. then put as numeric

recode_MIC <- as.vector(phenotype$daptomycin_mic)
recode_MIC <- gsub(">", "", recode_MIC); recode_MIC = gsub("<", "", recode_MIC); recode_MIC = gsub(">=", "", recode_MIC);  recode_MIC = gsub("<=", "", recode_MIC)

tmp <-  which(grepl("<", phenotype$daptomycin_mic)==T); if(length(tmp)>0){ recode_MIC[tmp] = as.character(as.numeric(recode_MIC[tmp])/2); }
tmp <- which(grepl(">", phenotype$daptomycin_mic)==T); if(length(tmp)>0){ recode_MIC[tmp] = as.character(as.numeric(recode_MIC[tmp])*2); }

recode_MIC <- replace(recode_MIC, 374, '1') #sorting troublesome values
recode_MIC <- replace(recode_MIC, 375, '1')
recode_MIC <- replace(recode_MIC, 382, '1')
recode_MIC <- replace(recode_MIC, 383, '1')
recode_MIC <- replace(recode_MIC, 550, '8')
recode_MIC <- replace(recode_MIC, 551, '8')
recode_MIC <- replace(recode_MIC, 552, '8')
recode_MIC <- replace(recode_MIC, 270, '0.5')

recode_MIC <- as.vector(recode_MIC)
phenotype$daptomycin_mic <- recode_MIC
View(phenotype)

phenotype$daptomycin_mic <- as.numeric(phenotype$daptomycin_mic)
phenotype$daptomycin_RIS <- ifelse(phenotype$daptomycin_mic > 8 | phenotype$daptomycin_mic == 8, "NS", "S") 
View(phenotype)

## DEALING WITH DUPLICATE PHENOTYPE DATA

n_occur <- data.frame(table(phenotype$sample_id)) #checking which sample_ids are duplicated
View(n_occur)

phenotype <- phenotype[-c(29,30,31,35),] #removed the E-test entries for these. CHECKED
View(phenotype)
dim(phenotype) #[1] 621 6

# Joining
phenotype_binary <- select(phenotype, sample_id, daptomycin_RIS)

known_table_nonsyn2 <- known_table_nonsyn2[,-(3:4)] #removing blank columns
colnames(known_table_nonsyn2) <- c('sample_id', 'gene', 'mean')
View(known_table_nonsyn2)
known_table_nonsyn2$sample_id <- as.character(known_table_nonsyn2$sample_id) #data wrangling to prepare to join
phenotype_binary$sample_id <- as.character(phenotype_binary$sample_id)
known_nonsyn_phenotype <- left_join(phenotype_binary, known_table_nonsyn2, by="sample_id") #left join to discount all samples without phenotype data
View(known_nonsyn_phenotype)

#Do something to change mean >0 to be binary

known_nonsyn_phenotype$mean <- ifelse(known_nonsyn_phenotype$mean > 0, "1", known_nonsyn_phenotype$mean)
known_nonsyn_phenotype <- arrange(known_nonsyn_phenotype, gene) #arrange by gene
known_nonsyn_phenotype$yes_no <- known_nonsyn_phenotype$mean
known_nonsyn_phenotype <- select(known_nonsyn_phenotype, -mean)

#Recode daptomycin_RIS to be binary (1=NS, 0=S)

known_nonsyn_phenotype$daptomycin_RIS <- ifelse(known_nonsyn_phenotype$daptomycin_RIS == "NS", 1, 0)
View(known_nonsyn_phenotype)

sum(is.na(known_nonsyn_phenotype$yes_no)) #0 NA values

known_nonsyn_phenotype_res <- unique(select(known_nonsyn_phenotype, sample_id, daptomycin_RIS))
as.character(known_nonsyn_phenotype_res$sample_id)

View(known_nonsyn_phenotype_res)
write.table(known_nonsyn_phenotype_res, file = "known_nonsyn_phenotype_res.txt", sep="\t", row.names=F)

#check resistance in clade B samples

clade_b <- read_csv("efm_dap_dataset_v1.cladeB_samples.csv")
clade_b_names <- colnames(clade_b)
clade_b <- rbind(clade_b_names, clade_b)
colnames(clade_b) <- "sample_id"
View(clade_b)
clade_b_vec <- pull(clade_b, sample_id)
View(clade_b_vec)

clade_b_phenotype <- left_join(clade_b, phenotype_binary) 
View(clade_b_phenotype) #31 samples in Clade B, 1/26 with pheno data is NS --> keep analysis with both clades

known_nonsyn_phenotype_res_A <- filter(known_nonsyn_phenotype_res, !(sample_id %in% clade_b_vec))
View(known_nonsyn_phenotype_res_A)

write.table(known_nonsyn_phenotype_res_A, file = "known_nonsyn_phenotype_res_A.txt", sep="\t", row.names=F)



################################################################################################################################

#UNIVARIATE LOGISTIC REGRESSION BY GENE (BINARY)

## Pivot_wider to get gene names across top

known_nonsyn_phenotype1 <- as_tibble(known_nonsyn_phenotype)
wide_known_nonsyn_phenotype <- known_nonsyn_phenotype1 %>%
  pivot_wider(names_from = gene, values_from = yes_no) %>%
  group_by(sample_id)
View(wide_known_nonsyn_phenotype)

## Generate for loop to automate univariate logistic regressions, and present coefficients in tibble form 

library(broom)

uni_log_bin <- list()
uni_log_bin_model <- list()
for(i in 3:ncol(wide_known_nonsyn_phenotype)) {
  uni_log_bin_names <- names(wide_known_nonsyn_phenotype)[[i]]
  uni_log_bin_model[[uni_log_bin_names]] <- glm(daptomycin_RIS ~ wide_known_nonsyn_phenotype[[i]], family="binomial", data=wide_known_nonsyn_phenotype)
  uni_log_bin[[uni_log_bin_names]] <- tidy(uni_log_bin_model[[uni_log_bin_names]])
}

str(uni_log_bin[[uni_log_bin_names]])

length(uni_log_bin) #[1] 16

uni_log_bin_est <- matrix(NA, nrow=length(uni_log_bin), ncol=2)

for(i in 1:length(uni_log_bin)){
  uni_log_bin_est[i,] <- c(uni_log_bin[[i]]$estimate[2], 
                           uni_log_bin[[i]]$p.value[2])
}

View(uni_log_bin_est)
genes <- as_tibble(colnames(wide_known_nonsyn_phenotype))
genes <- genes[-1:-2,] #remove sample/RIS column names
View(genes) #16 genes

uni_log_bin_p <- cbind(genes, uni_log_bin_est)
names(uni_log_bin_p) <- c('Gene', 'Estimate', 'p-value')
View(uni_log_bin_p)

#check by running uni_log_bin - this works

##filter out significant ones (can use on dataframe too)
log_sig_bin_est <-uni_log_bin_p[uni_log_bin_p[,3] <0.05,] 
dim(log_sig_bin_est) #[1] 10 3
View(log_sig_bin_est)
log_sig_genesbin <- as_tibble(log_sig_bin_est[,1])
View(log_sig_genesbin) #the significant genes

##subset significant ones
vectorgene <- pull(log_sig_genesbin, value) #vectorise
wide_known_nonsyn_phenotype_keep <- wide_known_nonsyn_phenotype[,vectorgene]
samplegene_RIS <- select(wide_known_nonsyn_phenotype, sample_id, daptomycin_RIS) #add back the annotation
wide_known_nonsyn_phenotype_keep <- cbind(samplegene_RIS, wide_known_nonsyn_phenotype_keep)

View(wide_known_nonsyn_phenotype_keep) #621x10

write.table(log_sig_bin_est, file = "log_sig_bin_est.txt", sep=",")

####################################################################################

##MULTIVARIATE LOGISTIC REGRESSION by gene binary - MODEL SELECTION

library(pscl)

wide_known_nonsyn_phenotype_keep1 <- wide_known_nonsyn_phenotype_keep[,-1] #remove sample name for regression
View(wide_known_nonsyn_phenotype_keep1)
multi_log_bin1 <- glm(daptomycin_RIS~., data = wide_known_nonsyn_phenotype_keep1, family="binomial")
summary(multi_log_bin1)

pR2(multi_log_bin1)['McFadden']

multi_log_bin2 <- step(multi_log_bin1, direction="both")
summary(multi_log_bin2) 

pR2(multi_log_bin2)['McFadden']

multi_log_bin3 <- step(multi_log_bin1, direction="backward")
summary(multi_log_bin3) 

multi_log_bin_null <- glm(daptomycin_RIS~1, data = wide_known_nonsyn_phenotype_keep1)
multi_log_bin4 <- step(multi_log_bin_null, direction="forward", scope=formula(multi_log_bin1), trace=0)
summary(multi_log_bin4)


########################################################################################################

##POPN STRUCTURE - CODE in "phylogenetic_clustering - Clade A.R" and "- Clade B.R"

#check resistance in clade B samples

clade_b <- read_csv("efm_dap_dataset_v1.cladeB_samples.csv")
clade_b_names <- colnames(clade_b)
clade_b <- rbind(clade_b_names, clade_b)
colnames(clade_b) <- "sample_id"

clade_b_phenotype <- left_join(clade_b, phenotype_binary) 
View(clade_b_phenotype) #31 samples in Clade B, 1/26 with pheno data is NS --> keep analysis with both clades

##Clade A

View(lmTable)

lmTable <- as_tibble(lmTable)
lmTable$sample_id <- lmTable$Taxa

#recode all singletons --> singleton
lmTable$Cluster.DSF <- as_tibble(ifelse(grepl('^singleton',lmTable$Cluster.DSF), 'singleton', lmTable$Cluster.DSF))
View(lmTable)

#recode sample_ids back
phylo_sample <- as.data.frame(lmTable$sample_id)
phylo_sample$sample_id2 <- gsub("_NA#NA$", "", phylo_sample$`lmTable$sample_id`)
phylo_sample$sample_id3 <- gsub("#NA$", "", phylo_sample$sample_id2)
View(phylo_sample)
phylo_sample2 <- as_tibble(phylo_sample[,3])
View(phylo_sample2)

#join to cluster data
phylo_table <- cbind(phylo_sample2, lmTable)
phylo_table <- phylo_table[,-c(2,4)]
names(phylo_table) <- c("sample_id", "Cluster.DSF")
View(phylo_table)

dim(phylo_table)
# [1] 639 2 -- i.e. before removing any samples without phenotype data

##Clade B 

View(lmTableB)

lmTableB1 <- as_tibble(lmTableB)
lmTableB1$sample_id <- lmTableB1$Taxa

#recode all singletons --> singleton
lmTableB1$Cluster.DSF <- ifelse(grepl('^singleton',lmTableB1$Cluster.DSF), 'singleton', lmTableB1$Cluster.DSF)
View(lmTableB1)

#recode sample_ids back
phylo_sample_b <- as.data.frame(lmTableB1$sample_id)
phylo_sample_b$sample_id2 <- gsub("_NA#NA$", "", phylo_sample_b$`lmTableB1$sample_id`)
phylo_sample_b$sample_id3 <- gsub("#NA$", "", phylo_sample_b$sample_id2)
phylo_sample_b2 <- as_tibble(phylo_sample_b[,3])
View(phylo_sample_b2)

#join to cluster data
phylo_table_b <- cbind(phylo_sample_b2, lmTableB1)
phylo_table_b <- phylo_table_b[,-c(2,4)]
names(phylo_table_b) <- c("sample_id", "Cluster.DSF")
View(phylo_table_b)

dim(phylo_table_b)
# [1] 30 2 -- again before removal of samples without phenotype data

#recode B samples to have B in front of them

class(phylo_table_b$Cluster.DSF)

phylo_table_b$Cluster.DSF <- paste0("b_", phylo_table_b$Cluster.DSF)
View(phylo_table_b)

##Combining the two

complete_phylo <- rbind(phylo_table, phylo_table_b)
View(complete_phylo)

##need to join to rest of data - several versions


#####################################################################################

## UNIVARIATE LOGISTIC REGRESSION BY GENE (BINARY) - POPN STRUCTURE

uni_log_bin_phylo <- inner_join(complete_phylo, wide_known_nonsyn_phenotype, join_by= "sample_id")
dim(uni_log_bin_phylo) #[1] 620 19
View(uni_log_bin_phylo) 
sum(is.na(uni_log_bin_phylo)) #0

uni_log_bin_phylo1 <- uni_log_bin_phylo[,-1] #remove sample IDs for regressions
uni_log_bin_phylo1$Cluster.DSF <- pull(uni_log_bin_phylo1$Cluster.DSF) #vectorise
View(uni_log_bin_phylo1)
pop_uni_log_bin_model <- list()
pop_uni_log_bin <- list()
for(i in 3:ncol(uni_log_bin_phylo1)) {
  pop_uni_log_bin_names <- names(uni_log_bin_phylo1)[[i]]
  pop_uni_log_bin_model[[pop_uni_log_bin_names]] <- glm(daptomycin_RIS ~ uni_log_bin_phylo1[[i]] + Cluster.DSF, family="binomial", data=uni_log_bin_phylo1)
  pop_uni_log_bin[[pop_uni_log_bin_names]] <- tidy(pop_uni_log_bin_model[[pop_uni_log_bin_names]])
}

head(pop_uni_log_bin) 
str(pop_uni_log_bin[[pop_uni_log_bin_names]])

length(pop_uni_log_bin)

pop_uni_log_bin_est <- matrix(NA, nrow=length(pop_uni_log_bin), ncol=2)

for(i in 1:length(pop_uni_log_bin)){
  pop_uni_log_bin_est[i,] <- c(pop_uni_log_bin[[i]]$estimate[2], 
                           pop_uni_log_bin[[i]]$p.value[2])
}

View(pop_uni_log_bin_est)

pop_uni_log_bin_p <- cbind(genes, pop_uni_log_bin_est)
names(pop_uni_log_bin_p) <- c('Gene', 'Estimate', 'p-value')
View(pop_uni_log_bin_p)

##filter out significant ones (can use on dataframe too)
pop_log_sig_bin_est <-pop_uni_log_bin_p[pop_uni_log_bin_p[,3] <0.05,] ##5 significant genes
View(pop_log_sig_bin_est)
pop_log_sig_genesbin <- as_tibble(pop_log_sig_bin_est[,1]) ##the genes which are significant
View(pop_log_sig_genesbin)

##subset significant ones
pop_vectorgene <- pull(pop_log_sig_genesbin, value)
uni_log_bin_phylo_keep <- uni_log_bin_phylo1[,pop_vectorgene]
pop_samplegene_RIS <- select(uni_log_bin_phylo, sample_id, daptomycin_RIS, Cluster.DSF)
uni_log_bin_phylo_keep <- cbind(pop_samplegene_RIS, uni_log_bin_phylo_keep)
View(uni_log_bin_phylo_keep)

write.table(pop_log_sig_bin_est, file = "pop_log_sig_bin_est.txt", sep=",")

################################################################################################

## MULTIVARIATE LOGISTIC REGRESSION by gene (binary) - POPN STRUCTURE

uni_log_bin_phylo_keep1 <- uni_log_bin_phylo_keep[,-1] 
uni_log_bin_phylo_keep1$Cluster.DSF <- pull(uni_log_bin_phylo_keep1$Cluster.DSF)

multi_log_bin_phylo1 <- glm(daptomycin_RIS~., data = uni_log_bin_phylo_keep1, family="binomial")
summary(multi_log_bin_phylo1)

pR2(multi_log_bin_phylo1)['McFadden']

multi_log_bin_phylo2 <- step(multi_log_bin_phylo1, direction="both")
summary(multi_log_bin_phylo2) 

pR2(multi_log_bin_phylo2)['McFadden']

##########################################################################################

## DATA WRANGLING - GENE COUNT (number of variants per gene)

known_table_nonsyn1$count <- ifelse(known_table_nonsyn1$count == 2, "1", "0")
View(known_table_nonsyn1)

known_table_nonsyn1$count <- as.numeric(known_table_nonsyn1$count)

known_table_nonsyn3 <- aggregate(known_table_nonsyn1$count, by=list(known_table_nonsyn1$sample_id, known_table_nonsyn1$gene), FUN=sum)
View(known_table_nonsyn3)
colnames(known_table_nonsyn3) <- c('sample_id', 'gene', 'number')
known_table_nonsyn3$sample_id <- as.character(known_table_nonsyn3$sample_id) #data wrangling to prepare to join

phenotype_binary$sample_id <- as.character(phenotype_binary$sample_id)
known_nonsyn_phenotype2 <- left_join(phenotype_binary, known_table_nonsyn3, by="sample_id")
dim(known_nonsyn_phenotype2) #[1] 9936 4
View(known_nonsyn_phenotype2)

#Recode daptomycin_RIS to be binary (1=NS, 0=S)

known_nonsyn_phenotype2$daptomycin_RIS <- ifelse(known_nonsyn_phenotype2$daptomycin_RIS == "NS", 1, 0)
View(known_nonsyn_phenotype2) 

#############################################################################################

## UNIVARIATE LOGISTIC REGRESSION - GENE COUNT

## Pivot_wider to get gene names across top

known_nonsyn_phenotype2 <- as_tibble(known_nonsyn_phenotype2)
wide_known_nonsyn_phenotype2 <- known_nonsyn_phenotype2 %>%
  pivot_wider(names_from = gene, values_from = number) %>%
  group_by(sample_id)
View(wide_known_nonsyn_phenotype2) #621x16

## Generate for loop to automate univariate logistic regressions, and present coefficients in tibble form 

library(broom)

uni_log_count <- list()
uni_log_count_model <- list()
for(i in 3:ncol(wide_known_nonsyn_phenotype2)) {
  uni_log_count_names <- names(wide_known_nonsyn_phenotype2)[[i]]
  uni_log_count_model[[uni_log_count_names]] <- glm(daptomycin_RIS ~ wide_known_nonsyn_phenotype2[[i]], family="binomial", data=wide_known_nonsyn_phenotype2)
  uni_log_count[[uni_log_count_names]] <- tidy(uni_log_count_model[[uni_log_count_names]])
}

str(uni_log_count[[uni_log_count_names]])

length(uni_log_count) #[1] 16

uni_log_count_est <- matrix(NA, nrow=length(uni_log_count), ncol=2)

for(i in 1:length(uni_log_count)){
  uni_log_count_est[i,] <- c(uni_log_count[[i]]$estimate[2], uni_log_count[[i]]$p.value[2])
}

View(uni_log_count_est)
uni_log_count_p <- cbind(genes, uni_log_count_est)
names(uni_log_count_p) <- c('Gene', 'Estimate', 'p-value')
View(uni_log_count_p) #CHECKED

##filter out significant ones (can use on dataframe too)
log_sig_count_est <- uni_log_count_p[uni_log_count_p[,3] <0.05,] ##14 significant genes
View(log_sig_count_est)
log_sig_genescount <- as_tibble(log_sig_count_est[,1])
View(log_sig_genescount)

##subset significant ones
vectorgenec <- pull(log_sig_genescount, value)
wide_known_nonsyn_phenotype2_keep <- wide_known_nonsyn_phenotype2[,vectorgenec]
samplecount_RIS <- select(wide_known_nonsyn_phenotype2, sample_id, daptomycin_RIS)
wide_known_nonsyn_phenotype2_keep <- cbind(samplecount_RIS, wide_known_nonsyn_phenotype2_keep)
View(wide_known_nonsyn_phenotype2_keep)
write.table(log_sig_count_est, file = "log_sig_count_est.txt", sep=",")

#############################################################################################

##MULTIVARIATE LOGISTIC REGRESSION - GENE COUNT - MODEL SELECTION

wide_known_nonsyn_phenotype2_keep <- wide_known_nonsyn_phenotype2_keep[,-1]
multi_log_count1 <- glm(daptomycin_RIS~., data = wide_known_nonsyn_phenotype2_keep, family="binomial")
summary(multi_log_count1)

pR2(multi_log_count1)['McFadden']

multi_log_count2 <- step(multi_log_count1, direction="both")
summary(multi_log_count2) #AIC 

pR2(multi_log_count2)['McFadden']

multi_log_count3 <- step(multi_log_count1, direction="backward")
summary(multi_log_count3) #AIC 

multi_log_count_null <- glm(daptomycin_RIS~1, data = wide_known_nonsyn_phenotype2_keep)
multi_log_count4 <- step(multi_log_count_null, direction="forward", scope=formula(multi_log_count1), trace=0)
summary(multi_log_count4) #AIC 

################################################################################################

##UNIVARIATE LOGISTIC REGRESSION - GENE COUNT - POPN STRUCTURE

uni_log_count_phylo <- inner_join(complete_phylo, wide_known_nonsyn_phenotype2, join_by= "sample_id")
uni_log_count_phylo1 <- uni_log_count_phylo[,-1]
dim(uni_log_bin_phylo1) #[1] 620 18
View(uni_log_count_phylo1)
uni_log_count_phylo1$Cluster.DSF <- pull(uni_log_count_phylo1$Cluster.DSF)

pop_uni_log_count_model <- list()
pop_uni_log_count <- list()
for(i in 3:ncol(uni_log_count_phylo1)) {
  pop_uni_log_count_names <- names(uni_log_count_phylo1)[[i]]
  pop_uni_log_count_model[[pop_uni_log_count_names]] <- glm(daptomycin_RIS ~ uni_log_count_phylo1[[i]] + Cluster.DSF, family="binomial", data=uni_log_count_phylo1)
  pop_uni_log_count[[pop_uni_log_count_names]] <- tidy(pop_uni_log_count_model[[pop_uni_log_count_names]])
}

head(pop_uni_log_count)

str(pop_uni_log_count[[pop_uni_log_count_names]])

length(pop_uni_log_count) #[1] 16

pop_uni_log_count_est <- matrix(NA, nrow=length(pop_uni_log_count), ncol=2)

for(i in 1:length(pop_uni_log_count)){
  pop_uni_log_count_est[i,] <- c(pop_uni_log_count[[i]]$estimate[2], 
                         pop_uni_log_count[[i]]$p.value[2])
}

View(pop_uni_log_count_est)

pop_uni_log_count_p <- cbind(genes, pop_uni_log_count_est) 
names(pop_uni_log_count_p) <- c('gene', 'estimate', 'p-value')
View(pop_uni_log_count_p)

##filter out significant ones (can use on dataframe too)
pop_log_sig_count_est <- pop_uni_log_count_p[pop_uni_log_count_p[,3] <0.05,] ##5 significant genes
dim(pop_log_sig_count_est) #[1] 5 3
View(pop_log_sig_count_est)
pop_log_sig_genescount <- as_tibble(pop_log_sig_count_est[,1])
View(pop_log_sig_genescount)

##subset significant ones
pop_vectorgenec <- pull(pop_log_sig_genescount, value)
uni_log_count_phylo_keep <- uni_log_count_phylo1[,pop_vectorgenec]
pop_samplecount_RIS <- select(uni_log_count_phylo, sample_id, daptomycin_RIS, Cluster.DSF)
uni_log_count_phylo_keep <- cbind(pop_samplecount_RIS, uni_log_count_phylo_keep)

write.table(pop_log_sig_count_est, file = "pop_log_sig_count_est.txt", sep=",")

#################################################################################################

##MULTIVARIATE LOGISTIC REGRESSION - GENE COUNT - POPN STRUCTURE

uni_log_count_phylo_keep1 <- uni_log_count_phylo_keep[,-1] 
uni_log_count_phylo_keep1$Cluster.DSF <- pull(uni_log_count_phylo_keep1$Cluster.DSF)
View(uni_log_count_phylo_keep1)

multi_log_count_phylo1 <- glm(daptomycin_RIS~., data = uni_log_count_phylo_keep1, family="binomial")
summary(multi_log_count_phylo1)

pR2(multi_log_count_phylo1)['McFadden']

multi_log_count_phylo2 <- step(multi_log_count_phylo1, direction="both")
summary(multi_log_count_phylo2) 

pR2(multi_log_count_phylo2)['McFadden']


#################################################################################################

# DATA WRANGLING - LOGISTIC REGRESSION BY VARIANT - start with known_full_table

full_variant <- rownames_to_column(known_full_table)
transp_full_variant <- t(full_variant)
transp_full_variant <- as.data.frame(transp_full_variant)
transp_full_variant <- rownames_to_column(transp_full_variant)
colnames(transp_full_variant) <- transp_full_variant[1,]
transp_full_variant <- transp_full_variant[-1,]
View(transp_full_variant)
transp_full_variant <- as_tibble(transp_full_variant)
transp_variant_nonsyn <- filter(transp_full_variant, severity!="LOW")
transp_variant_nonsyn <- transp_variant_nonsyn[,-(2:5)]
View(transp_variant_nonsyn) #548 variants

transp_variant_nonsyn <- transp_variant_nonsyn[-167,] #removing variant in 02872 gene again
View(transp_variant_nonsyn) #547x670
variant_nonsyn <- t(transp_variant_nonsyn)
  
colnames(variant_nonsyn) <- variant_nonsyn[1,]
variant_nonsyn <- variant_nonsyn %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(sample_id = rowname)
View(variant_nonsyn)
variant_nonsyn <- variant_nonsyn[-1,]
variant_nonsyn <- variant_nonsyn[-671,]
View(variant_nonsyn) #670x547 variants plus sample_id column

##recode 2 to 1 (proper binary)

variant_nonsyn[variant_nonsyn == "2"] <- "1"
View(variant_nonsyn)

## add phenotype_binary
variant_nonsyn$sample_id <- as.character(variant_nonsyn$sample_id) #data wrangling to prepare to join
phenotype_binary$sample_id <- as.character(phenotype_binary$sample_id)
variant_nonsyn_pheno <- left_join(phenotype_binary, variant_nonsyn, by="sample_id")
View(variant_nonsyn_pheno)
dim(variant_nonsyn_pheno) #[1] 621 549
sum(is.na(variant_nonsyn_pheno)) #[1] 0

##Recode daptomycin_RIS to be binary (1=NS, 0=S)

variant_nonsyn_pheno$daptomycin_RIS <- ifelse(variant_nonsyn_pheno$daptomycin_RIS == "NS", 1, 0)
View(variant_nonsyn_pheno)

head(sapply(lapply(variant_nonsyn_pheno, unique), length)) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these

uniquelength <- sapply(variant_nonsyn_pheno,function(x) length(unique(x)))
variant_nonsyn_pheno <- subset(variant_nonsyn_pheno, select=uniquelength>1)
View(variant_nonsyn_pheno) #621 samples, 525 variants

####################################################################################

##UNIVARIATE LOGISTIC REGRESSION - VARIANT/BINARY

library(broom)

uni_log_var <- list()
uni_log_var_model <- list()
for(i in 3:ncol(variant_nonsyn_pheno)) {
  uni_log_var_names <- names(variant_nonsyn_pheno)[[i]]
  uni_log_var_model[[uni_log_var_names]] <- glm(daptomycin_RIS ~ variant_nonsyn_pheno[[i]], family="binomial", data=variant_nonsyn_pheno)
  uni_log_var[[uni_log_var_names]] <- tidy(uni_log_var_model[[uni_log_var_names]])
    }

str(uni_log_var[[uni_log_var_names]])

length(uni_log_var)

uni_log_var_est <- matrix(NA, nrow=length(uni_log_var), ncol=2)

for(i in 1:length(uni_log_var)){
  uni_log_var_est[i,] <- c(uni_log_var[[i]]$estimate[2], 
                           uni_log_var[[i]]$p.value[2])
}

dim(uni_log_var_est) #[1] 525 2
variants <- as_tibble(colnames(variant_nonsyn_pheno))
variants <- variants[-1:-2,]
View(variants)

uni_log_var_p <- cbind(variants, uni_log_var_est)
names(uni_log_var_p) <- c('Variant', 'Estimate', 'p-value')
View(uni_log_var_p)

##filter out significant ones (can use on dataframe too)
uni_log_var_p[,3] < 0.05 ##returns true or false
log_sig_var_est <-uni_log_var_p[uni_log_var_p[,3] <0.05,] ##82 significant variants
View(log_sig_var_est)
log_sig_variants <- as_tibble(log_sig_var_est[,1])
View(log_sig_variants)

##subset variant_nonsyn_pheno
vector <- pull(log_sig_variants, value)
variant_nonsyn_pheno1 <- variant_nonsyn_pheno[,vector, with=FALSE]
View(variant_nonsyn_pheno1)
sample_RIS <- select(variant_nonsyn_pheno, sample_id, daptomycin_RIS)
variant_nonsyn_pheno1 <- cbind(sample_RIS, variant_nonsyn_pheno1)
View(variant_nonsyn_pheno1)

write.table(log_sig_var_est, file = "log_sig_var_est.txt", sep=",")

#######################################################################################

##MULTIVARIATE LOGISTIC REGRESSION  - VARIANT - MODEL SELECTION

variant_nonsyn_pheno1 <- variant_nonsyn_pheno1[,-1] 
View(variant_nonsyn_pheno1)
multi_log_variant1 <- glm(daptomycin_RIS~., data = variant_nonsyn_pheno1, family="binomial") ##some predictors have perfect correlation/linear relationship?
summary(multi_log_variant1)

pR2(multi_log_variant1)['McFadden']

multi_log_variant2 <- step(multi_log_variant1, direction="both")
summary(multi_log_variant2) #AIC 

pR2(multi_log_variant2)['McFadden']

multi_log_variant3 <- step(multi_log_variant1, direction="backward")
summary(multi_log_variant3) #AIC 

multi_log_variant_null <- glm(daptomycin_RIS~1, data = variant_nonsyn_pheno1)
multi_log_variant4 <- step(multi_log_variant_null, direction="forward", scope=formula(multi_log_variant1), trace=0)
summary(multi_log_variant4) #AIC 

############################################################################################

##UNIVARIATE LOGISTIC REGRESSION - VARIANT - POPN STRUCTURE

uni_log_var_phylo <- inner_join(complete_phylo, variant_nonsyn_pheno, join_by= "sample_id")
View(uni_log_var_phylo) #620x525+cluster

clustervector <- uni_log_var_phylo$Cluster.DSF

sapply(lapply(uni_log_var_phylo, unique), length) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these
uniquelength2 <- sapply(uni_log_var_phylo,function(x) length(unique(x)))
uni_log_var_phylo <- subset(uni_log_var_phylo, select=uniquelength2>1)
dim(uni_log_var_phylo) #[1] 620 521 variants. loses Cluster column so add this back in

uni_log_var_phylo <- cbind(clustervector, uni_log_var_phylo)
View(uni_log_var_phylo)
uni_log_var_phylo <- rename(uni_log_var_phylo, Cluster.DSF = value)
uni_log_var_phylo1 <- uni_log_var_phylo[,-2]  
View(uni_log_var_phylo)

pop_uni_log_var <- list()
pop_uni_log_var_model <- list()
for(i in 3:ncol(uni_log_var_phylo1)) {
  pop_uni_log_var_names <- names(uni_log_var_phylo1)[[i]]
  pop_uni_log_var_model[[pop_uni_log_var_names]] <- glm(daptomycin_RIS ~ uni_log_var_phylo1[[i]]+ Cluster.DSF, family="binomial", data=uni_log_var_phylo1)
  pop_uni_log_var[[pop_uni_log_var_names]] <- tidy(pop_uni_log_var_model[[pop_uni_log_var_names]])
}

head(pop_uni_log_var)

str(pop_uni_log_var[[pop_uni_log_var_names]])

length(pop_uni_log_var)#520 - 519 variants + 1 cluster term

pop_uni_log_var_est <- matrix(NA, nrow=length(pop_uni_log_var), ncol=2)

for(i in 1:length(pop_uni_log_var)){
  pop_uni_log_var_est[i,] <- c(pop_uni_log_var[[i]]$estimate[2], 
                           pop_uni_log_var[[i]]$p.value[2])
}

View(pop_uni_log_var_est)
variants2 <- as_tibble(colnames(uni_log_var_phylo1))
variants2 <- variants2[-1:-2,]
View(variants2)

pop_uni_log_var_p <- cbind(variants2, pop_uni_log_var_est)
names(pop_uni_log_var_p) <- c('Variant', 'Estimate', 'p-value')
View(pop_uni_log_var_p)

##filter out significant variants
pop_uni_log_var_p[,3] < 0.05 ##returns true or false
pop_log_sig_var_est <- pop_uni_log_var_p[pop_uni_log_var_p[,3] <0.05,] ##32 significant variants
View(pop_log_sig_var_est)
write.table(pop_log_sig_var_est, file = "pop_log_sig_var_est.txt", sep=",")

pop_log_sig_variants <- as_tibble(pop_log_sig_var_est[,1])
View(pop_log_sig_variants)

##subset significant ones
vector_pop_var <- pull(pop_log_sig_variants, value)
uni_log_var_phylo_keep <- uni_log_var_phylo[,vector_pop_var]
View(uni_log_var_phylo_keep)
cluster_RIS <- select(uni_log_var_phylo, sample_id, daptomycin_RIS, Cluster.DSF)
uni_log_var_phylo_keep <- cbind(cluster_RIS, uni_log_var_phylo_keep)
View(uni_log_var_phylo_keep)


#############################################################################################

##MULTIVARIATE LOGISTIC REGRESSION - VARIANT - POPN STRUCTURE

uni_log_var_phylo_keep1 <- uni_log_var_phylo_keep[,-1] #remove sample_id
View(uni_log_var_phylo_keep1)

multi_log_var_phylo1 <- glm(daptomycin_RIS~., data = uni_log_var_phylo_keep1, family="binomial") ##some exact fitting of probabilities?
summary(multi_log_var_phylo1)

pR2(multi_log_var_phylo1)['McFadden']

multi_log_var_phylo2 <- step(multi_log_var_phylo1, direction="both")
summary(multi_log_var_phylo2) 

pR2(multi_log_var_phylo2)['McFadden']


#############################################################################################
## PART 2 - MIC (linear regressions)
#############################################################################################

## DATA WRANGLING - MIC

phenotype_MIC <- select(phenotype, sample_id, daptomycin_mic)
phenotype_MIC$trans_MIC <- log2(phenotype_MIC$daptomycin_mic)
phenotype_MIC <- select(phenotype_MIC, -daptomycin_mic)
View(phenotype_MIC)

###############################################################################################

## UNIVARIATE LINEAR REGRESSION - GENE/BINARY

View(known_table_nonsyn2)
known_table_nonsyn2$sample_id <- as.character(known_table_nonsyn2$sample_id)
MIC_full_table <- left_join(phenotype_MIC, known_table_nonsyn2)
sum(is.na(MIC_full_table)) #[1] 0
View(MIC_full_table) #9936 entries

MIC_full_table$mean <- ifelse(MIC_full_table$mean > 0, "1", MIC_full_table$mean)
MIC_full_table <- arrange(MIC_full_table, gene) #arrange by gene
MIC_full_table$yes_no <- MIC_full_table$mean
MIC_full_table <- select(MIC_full_table, -mean)
View(MIC_full_table)

pheno_linear <- unique(select(MIC_full_table, sample_id, trans_MIC))
as.character(pheno_linear$sample_id)

View(pheno_linear)
write.table(pheno_linear, file = "pheno_linear.txt", sep="\t", row.names=F)

pheno_linear_A <- filter(pheno_linear, !(sample_id %in% clade_b_vec))
View(pheno_linear_A)
as.character(pheno_linear_A$sample_id)
write.table(pheno_linear_A, file="pheno_linear_A.txt", sep="\t", row.names = F)

## Pivot_wider to get gene names across top

MIC_full_table1 <- as_tibble(MIC_full_table)
wide_MIC_full_table <- MIC_full_table1 %>%
  pivot_wider(names_from = gene, values_from = yes_no) %>%
  group_by(sample_id)
dim(wide_MIC_full_table) #[1] 621 18
View(wide_MIC_full_table)

## Generate for loop to automate univariate linear regressions, and present coefficients in tibble form 

library(broom)
wide_MIC_full_table$trans_MIC <- as.numeric(wide_MIC_full_table$trans_MIC)
uni_lin_bin_model <- list()
uni_lin_bin <- list()
for(i in 3:ncol(wide_MIC_full_table)) {
  uni_lin_bin_names <- names(wide_MIC_full_table)[[i]]
  uni_lin_bin_model[[uni_lin_bin_names]] <- glm(trans_MIC ~ wide_MIC_full_table[[i]], data=wide_MIC_full_table)
  uni_lin_bin[[uni_lin_bin_names]] <- tidy(uni_lin_bin_model[[uni_lin_bin_names]])
}

str(uni_lin_bin[[uni_lin_bin_names]])

length(uni_lin_bin) #[1] 16

uni_lin_bin_est <- matrix(NA, nrow=length(uni_lin_bin), ncol=2)

for(i in 1:length(uni_lin_bin)){
  uni_lin_bin_est[i,] <- c(uni_lin_bin[[i]]$estimate[2], 
                    uni_lin_bin[[i]]$p.value[2])
}

View(uni_lin_bin_est)

uni_lin_bin_p <- cbind(genes, uni_lin_bin_est)
names(uni_lin_bin_p) <- c('gene', 'Estimate', 'p-value')
View(uni_lin_bin_p)

##filter out significant ones (can use on dataframe too)
lin_sig_bin_est <-uni_lin_bin_p[uni_lin_bin_p[,3] <0.05,] ##10 significant genes
dim(lin_sig_bin_est) #[1] 10 3
View(lin_sig_bin_est)
lin_sig_genesbin <- as_tibble(lin_sig_bin_est[,1])
View(lin_sig_genesbin)

##subset significant ones
lin_vectorgene <- pull(lin_sig_genesbin, value)
wide_MIC_full_table_keep <- wide_MIC_full_table[,lin_vectorgene]
samplebin_MIC <- select(wide_MIC_full_table, sample_id, trans_MIC)
wide_MIC_full_table_keep <- cbind(samplebin_MIC, wide_MIC_full_table_keep)
View(wide_MIC_full_table_keep)
write.table(lin_sig_bin_est, file = "lin_sig_bin_est.txt", sep=",")


##########################################################################################################

## MULTIVARIATE LINEAR REGRESSION - GENE/BINARY - MODEL SELECTION

wide_MIC_full_table_keep1 <- wide_MIC_full_table_keep[,-1] #remove sample ID column

multi_lin_bin1 <- lm(trans_MIC~., data = wide_MIC_full_table_keep1) ##full model
summary(multi_lin_bin1)


multi_lin_bin2 <- step(multi_lin_bin1, direction="both")
summary(multi_lin_bin2) #AIC - R2 values?


multi_lin_bin3 <- step(multi_lin_bin1, direction="backward")
summary(multi_lin_bin3) #AIC

multi_lin_bin_null <- glm(trans_MIC~1, data = wide_MIC_full_table1)
multi_lin_bin4 <- step(multi_lin_bin_null, direction="forward", scope=formula(multi_lin_bin1), trace=0)
summary(multi_lin_bin4) #AIC

################################################################################################################

## UNIVARIATE LINEAR REGRESSION - GENE/BINARY - POPN STRUCTURE

uni_lin_bin_phylo <- inner_join(complete_phylo, wide_MIC_full_table, join_by= "sample_id")
dim(uni_lin_bin_phylo) #[1] 620 19

uni_lin_bin_phylo1 <- uni_lin_bin_phylo[,-1]
uni_lin_bin_phylo1$Cluster.DSF <- pull(uni_lin_bin_phylo1$Cluster.DSF)
View(uni_lin_bin_phylo1)

pop_uni_lin_bin_model <- list()
pop_uni_lin_bin <- list()
for(i in 3:ncol(uni_lin_bin_phylo1)) {
  pop_uni_lin_bin_names <- names(uni_lin_bin_phylo1)[[i]]
  pop_uni_lin_bin_model[[pop_uni_lin_bin_names]] <- glm(trans_MIC ~ uni_lin_bin_phylo1[[i]] + Cluster.DSF, data=uni_lin_bin_phylo1)
  pop_uni_lin_bin[[pop_uni_lin_bin_names]] <- tidy(pop_uni_lin_bin_model[[pop_uni_lin_bin_names]])
}

pop_uni_lin_bin

str(pop_uni_lin_bin[[pop_uni_lin_bin_names]])

length(pop_uni_lin_bin) #[1] 16

pop_uni_lin_bin_est <- matrix(NA, nrow=length(pop_uni_lin_bin), ncol=2)

for(i in 1:length(pop_uni_lin_bin)){
  pop_uni_lin_bin_est[i,] <- c(pop_uni_lin_bin[[i]]$estimate[2], 
                       pop_uni_lin_bin[[i]]$p.value[2])
}

View(pop_uni_lin_bin_est)

pop_uni_lin_bin_p <- cbind(genes, pop_uni_lin_bin_est)
names(pop_uni_lin_bin_p) <- c('gene', 'Estimate', 'p-value')
View(pop_uni_lin_bin_p)

##filter out significant ones (can use on dataframe too)
pop_lin_sig_bin_est <- pop_uni_lin_bin_p[pop_uni_lin_bin_p[,3] <0.05,] 
dim(pop_lin_sig_bin_est) #[1] 7 3
View(pop_lin_sig_bin_est)
pop_lin_sig_genesbin <- as_tibble(pop_lin_sig_bin_est[,1])
View(pop_lin_sig_genesbin)

##subset significant ones
pop_lin_vectorgene <- pull(pop_lin_sig_genesbin, value)
uni_lin_bin_phylo_keep <- uni_lin_bin_phylo1[,pop_lin_vectorgene]
pop_samplebin_MIC <- select(uni_lin_bin_phylo, sample_id, trans_MIC, Cluster.DSF)
uni_lin_bin_phylo_keep <- cbind(pop_samplebin_MIC, uni_lin_bin_phylo_keep)
View(uni_lin_bin_phylo_keep)

write.table(pop_lin_sig_bin_est, file = "pop_lin_sig_bin_est.txt", sep=",")

################################################################################################################

## MULTIVARIATE LINEAR REGRESSION - GENE/BINARY - POPN STRUCTURE

uni_lin_bin_phylo_keep1 <- uni_lin_bin_phylo_keep[,-1] 
uni_lin_bin_phylo_keep1$Cluster.DSF <- pull(uni_lin_bin_phylo_keep1$Cluster.DSF)

multi_lin_bin_phylo1 <- glm(trans_MIC~., data = uni_lin_bin_phylo_keep1)
summary(multi_lin_bin_phylo1)

multi_lin_bin_phylo2 <- step(multi_lin_bin_phylo1, direction="both")
summary(multi_lin_bin_phylo2) 

################################################################################################################

## UNIVARIATE LINEAR REGRESSION (MIC) - COUNT

View(known_table_nonsyn3)

known_table_nonsyn3$sample_id <- as.character(known_table_nonsyn3$sample_id) #data wrangling to prepare to join
phenotype_MIC$sample_id <- as.character(phenotype_MIC$sample_id)
known_nonsyn_phenotype3 <- left_join(phenotype_MIC, known_table_nonsyn3, by="sample_id")
View(known_nonsyn_phenotype3) #9936 entries

## Pivot_wider to get gene names across top

known_nonsyn_phenotype4 <- distinct(known_nonsyn_phenotype3)
known_nonsyn_phenotype4 <- as_tibble(known_nonsyn_phenotype4)
View(known_nonsyn_phenotype4)
wide_known_nonsyn_phenotype3 <- known_nonsyn_phenotype4 %>%
  pivot_wider(names_from = gene, values_from = number) %>%
  group_by(sample_id)
View(wide_known_nonsyn_phenotype3) #621 by 16 genes

## Generate for loop to automate univariate linear regressions, and present coefficients in tibble form 

library(broom)

uni_lin_count_model <- list()
uni_lin_count <- list()
for(i in 3:ncol(wide_known_nonsyn_phenotype3)) {
  uni_lin_count_names <- names(wide_known_nonsyn_phenotype3)[[i]]
  uni_lin_count_model[[uni_lin_count_names]] <- glm(trans_MIC ~ wide_known_nonsyn_phenotype3[[i]], data=wide_known_nonsyn_phenotype3)
  uni_lin_count[[uni_lin_count_names]] <- tidy(uni_lin_count_model[[uni_lin_count_names]])
}

uni_lin_count

str(uni_lin_count[[uni_lin_count_names]])

length(uni_lin_count) #[1] 16

uni_lin_count_est <- matrix(NA, nrow=length(uni_lin_count), ncol=2)

for(i in 1:length(uni_lin_count)){
  uni_lin_count_est[i,] <- c(uni_lin_count[[i]]$estimate[2], 
                         uni_lin_count[[i]]$p.value[2])
}

View(uni_lin_count_est)

uni_lin_count_p <- cbind(genes, uni_lin_count_est)
names(uni_lin_count_p) <- c('gene', 'Estimate', 'p-value')
View(uni_lin_count_p)

##filter out significant ones (can use on dataframe too)
lin_sig_count_est <-uni_lin_count_p[uni_lin_count_p[,3] <0.05,] 
dim(lin_sig_count_est) #[1] 9 3
View(lin_sig_count_est)
lin_sig_genescount <- as_tibble(lin_sig_count_est[,1])
View(lin_sig_genescount)

##subset significant ones
lin_vectorgenec <- pull(lin_sig_genescount, value)
wide_known_nonsyn_phenotype3_keep <- wide_known_nonsyn_phenotype3[,lin_vectorgenec]
samplecount_MIC <- select(wide_known_nonsyn_phenotype3, sample_id, trans_MIC)
wide_known_nonsyn_phenotype3_keep <- cbind(samplecount_MIC, wide_known_nonsyn_phenotype3_keep)
View(wide_known_nonsyn_phenotype3_keep)
write.table(lin_sig_count_est, file = "lin_sig_count_est.txt", sep=",")


###########################################################################################

## MULTIVARIATE LINEAR REGRESSION - GENE COUNT - MODEL SELECTION

wide_known_nonsyn_phenotype3_keep <- wide_known_nonsyn_phenotype3_keep[,-1]

multi_lin_count1 <- glm(trans_MIC~., data = wide_known_nonsyn_phenotype3_keep)
summary(multi_lin_count1)

multi_lin_count2 <- step(multi_lin_count1, direction="both")
summary(multi_lin_count2) #AIC 

multi_lin_count3 <- step(multi_lin_count1, direction="backward")
summary(multi_lin_count3) #AIC 

multi_lin_count_null <- glm(trans_MIC~1, data = wide_known_nonsyn_phenotype5)
multi_lin_count4 <- step(multi_lin_count_null, direction="forward", scope=formula(multi_lin_count1), trace=0)
summary(multi_lin_count4) #AIC 

############################################################################################

## UNIVARIATE LINEAR REGRESSION - GENE COUNT - POPN STRUCTURE

uni_lin_count_phylo <- inner_join(complete_phylo, wide_known_nonsyn_phenotype3, join_by= "sample_id")
uni_lin_count_phylo1 <- uni_lin_count_phylo[,-1]
View(uni_lin_count_phylo1)
uni_lin_count_phylo1$Cluster.DSF <- pull(uni_lin_count_phylo1$Cluster.DSF)

pop_uni_lin_count_model <- list()
pop_uni_lin_count <- list()
for(i in 3:ncol(uni_lin_count_phylo1)) {
  pop_uni_lin_count_names <- names(uni_lin_count_phylo1)[[i]]
  pop_uni_lin_count_model[[pop_uni_lin_count_names]] <- glm(trans_MIC ~ uni_lin_count_phylo1[[i]]+ Cluster.DSF, data=uni_lin_count_phylo1)
  pop_uni_lin_count[[pop_uni_lin_count_names]] <- tidy(pop_uni_lin_count_model[[pop_uni_lin_count_names]])
}

pop_uni_lin_count

str(pop_uni_lin_count[[pop_uni_lin_count_names]])

length(pop_uni_lin_count) #[1] 16

pop_uni_lin_count_est <- matrix(NA, nrow=length(pop_uni_lin_count), ncol=2)

for(i in 1:length(pop_uni_lin_count)){
  pop_uni_lin_count_est[i,] <- c(pop_uni_lin_count[[i]]$estimate[2], 
                         pop_uni_lin_count[[i]]$p.value[2])
}

View(pop_uni_lin_count_est)

pop_uni_lin_count_p <- cbind(genes, pop_uni_lin_count_est)
names(pop_uni_lin_count_p) <- c('gene', 'Estimate', 'p-value')
View(pop_uni_lin_count_p)

##filter out significant ones (can use on dataframe too)
pop_lin_sig_count_est <- pop_uni_lin_count_p[pop_uni_lin_count_p[,3] <0.05,] ##how many significant variants
dim(pop_lin_sig_count_est) #[1] 5 3
View(pop_lin_sig_count_est)
pop_lin_sig_genescount <- as_tibble(pop_lin_sig_count_est[,1])
View(pop_lin_sig_genescount)

##subset significant ones
pop_lin_vectorgenec <- pull(pop_lin_sig_genescount, value)
uni_lin_count_phylo1_keep <- uni_lin_count_phylo1[,pop_lin_vectorgenec]
pop_samplecount_MIC <- select(uni_lin_count_phylo, sample_id, trans_MIC, Cluster.DSF)
uni_lin_count_phylo1_keep <- cbind(pop_samplecount_MIC, uni_lin_count_phylo1_keep)
View(uni_lin_count_phylo1_keep)
write.table(pop_lin_sig_count_est, file = "pop_lin_sig_count_est.txt", sep=",")

############################################################################################

##MULTIVARIATE LINEAR REGRESSION - GENE COUNT - POPN STRUCTURE

uni_lin_count_phylo1_keep1 <- uni_lin_count_phylo1_keep[,-1] 
uni_lin_count_phylo1_keep1$Cluster.DSF <- pull(uni_lin_count_phylo1_keep1$Cluster.DSF)
View(uni_lin_count_phylo1_keep1)

multi_lin_count_phylo1 <- glm(trans_MIC~., data = uni_lin_count_phylo1_keep1)
summary(multi_lin_count_phylo1)

multi_lin_count_phylo2 <- step(multi_lin_count_phylo1, direction="both")
summary(multi_lin_count_phylo2) 

#############################################################################################

## DATA WRANGLING - MIC BY VARIANT

variant_nonsyn$sample_id <- as.character(variant_nonsyn$sample_id) #data wrangling to prepare to join
phenotype_MIC$sample_id <- as.character(phenotype_MIC$sample_id)
variant_nonsyn_MIC <- left_join(phenotype_MIC, variant_nonsyn, by="sample_id")
View(variant_nonsyn_MIC) #547 variants, 621 samples

sapply(lapply(variant_nonsyn_MIC, unique), length) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these

uniquelengthlin <- sapply(variant_nonsyn_MIC,function(x) length(unique(x)))
variant_nonsyn_MIC <- subset(variant_nonsyn_MIC, select=uniquelengthlin>1)
View(variant_nonsyn_MIC) ##525 eligible variants

################################################################################################

## UNIVARIATE LINEAR REGRESSION - VARIANT

library(broom)

variant_nonsyn_MIC$trans_MIC <- as.numeric(variant_nonsyn_MIC$trans_MIC)

uni_lin_var <- list()
uni_lin_var_model <- list()
for(i in 3:ncol(variant_nonsyn_MIC)) {
  uni_lin_var_names <- names(variant_nonsyn_MIC)[[i]]
  uni_lin_var_model[[uni_lin_var_names]] <- glm(trans_MIC ~ variant_nonsyn_MIC[[i]], data=variant_nonsyn_MIC)
  uni_lin_var[[uni_lin_var_names]] <- tidy(uni_lin_var_model[[uni_lin_var_names]])
}

head(uni_lin_var)

str(uni_lin_var[[uni_lin_var_names]])

length(uni_lin_var) ##[1] 525

uni_lin_var_est <- matrix(NA, nrow=length(uni_lin_var), ncol=2) ##create empty matrix

for(i in 1:length(uni_lin_var)){
  uni_lin_var_est[i,] <- c(uni_lin_var[[i]]$estimate[2],
                           uni_lin_var[[i]]$p.value[2])
}

View(uni_lin_var_est)
variants2 <- as_tibble(colnames(variant_nonsyn_MIC))
variants2 <- variants2[-1:-2,]
View(variants2)

uni_lin_var_p <- cbind(variants2, uni_lin_var_est)
View(uni_lin_var_p)

##filter out significant ones (can use on dataframe too)
uni_lin_var_p[,3] < 0.05 ##returns true or false
lin_sig_var_est <-uni_lin_var_p[uni_lin_var_p[,3] <0.05,] ##103 significant variants
names(lin_sig_var_est) <- c('Variant', 'Estimate', 'p-value')
View(lin_sig_var_est)
write.table(lin_sig_var_est, file = "lin_sig_var_est.txt", sep=",")

lin_sig_variants <- as_tibble(lin_sig_var_est[,1])
View(lin_sig_variants)

##subset variant_nonsyn_MIC
vector2 <- pull(lin_sig_variants, value)
variant_nonsyn_MIC1 <- variant_nonsyn_MIC[,vector2, with=FALSE]
sample_MIC <- select(variant_nonsyn_MIC, sample_id, trans_MIC)
variant_nonsyn_MIC1 <- cbind(sample_MIC, variant_nonsyn_MIC1)
View(variant_nonsyn_MIC1)

#############################################################################################

## MULTIVARIATE LINEAR REGRESSION - VARIANT - MODEL SELECTION

variant_nonsyn_MIC1 <- variant_nonsyn_MIC1[,-1] 

multi_lin_variant1 <- glm(trans_MIC~., data = variant_nonsyn_MIC1)
summary(multi_lin_variant1)

multi_lin_variant2 <- step(multi_lin_variant1, direction="both")
summary(multi_lin_variant2) #AIC - R2

multi_lin_variant3 <- step(multi_lin_variant1, direction="backward")
summary(multi_lin_variant3) #AIC 

multi_lin_variant_null <- glm(trans_MIC~1, data = variant_nonsyn_MIC1)
multi_lin_variant4 <- step(multi_lin_variant_null, direction="forward", scope=formula(multi_lin_variant1), trace=0)
summary(multi_lin_variant4) #AIC 

##############################################################################################

## UNIVARIATE LINEAR REGRESSION - VARIANT - POPN STRUCTURE

uni_lin_var_phylo <- inner_join(complete_phylo, variant_nonsyn_MIC, join_by= "sample_id")
View(uni_lin_var_phylo)

clustervector2 <- uni_lin_var_phylo$Cluster.DSF

sapply(lapply(uni_lin_var_phylo, unique), length) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these
uniquelength3 <- sapply(uni_lin_var_phylo,function(x) length(unique(x)))
uni_lin_var_phylo <- subset(uni_lin_var_phylo, select=uniquelength3>1)
View(uni_lin_var_phylo) #519 variants. loses Cluster column so add this back in

uni_lin_var_phylo <- cbind(clustervector2, uni_lin_var_phylo)
uni_lin_var_phylo <- rename(uni_lin_var_phylo, Cluster.DSF = value)
uni_lin_var_phylo1 <- uni_lin_var_phylo[,-2] #remove sample ID column  
View(uni_lin_var_phylo1)

pop_uni_lin_var <- list()
pop_uni_lin_var_model <- list()
for(i in 3:ncol(uni_lin_var_phylo1)) {
  pop_uni_lin_var_names <- names(uni_lin_var_phylo1)[[i]]
  pop_uni_lin_var_model[[pop_uni_lin_var_names]] <- glm(trans_MIC ~ uni_lin_var_phylo1[[i]] + Cluster.DSF, data=uni_lin_var_phylo1)
  pop_uni_lin_var[[pop_uni_lin_var_names]] <- tidy(pop_uni_lin_var_model[[pop_uni_lin_var_names]])
}

str(pop_uni_lin_var[[pop_uni_lin_var_names]])

length(pop_uni_lin_var) #[1] 519

pop_uni_lin_var_est <- matrix(NA, nrow=length(pop_uni_lin_var), ncol=2)

for(i in 1:length(pop_uni_lin_var)){
  pop_uni_lin_var_est[i,] <- c(pop_uni_lin_var[[i]]$estimate[2], 
                               pop_uni_lin_var[[i]]$p.value[2])
}

variants3 <- as_tibble(colnames(uni_lin_var_phylo1))
variants3 <- variants3[-1:-2,]
View(variants3)

pop_uni_lin_var_p <- cbind(variants3, pop_uni_lin_var_est)
names(pop_uni_lin_var_p) <- c('Variant', 'Estimate', 'p-value')
View(pop_uni_lin_var_p)

##filter out significant variants
pop_uni_lin_var_p[,3] < 0.05 ##returns true or false
pop_lin_sig_var_est <- pop_uni_lin_var_p[pop_uni_lin_var_p[,3] <0.05,] ##52 significant variants
View(pop_lin_sig_var_est)
write.table(pop_lin_sig_var_est, file = "pop_lin_sig_var_est.txt", sep=",")

pop_lin_sig_variants <- as_tibble(pop_lin_sig_var_est[,1])

##subset significant ones
vector_pop_var_lin <- pull(pop_lin_sig_variants, value)
uni_lin_var_phylo_keep <- uni_lin_var_phylo[,vector_pop_var_lin]
View(uni_lin_var_phylo_keep)
cluster_MIC <- select(uni_lin_var_phylo, sample_id, trans_MIC, Cluster.DSF)
uni_lin_var_phylo_keep <- cbind(cluster_MIC, uni_lin_var_phylo_keep)
View(uni_lin_var_phylo_keep)

#############################################################################################

##MULTIVARIATE LINEAR REGRESSION - VARIANT - POPN STRUCTURE

uni_lin_var_phylo_keep1 <- uni_lin_var_phylo_keep[,-1] #remove sample_id
View(uni_lin_var_phylo_keep1)

multi_lin_var_phylo1 <- glm(trans_MIC~., data = uni_lin_var_phylo_keep1)
summary(multi_lin_var_phylo1)

multi_lin_var_phylo2 <- step(multi_lin_var_phylo1, direction="both")
summary(multi_lin_var_phylo2) 


########################################################################################################
