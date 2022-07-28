#######################################################################################
##### VARIANT CANDIDATE GENE ANALYSIS USING ONLY VARIANTS IN INDIVIDUALLY SIG GENES ###
#######################################################################################

## This analysis repeats the variant candidate gene analysis using only variants 
## in genes found to be individually significant from the following 
## LOGISTIC --> 
  ##univariate logistic population structure - gene binary
  ##univariate logistic population structure length - gene count
  ##combine these --> logistic variant analyses x4 (U, UP, M, MP)
##LINEAR -->
  ##univariate linear population structure - gene binary
  ##univariate linear population structure length - gene count
  ##combine these --> linear variant analyses x4 (U, UP, M, MP)

########################################################################################

## LOGISTIC 

##combined genes significant from gene binary (using popn structure) and gene count (using popn structure and length adjustment)
##remove duplicates
##use this to select the variants to run/use

# run 4 analyses

sig_genes <- distinct(rbind(pop_log_sig_genesbin, length_pop_log_sig_genescount)) #8 genes 
View(sig_genes)

# subset 

View(known_table_nonsyn)
dim(known_table_nonsyn)
sig_genes_vec <- pull(sig_genes, value)
View(sig_genes_vec)

log_selected_variants <- transp_full_variant[is.element(transp_full_variant$gene, sig_genes_vec),]
View(log_selected_variants)
dim(log_selected_variants) #[1] 940 676

##need to add to phenotype data and remove certain samples etc.

log_selected_nonsyn <- filter(log_selected_variants, severity!="LOW") #02872 variant already excluded by gene
View(log_selected_nonsyn) #315 variants
log_selected_nonsyn <- log_selected_nonsyn[,-(2:5)]

t_log_selected_nonsyn <- t(log_selected_nonsyn)
View(t_log_selected_nonsyn)

colnames(t_log_selected_nonsyn) <- t_log_selected_nonsyn[1,]
t_log_selected_nonsyn <- t_log_selected_nonsyn %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(sample_id = rowname)
View(t_log_selected_nonsyn)
t_log_selected_nonsyn <- t_log_selected_nonsyn[-1,]
t_log_selected_nonsyn <- t_log_selected_nonsyn[-671,]
View(t_log_selected_nonsyn) #670x547 variants plus sample_id column

##recode 2 to 1 (proper binary)

t_log_selected_nonsyn[t_log_selected_nonsyn == "2"] <- "1"
View(t_log_selected_nonsyn)

## add phenotype_binary
t_log_selected_nonsyn$sample_id <- as.character(t_log_selected_nonsyn$sample_id) #data wrangling to prepare to join
phenotype_binary$sample_id <- as.character(phenotype_binary$sample_id)
sel_variant_nonsyn_pheno <- left_join(phenotype_binary, t_log_selected_nonsyn, by="sample_id")
View(sel_variant_nonsyn_pheno)
dim(sel_variant_nonsyn_pheno) #[1] 621 317
sum(is.na(sel_variant_nonsyn_pheno)) #[1] 0

##Recode daptomycin_RIS to be binary (1=NS, 0=S)

sel_variant_nonsyn_pheno$daptomycin_RIS <- ifelse(sel_variant_nonsyn_pheno$daptomycin_RIS == "NS", 1, 0)
View(sel_variant_nonsyn_pheno)

head(sapply(lapply(sel_variant_nonsyn_pheno, unique), length)) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these

sel_uniquelength <- sapply(sel_variant_nonsyn_pheno,function(x) length(unique(x)))
sel_variant_nonsyn_pheno <- subset(sel_variant_nonsyn_pheno, select=sel_uniquelength>1)
dim(sel_variant_nonsyn_pheno) #[1] 621 310

##UNIVARIATE LOGISTIC REGRESSION - VARIANT/BINARY

library(broom)

View(sel_variant_nonsyn_pheno)

sel_uni_log_var <- list()
sel_uni_log_var_model <- list()
for(i in 3:ncol(sel_variant_nonsyn_pheno)) {
  sel_uni_log_var_names <- names(sel_variant_nonsyn_pheno)[[i]]
  sel_uni_log_var_model[[sel_uni_log_var_names]] <- glm(daptomycin_RIS ~ sel_variant_nonsyn_pheno[[i]], family="binomial", data=sel_variant_nonsyn_pheno)
  sel_uni_log_var[[sel_uni_log_var_names]] <- tidy(sel_uni_log_var_model[[sel_uni_log_var_names]])
}

str(sel_uni_log_var[[sel_uni_log_var_names]])

length(sel_uni_log_var) #[1] 308

sel_uni_log_var_est <- matrix(NA, nrow=length(sel_uni_log_var), ncol=2)

for(i in 1:length(sel_uni_log_var)){
  sel_uni_log_var_est[i,] <- c(sel_uni_log_var[[i]]$estimate[2], 
                           sel_uni_log_var[[i]]$p.value[2])
}

dim(sel_uni_log_var_est) #[1] 308 2
sel_variants <- as_tibble(colnames(sel_variant_nonsyn_pheno))
sel_variants <- sel_variants[-1:-2,]
View(sel_variants)

sel_uni_log_var_p <- cbind(sel_variants, sel_uni_log_var_est)
names(sel_uni_log_var_p) <- c('Variant', 'Estimate', 'p-value')
View(sel_uni_log_var_p)

##filter out significant ones (can use on dataframe too)
sel_log_sig_var_est <-sel_uni_log_var_p[sel_uni_log_var_p[,3] <0.05,] ##50 significant variants
View(sel_log_sig_var_est)
sel_log_sig_variants <- as_tibble(sel_log_sig_var_est[,1])
View(sel_log_sig_variants)

##subset variant_nonsyn_pheno
sel_vector <- pull(sel_log_sig_variants, value)
sel_variant_nonsyn_pheno1 <- sel_variant_nonsyn_pheno[,sel_vector, with=FALSE]
View(sel_variant_nonsyn_pheno1)
sel_sample_RIS <- select(sel_variant_nonsyn_pheno, sample_id, daptomycin_RIS)
sel_variant_nonsyn_pheno1 <- cbind(sel_sample_RIS, sel_variant_nonsyn_pheno1)
View(sel_variant_nonsyn_pheno1)

write.table(sel_log_sig_var_est, file = "sel_log_sig_var_est.txt", sep=",")

#######################################################################################

##MULTIVARIATE LOGISTIC REGRESSION  - VARIANT - MODEL SELECTION

sel_variant_nonsyn_pheno1 <- sel_variant_nonsyn_pheno1[,-1] 
View(sel_variant_nonsyn_pheno1)
sel_multi_log_variant1 <- glm(daptomycin_RIS~., data = sel_variant_nonsyn_pheno1, family="binomial") ##some predictors have perfect correlation/linear relationship?
summary(sel_multi_log_variant1)

pR2(sel_multi_log_variant1)['McFadden']

sel_multi_log_variant2 <- step(sel_multi_log_variant1, direction="both")
summary(sel_multi_log_variant2) #AIC 

pR2(sel_multi_log_variant2)['McFadden']


############################################################################################

##UNIVARIATE LOGISTIC REGRESSION - VARIANT - POPN STRUCTURE

sel_uni_log_var_phylo <- inner_join(complete_phylo, sel_variant_nonsyn_pheno, join_by= "sample_id")
View(sel_uni_log_var_phylo) #620x308 var

sel_clustervector <- sel_uni_log_var_phylo$Cluster.DSF

sapply(lapply(sel_uni_log_var_phylo, unique), length) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these
sel_uniquelength2 <- sapply(sel_uni_log_var_phylo,function(x) length(unique(x)))
sel_uni_log_var_phylo <- subset(sel_uni_log_var_phylo, select=sel_uniquelength2>1)
dim(sel_uni_log_var_phylo) #[1] 620 307

sel_uni_log_var_phylo <- cbind(sel_clustervector, sel_uni_log_var_phylo)
sel_uni_log_var_phylo <- rename(sel_uni_log_var_phylo, Cluster.DSF = value)
sel_uni_log_var_phylo1 <- sel_uni_log_var_phylo[,-2]  
View(sel_uni_log_var_phylo1)

sel_pop_uni_log_var <- list()
sel_pop_uni_log_var_model <- list()
for(i in 3:ncol(sel_uni_log_var_phylo1)) {
  sel_pop_uni_log_var_names <- names(sel_uni_log_var_phylo1)[[i]]
  sel_pop_uni_log_var_model[[sel_pop_uni_log_var_names]] <- glm(daptomycin_RIS ~ sel_uni_log_var_phylo1[[i]]+ Cluster.DSF, family="binomial", data=sel_uni_log_var_phylo1)
  sel_pop_uni_log_var[[sel_pop_uni_log_var_names]] <- tidy(sel_pop_uni_log_var_model[[sel_pop_uni_log_var_names]])
}

head(sel_pop_uni_log_var)

str(sel_pop_uni_log_var[[sel_pop_uni_log_var_names]])

length(sel_pop_uni_log_var)#[1] 305

sel_pop_uni_log_var_est <- matrix(NA, nrow=length(sel_pop_uni_log_var), ncol=2)

for(i in 1:length(sel_pop_uni_log_var)){
  sel_pop_uni_log_var_est[i,] <- c(sel_pop_uni_log_var[[i]]$estimate[2], 
                               sel_pop_uni_log_var[[i]]$p.value[2])
}

View(sel_pop_uni_log_var_est)
sel_variants2 <- as_tibble(colnames(sel_uni_log_var_phylo1))
sel_variants2 <- sel_variants2[-1:-2,]
View(sel_variants2)

sel_pop_uni_log_var_p <- cbind(sel_variants2, sel_pop_uni_log_var_est)
names(sel_pop_uni_log_var_p) <- c('Variant', 'Estimate', 'p-value')
View(sel_pop_uni_log_var_p)

##filter out significant variants
sel_pop_log_sig_var_est <- sel_pop_uni_log_var_p[sel_pop_uni_log_var_p[,3] <0.05,] ##30 significant variants
View(sel_pop_log_sig_var_est)
write.table(sel_pop_log_sig_var_est, file = "sel_pop_log_sig_var_est.txt", sep=",")

sel_pop_log_sig_variants <- as_tibble(sel_pop_log_sig_var_est[,1])
View(sel_pop_log_sig_variants)

##subset significant ones
sel_vector_pop_var <- pull(sel_pop_log_sig_variants, value)
sel_uni_log_var_phylo_keep <- sel_uni_log_var_phylo[,sel_vector_pop_var]
View(sel_uni_log_var_phylo_keep)
sel_cluster_RIS <- select(sel_uni_log_var_phylo, sample_id, daptomycin_RIS, Cluster.DSF)
sel_uni_log_var_phylo_keep <- cbind(sel_cluster_RIS, sel_uni_log_var_phylo_keep)
View(sel_uni_log_var_phylo_keep)


#############################################################################################

##MULTIVARIATE LOGISTIC REGRESSION - VARIANT - POPN STRUCTURE

sel_uni_log_var_phylo_keep1 <- sel_uni_log_var_phylo_keep[,-1] #remove sample_id
View(sel_uni_log_var_phylo_keep1)

sel_multi_log_var_phylo1 <- glm(daptomycin_RIS~., data = sel_uni_log_var_phylo_keep1, family="binomial") ##some exact fitting of probabilities?
summary(sel_multi_log_var_phylo1)

pR2(sel_multi_log_var_phylo1)['McFadden']

sel_multi_log_var_phylo2 <- step(sel_multi_log_var_phylo1, direction="both")
summary(sel_multi_log_var_phylo2) 

pR2(sel_multi_log_var_phylo2)['McFadden']

###########################################################################################
###########################################################################################


## LINEAR 

##combined genes significant from gene binary and gene count
##remove duplicates
##use this to select the variants to run/use

# run 4 analyses. NEED TO EDIT THIS ONCE RUN !!

lin_sig_genes <- distinct(rbind(pop_lin_sig_genesbin, length_pop_lin_sig_genescount)) #9 genes 
View(lin_sig_genes)

lin_sig_genes_vec <- pull(lin_sig_genes, value)
View(lin_sig_genes_vec)

lin_selected_variants <- transp_full_variant[is.element(transp_full_variant$gene, lin_sig_genes_vec),]
View(lin_selected_variants)
dim(lin_selected_variants) #[1] 894 676

##need to add to phenotype data and remove certain samples etc.

lin_selected_nonsyn <- filter(lin_selected_variants, severity!="LOW") #02872 variant already excluded by gene
View(lin_selected_nonsyn) #300 variants
lin_selected_nonsyn <- lin_selected_nonsyn[,-(2:5)]

t_lin_selected_nonsyn <- t(lin_selected_nonsyn)
View(t_lin_selected_nonsyn)

colnames(t_lin_selected_nonsyn) <- t_lin_selected_nonsyn[1,]
t_lin_selected_nonsyn <- t_lin_selected_nonsyn %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(sample_id = rowname)
View(t_lin_selected_nonsyn)
t_lin_selected_nonsyn <- t_lin_selected_nonsyn[-1,]
t_lin_selected_nonsyn <- t_lin_selected_nonsyn[-671,]
View(t_lin_selected_nonsyn) #670x300 variants plus sample_id column

##recode 2 to 1 (proper binary)

t_lin_selected_nonsyn[t_lin_selected_nonsyn == "2"] <- "1"
View(t_lin_selected_nonsyn)

## add phenotype_binary

t_lin_selected_nonsyn$sample_id <- as.character(t_lin_selected_nonsyn$sample_id) #data wrangling to prepare to join
phenotype_MIC$sample_id <- as.character(phenotype_MIC$sample_id)
lin_sel_variant_MIC <- left_join(phenotype_MIC, t_lin_selected_nonsyn, by="sample_id")
View(lin_sel_variant_MIC)
dim(lin_sel_variant_MIC) #[1] 621 302
sum(is.na(lin_sel_variant_MIC)) #[1] 0

head(sapply(lapply(lin_sel_variant_MIC, unique), length)) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these

lin_sel_uniquelength <- sapply(lin_sel_variant_MIC,function(x) length(unique(x)))
lin_sel_variant_MIC <- subset(lin_sel_variant_MIC, select=lin_sel_uniquelength>1)
dim(lin_sel_variant_MIC) #[1] 621 294

## UNIVARIATE LINEAR REGRESSION - VARIANT

library(broom)

lin_sel_variant_MIC$trans_MIC <- as.numeric(lin_sel_variant_MIC$trans_MIC)

sel_uni_lin_var <- list()
sel_uni_lin_var_model <- list()
for(i in 3:ncol(lin_sel_variant_MIC)) {
  sel_uni_lin_var_names <- names(lin_sel_variant_MIC)[[i]]
  sel_uni_lin_var_model[[sel_uni_lin_var_names]] <- glm(trans_MIC ~ lin_sel_variant_MIC[[i]], data=lin_sel_variant_MIC)
  sel_uni_lin_var[[sel_uni_lin_var_names]] <- tidy(sel_uni_lin_var_model[[sel_uni_lin_var_names]])
}

head(sel_uni_lin_var)

str(sel_uni_lin_var[[sel_uni_lin_var_names]])

length(sel_uni_lin_var) ##[1] 292

sel_uni_lin_var_est <- matrix(NA, nrow=length(sel_uni_lin_var), ncol=2) ##create empty matrix

for(i in 1:length(sel_uni_lin_var)){
  sel_uni_lin_var_est[i,] <- c(sel_uni_lin_var[[i]]$estimate[2],
                           sel_uni_lin_var[[i]]$p.value[2])
}

View(sel_uni_lin_var_est)
sel_lin_variants2 <- as_tibble(colnames(lin_sel_variant_MIC))
sel_lin_variants2 <- sel_lin_variants2[-1:-2,]
View(sel_lin_variants2)

sel_uni_lin_var_p <- cbind(sel_lin_variants2, sel_uni_lin_var_est)
View(sel_uni_lin_var_p)

##filter out significant ones (can use on dataframe too)
sel_lin_sig_var_est <- sel_uni_lin_var_p[sel_uni_lin_var_p[,3] <0.05,] ##80 significant variants
names(sel_lin_sig_var_est) <- c('Variant', 'Estimate', 'p-value')
View(sel_lin_sig_var_est)
write.table(sel_lin_sig_var_est, file = "sel_lin_sig_var_est.txt", sep=",")

sel_lin_sig_variants <- as_tibble(sel_lin_sig_var_est[,1])
View(sel_lin_sig_variants)

##subset variant_nonsyn_MIC
sel_vector2 <- pull(sel_lin_sig_variants, value)
lin_sel_variant_MIC1 <- lin_sel_variant_MIC[,sel_vector2, with=FALSE]
sel_sample_MIC <- select(lin_sel_variant_MIC, sample_id, trans_MIC)
lin_sel_variant_MIC1 <- cbind(sel_sample_MIC, lin_sel_variant_MIC1)
View(lin_sel_variant_MIC1)

#############################################################################################

## MULTIVARIATE LINEAR REGRESSION - VARIANT - MODEL SELECTION

lin_sel_variant_MIC1 <- lin_sel_variant_MIC1[,-1] 

lin_multi_lin_variant1 <- glm(trans_MIC~., data = lin_sel_variant_MIC1)
summary(lin_multi_lin_variant1)

lin_multi_lin_variant2 <- step(lin_multi_lin_variant1, direction="both")
summary(lin_multi_lin_variant2) #AIC - R2
 

##############################################################################################

## UNIVARIATE LINEAR REGRESSION - VARIANT - POPN STRUCTURE

sel_uni_lin_var_phylo <- inner_join(complete_phylo, lin_sel_variant_MIC, join_by= "sample_id")
View(sel_uni_lin_var_phylo)

sel_clustervector2 <- sel_uni_lin_var_phylo$Cluster.DSF

sapply(lapply(sel_uni_lin_var_phylo, unique), length) ##indicates that after removal of certain samples, some variant columns no longer have any variation. need to remove these
sel_uniquelength3 <- sapply(sel_uni_lin_var_phylo,function(x) length(unique(x)))
sel_uni_lin_var_phylo <- subset(sel_uni_lin_var_phylo, select=sel_uniquelength3>1)
View(sel_uni_lin_var_phylo) #291 variants. loses Cluster column so add this back in

sel_uni_lin_var_phylo <- cbind(sel_clustervector2, sel_uni_lin_var_phylo)
sel_uni_lin_var_phylo <- rename(sel_uni_lin_var_phylo, Cluster.DSF = value)
sel_uni_lin_var_phylo1 <- sel_uni_lin_var_phylo[,-2] #remove sample ID column  
View(sel_uni_lin_var_phylo1)

sel_pop_uni_lin_var <- list()
sel_pop_uni_lin_var_model <- list()
for(i in 3:ncol(sel_uni_lin_var_phylo1)) {
  sel_pop_uni_lin_var_names <- names(sel_uni_lin_var_phylo1)[[i]]
  sel_pop_uni_lin_var_model[[sel_pop_uni_lin_var_names]] <- glm(trans_MIC ~ sel_uni_lin_var_phylo1[[i]] + Cluster.DSF, data=sel_uni_lin_var_phylo1)
  sel_pop_uni_lin_var[[sel_pop_uni_lin_var_names]] <- tidy(sel_pop_uni_lin_var_model[[sel_pop_uni_lin_var_names]])
}

str(sel_pop_uni_lin_var[[sel_pop_uni_lin_var_names]])

length(sel_pop_uni_lin_var) #[1] 291

sel_pop_uni_lin_var_est <- matrix(NA, nrow=length(sel_pop_uni_lin_var), ncol=2)

for(i in 1:length(sel_pop_uni_lin_var)){
  sel_pop_uni_lin_var_est[i,] <- c(sel_pop_uni_lin_var[[i]]$estimate[2], 
                               sel_pop_uni_lin_var[[i]]$p.value[2])
}

sel_variants3 <- as_tibble(colnames(sel_uni_lin_var_phylo1))
sel_variants3 <- sel_variants3[-1:-2,]
View(sel_variants3)

sel_pop_uni_lin_var_p <- cbind(sel_variants3, sel_pop_uni_lin_var_est)
names(sel_pop_uni_lin_var_p) <- c('Variant', 'Estimate', 'p-value')
View(sel_pop_uni_lin_var_p)

##filter out significant variants
sel_pop_lin_sig_var_est <- sel_pop_uni_lin_var_p[sel_pop_uni_lin_var_p[,3] <0.05,] ##38 significant variants
View(sel_pop_lin_sig_var_est)
write.table(sel_pop_lin_sig_var_est, file = "sel_pop_lin_sig_var_est.txt", sep=",")

sel_pop_lin_sig_variants <- as_tibble(sel_pop_lin_sig_var_est[,1])

##subset significant ones
sel_vector_pop_var_lin <- pull(sel_pop_lin_sig_variants, value)
sel_uni_lin_var_phylo_keep <- sel_uni_lin_var_phylo[,sel_vector_pop_var_lin]
View(sel_uni_lin_var_phylo_keep)
sel_cluster_MIC <- select(sel_uni_lin_var_phylo, sample_id, trans_MIC, Cluster.DSF)
sel_uni_lin_var_phylo_keep <- cbind(cluster_MIC, sel_uni_lin_var_phylo_keep)
View(sel_uni_lin_var_phylo_keep)

#############################################################################################

##MULTIVARIATE LINEAR REGRESSION - VARIANT - POPN STRUCTURE

sel_uni_lin_var_phylo_keep1 <- sel_uni_lin_var_phylo_keep[,-1] #remove sample_id
View(sel_uni_lin_var_phylo_keep1)

sel_multi_lin_var_phylo1 <- glm(trans_MIC~., data = sel_uni_lin_var_phylo_keep1)
summary(sel_multi_lin_var_phylo1)

sel_multi_lin_var_phylo2 <- step(sel_multi_lin_var_phylo1, direction="both")
summary(sel_multi_lin_var_phylo2) 


########################################################################################################



