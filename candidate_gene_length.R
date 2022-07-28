########################################################################################
#### CANDIDATE GENE ANALYSIS ADJUSTING FOR GENE LENGTH #################################
########################################################################################

## This analysis repeats the 'gene count' analyses (univariate and multivariate, 
## with and without population structure, logistic and linear) adjusting for the 
## length of the gene.

########################################################################################

##Calculate gene length

bed_file <- read_csv("Known Genes.csv")
View(bed_file) ##bed file doesn't have 2872 on it

bed_file_name <- colnames(bed_file)
bed_file2 <- rbind(bed_file, bed_file_name)
colnames(bed_file2)<-c(1)
ncolsbed <- max(str_count(bed_file2[17,], "\t"))+1
ncolsbed #[1] 6
colmnbed <- paste(1:ncolsbed)
colmnbed <- as.character(colmnbed)
ncol(bed_file2) #[1] 1
bed_file2 <- separate(bed_file2,
                      col=1,
                      into=colmnbed,
                      sep="\t")
bed_file2 <- select(bed_file2, 2, 3, 4, 5)
colnames(bed_file2) <- c('Start', 'End', 'Gene_name', 'gene')
bed_file2 <- as_tibble(bed_file2)
bed_file2$Start <- as.numeric(bed_file2$Start)
bed_file2$End <- as.numeric(bed_file2$End)

bed_file2 <- mutate(bed_file2, length = End - Start)
View(bed_file2)

bed_file3 <- select(bed_file2, gene, length)
bed_file3 <- unique(bed_file3)
dim(bed_file3) #[1] 16 2
View(bed_file3)

#########################################################################################

#UNIVARIATE LOGISTIC - count, with length

View(known_nonsyn_phenotype2) #tibble for count
length_count <- full_join(known_nonsyn_phenotype2, bed_file3, by="gene")
View(length_count) #9936 entries (16x621)

length_count <- mutate(length_count, adj_count = number/length)
View(length_count)

length_count2 <- select(length_count, sample_id, daptomycin_RIS, gene, adj_count)

## Pivot_wider to get gene names across top

length_count2 <- as_tibble(length_count2)
wide_length_count <- length_count2 %>%
  pivot_wider(names_from = gene, values_from = adj_count) %>%
  group_by(sample_id)
dim(wide_length_count) #[1] 621 18
View(wide_length_count)

## Generate for loop to automate univariate logistic regressions, and present coefficients in tibble form 

library(broom)

uni_log_count_length <- list()
uni_log_count_length_model <- list()
for(i in 3:ncol(wide_length_count)) {
  uni_log_count_length_names <- names(wide_length_count)[[i]]
  uni_log_count_length_model[[uni_log_count_length_names]] <- glm(daptomycin_RIS ~ wide_length_count[[i]], family="binomial", data=wide_length_count)
  uni_log_count_length[[uni_log_count_length_names]] <- tidy(uni_log_count_length_model[[uni_log_count_length_names]])
}

str(uni_log_count_length[[uni_log_count_length_names]])

length(uni_log_count_length)

uni_log_count_length_est <- matrix(NA, nrow=length(uni_log_count_length), ncol=2)

for(i in 1:length(uni_log_count_length)){
  uni_log_count_length_est[i,] <- c(uni_log_count_length[[i]]$estimate[2], 
                         uni_log_count_length[[i]]$p.value[2])
}

View(uni_log_count_length_est)
genes
uni_log_count_length_p <- cbind(genes, uni_log_count_length_est)
names(uni_log_count_length_p) <- c('Gene', 'Estimate', 'p-value')
View(uni_log_count_length_p)

##filter out significant ones (can use on dataframe too)
log_sig_count_length_est <-uni_log_count_length_p[uni_log_count_length_p[,3] <0.05,]
dim(log_sig_count_length_est) #14 sig genes
View(log_sig_count_length_est)

write.table(log_sig_count_length_est, file = "length_log_sig_count_est.txt", sep=",")
log_sig_genescount_length <- as_tibble(log_sig_count_length_est[,1])
View(log_sig_genescount_length)

##subset significant ones
vectorgenec_l <- pull(log_sig_genescount_length, value)
wide_length_count_keep <- wide_length_count[,vectorgenec_l]
samplecount_RIS_length <- select(wide_length_count, sample_id, daptomycin_RIS)
wide_length_count_keep <- cbind(samplecount_RIS_length, wide_length_count_keep)
View(wide_length_count_keep)

################################################################################

##MULTIVARIATE LOGISTIC REGRESSION - GENE COUNT - MODEL SELECTION

wide_length_count_keep1 <- wide_length_count_keep[,-1]
length_multi_log_count1 <- glm(daptomycin_RIS~., data = wide_length_count_keep1, family="binomial")
summary(length_multi_log_count1)

pR2(length_multi_log_count1)['McFadden']

length_multi_log_count2 <- step(length_multi_log_count1, direction="both")
summary(length_multi_log_count2) #AIC 

pR2(length_multi_log_count2)['McFadden']

################################################################################################

##UNIVARIATE LOGISTIC REGRESSION - GENE COUNT - POPN STRUCTURE

length_uni_log_count_phylo <- inner_join(complete_phylo, wide_length_count, join_by= "sample_id")
View(length_uni_log_count_phylo)
length_uni_log_count_phylo1 <- length_uni_log_count_phylo[,-1]
View(length_uni_log_count_phylo1) #620 samples
length_uni_log_count_phylo1$Cluster.DSF <- pull(length_uni_log_count_phylo1$Cluster.DSF)

length_pop_uni_log_count_model <- list()
length_pop_uni_log_count <- list()
for(i in 3:ncol(length_uni_log_count_phylo1)) {
  length_pop_uni_log_count_names <- names(length_uni_log_count_phylo1)[[i]]
  length_pop_uni_log_count_model[[length_pop_uni_log_count_names]] <- glm(daptomycin_RIS ~ length_uni_log_count_phylo1[[i]] + Cluster.DSF, family="binomial", data=length_uni_log_count_phylo1)
  length_pop_uni_log_count[[length_pop_uni_log_count_names]] <- tidy(length_pop_uni_log_count_model[[length_pop_uni_log_count_names]])
}

head(length_pop_uni_log_count)

str(length_pop_uni_log_count[[length_pop_uni_log_count_names]])

length(length_pop_uni_log_count) #[1] 16

length_pop_uni_log_count_est <- matrix(NA, nrow=length(length_pop_uni_log_count), ncol=2)

for(i in 1:length(length_pop_uni_log_count)){
  length_pop_uni_log_count_est[i,] <- c(length_pop_uni_log_count[[i]]$estimate[2], 
                             length_pop_uni_log_count[[i]]$p.value[2])
}

View(length_pop_uni_log_count_est)

length_pop_uni_log_count_p <- cbind(genes, length_pop_uni_log_count_est)
names(length_pop_uni_log_count_p) <- c('gene', 'Estimate', 'p-value')
View(length_pop_uni_log_count_p)

##filter out significant ones (can use on dataframe too)
length_pop_log_sig_count_est <- length_pop_uni_log_count_p[length_pop_uni_log_count_p[,3] <0.05,] 
dim(length_pop_log_sig_count_est) #5 significant genes
View(length_pop_log_sig_count_est)
length_pop_log_sig_genescount <- as_tibble(length_pop_log_sig_count_est[,1])
View(length_pop_log_sig_genescount)

##subset significant ones
length_pop_vectorgenec <- pull(length_pop_log_sig_genescount, value)
length_uni_log_count_phylo_keep <- length_uni_log_count_phylo1[,length_pop_vectorgenec]
length_pop_samplecount_RIS <- select(length_uni_log_count_phylo, sample_id, daptomycin_RIS, Cluster.DSF)
length_uni_log_count_phylo_keep <- cbind(length_pop_samplecount_RIS, length_uni_log_count_phylo_keep)
View(length_uni_log_count_phylo_keep)

write.table(length_pop_log_sig_count_est, file = "length_pop_log_sig_count_est.txt", sep=",")

#################################################################################################

##MULTIVARIATE LOGISTIC REGRESSION - GENE COUNT - POPN STRUCTURE

length_uni_log_count_phylo_keep1 <- length_uni_log_count_phylo_keep[,-1] 
length_uni_log_count_phylo_keep1$Cluster.DSF <- pull(length_uni_log_count_phylo_keep1$Cluster.DSF)
View(length_uni_log_count_phylo_keep1)

length_multi_log_count_phylo1 <- glm(daptomycin_RIS~., data = length_uni_log_count_phylo_keep1, family="binomial")
summary(length_multi_log_count_phylo1)

pR2(length_multi_log_count_phylo1)['McFadden']

length_multi_log_count_phylo2 <- step(length_multi_log_count_phylo1, direction="both")
summary(length_multi_log_count_phylo2) 

pR2(length_multi_log_count_phylo2)['McFadden']

#######################################################################################################

################################################################################################################

## UNIVARIATE LINEAR REGRESSION (MIC) - COUNT

View(known_table_nonsyn3)

known_table_nonsyn3$sample_id <- as.character(known_table_nonsyn3$sample_id) #data wrangling to prepare to join
phenotype_MIC$sample_id <- as.character(phenotype_MIC$sample_id)
known_nonsyn_phenotype3 <- left_join(phenotype_MIC, known_table_nonsyn3, by="sample_id")
View(known_nonsyn_phenotype3)

lin_length_count <- full_join(known_nonsyn_phenotype3, bed_file3, by="gene")
View(lin_length_count)

lin_length_count <- mutate(lin_length_count, adj_count = number/length)
View(lin_length_count)

lin_length_count2 <- select(lin_length_count, sample_id, trans_MIC, gene, adj_count)

View(lin_length_count2)

## Pivot_wider to get gene names across top

lin_length_count2 <- distinct(lin_length_count2)
lin_length_count2 <- as_tibble(lin_length_count2)
View(lin_length_count2)
wide_lin_length_count <- lin_length_count2 %>%
  pivot_wider(names_from = gene, values_from = adj_count) %>%
  group_by(sample_id)
View(wide_lin_length_count)

## Generate for loop to automate univariate linear regressions, and present coefficients in tibble form 

library(broom)

length_uni_lin_count_model <- list()
length_uni_lin_count <- list()
for(i in 3:ncol(wide_lin_length_count)) {
  length_uni_lin_count_names <- names(wide_lin_length_count)[[i]]
  length_uni_lin_count_model[[length_uni_lin_count_names]] <- glm(trans_MIC ~ wide_lin_length_count[[i]], data=wide_lin_length_count)
  length_uni_lin_count[[length_uni_lin_count_names]] <- tidy(length_uni_lin_count_model[[length_uni_lin_count_names]])
}

length_uni_lin_count

str(length_uni_lin_count[[length_uni_lin_count_names]])

length(length_uni_lin_count)

length_uni_lin_count_est <- matrix(NA, nrow=length(length_uni_lin_count), ncol=2)

for(i in 1:length(length_uni_lin_count)){
  length_uni_lin_count_est[i,] <- c(length_uni_lin_count[[i]]$estimate[2], 
                         length_uni_lin_count[[i]]$p.value[2])
}

View(length_uni_lin_count_est)

length_uni_lin_count_p <- cbind(genes, length_uni_lin_count_est)
names(length_uni_lin_count_p) <- c('gene', 'Estimate', 'p-value')
View(length_uni_lin_count_p)

##filter out significant ones (can use on dataframe too)
length_lin_sig_count_est <- length_uni_lin_count_p[length_uni_lin_count_p[,3] <0.05,] 
dim(length_lin_sig_count_est) #[1] 9 3
View(length_lin_sig_count_est)
length_lin_sig_genescount <- as_tibble(length_lin_sig_count_est[,1])
View(length_lin_sig_genescount)

##subset significant ones
length_lin_vectorgenec <- pull(length_lin_sig_genescount, value)
wide_lin_length_count_keep <- wide_lin_length_count[,length_lin_vectorgenec]
length_samplecount_MIC <- select(wide_lin_length_count, sample_id, trans_MIC)
wide_lin_length_count_keep <- cbind(length_samplecount_MIC, wide_lin_length_count_keep)

write.table(length_lin_sig_count_est, file = "length_lin_sig_count_est.txt", sep=",")

###########################################################################################

## MULTIVARIATE LINEAR REGRESSION - GENE COUNT - MODEL SELECTION

wide_lin_length_count_keep <- wide_lin_length_count_keep[,-1]

length_multi_lin_count1 <- glm(trans_MIC~., data = wide_lin_length_count_keep)
summary(length_multi_lin_count1)

length_multi_lin_count2 <- step(length_multi_lin_count1, direction="both")
summary(length_multi_lin_count2) #AIC 

############################################################################################

## UNIVARIATE LINEAR REGRESSION - GENE COUNT - POPN STRUCTURE

length_uni_lin_count_phylo <- inner_join(complete_phylo, wide_lin_length_count, join_by= "sample_id")
length_uni_lin_count_phylo1 <- length_uni_lin_count_phylo[,-1]
View(length_uni_lin_count_phylo1)
length_uni_lin_count_phylo1$Cluster.DSF <- pull(length_uni_lin_count_phylo1$Cluster.DSF)

length_pop_uni_lin_count_model <- list()
length_pop_uni_lin_count <- list()
for(i in 3:ncol(length_uni_lin_count_phylo1)) {
  length_pop_uni_lin_count_names <- names(length_uni_lin_count_phylo1)[[i]]
  length_pop_uni_lin_count_model[[length_pop_uni_lin_count_names]] <- glm(trans_MIC ~ length_uni_lin_count_phylo1[[i]]+ Cluster.DSF, data=length_uni_lin_count_phylo1)
  length_pop_uni_lin_count[[length_pop_uni_lin_count_names]] <- tidy(length_pop_uni_lin_count_model[[length_pop_uni_lin_count_names]])
}

length_pop_uni_lin_count

str(length_pop_uni_lin_count[[length_pop_uni_lin_count_names]])

length(length_pop_uni_lin_count)

length_pop_uni_lin_count_est <- matrix(NA, nrow=length(length_pop_uni_lin_count), ncol=2)

for(i in 1:length(length_pop_uni_lin_count)){
  length_pop_uni_lin_count_est[i,] <- c(length_pop_uni_lin_count[[i]]$estimate[2], 
                             length_pop_uni_lin_count[[i]]$p.value[2])
}

View(length_pop_uni_lin_count_est)

length_pop_uni_lin_count_p <- cbind(genes, length_pop_uni_lin_count_est)
names(length_pop_uni_lin_count_p) <- c('gene', 'Estimate', 'p-value')
View(length_pop_uni_lin_count_p)

##filter out significant ones (can use on dataframe too)
length_pop_lin_sig_count_est <- length_pop_uni_lin_count_p[length_pop_uni_lin_count_p[,3] <0.05,] 
dim(length_pop_lin_sig_count_est) #[1] 5 3
View(length_pop_lin_sig_count_est)
length_pop_lin_sig_genescount <- as_tibble(length_pop_lin_sig_count_est[,1])
View(length_pop_lin_sig_genescount)

##subset significant ones
length_pop_lin_vectorgenec <- pull(length_pop_lin_sig_genescount, value)
length_uni_lin_count_phylo1_keep <- length_uni_lin_count_phylo1[,length_pop_lin_vectorgenec]
length_pop_samplecount_MIC <- select(length_uni_lin_count_phylo, sample_id, trans_MIC, Cluster.DSF)
length_uni_lin_count_phylo1_keep <- cbind(length_pop_samplecount_MIC, length_uni_lin_count_phylo1_keep)

write.table(length_pop_lin_sig_count_est, file = "length_pop_lin_sig_count_est.txt", sep=",")

############################################################################################

##MULTIVARIATE LINEAR REGRESSION - GENE COUNT - POPN STRUCTURE

length_uni_lin_count_phylo1_keep1 <- length_uni_lin_count_phylo1_keep[,-1] 
length_uni_lin_count_phylo1_keep1$Cluster.DSF <- pull(length_uni_lin_count_phylo1_keep1$Cluster.DSF)
View(length_uni_lin_count_phylo1_keep1)

length_multi_lin_count_phylo1 <- lm(trans_MIC~., data = length_uni_lin_count_phylo1_keep1)
summary(length_multi_lin_count_phylo1)

length_multi_lin_count_phylo2 <- step(length_multi_lin_count_phylo1, direction="both")
summary(length_multi_lin_count_phylo2) 
