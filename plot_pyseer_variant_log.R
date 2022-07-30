#####################################################################################
### SCRIPT FOR GENERATING PLOT ######################################################
#####################################################################################

## This script is used to plot the effect sizes and allele frequencies of variants
# identified by GWAS


#####	READING INPUT FILES 

pyseer_output <- read.table("sig_dap_whole.txt", header=TRUE, sep="\t", dec=".")
View(pyseer_output)
snpeff_ann <- read.csv("efm_dap_dataset_v1.snippy.vcf.ann_table.csv")

lmPhen <- read.table("known_nonsyn_phenotype_res.txt", header=FALSE, sep="\t", dec = ".")

logpvalco <- -log10(as.numeric(1.29E-06))


num0 = length(which(lmPhen[,2]==0)) ##binary phenotype
num0 #422
num1 = length(which(lmPhen[,2]==1))
num1 #199

lmA <- snpeff_ann
lmT <- pyseer_output
View(lmT)

# Adding variant chromosome position

getposvariant = function(x){ pos = unlist(strsplit(x,"_"))[2]; return(pos); }

laPos = as.vector(sapply(as.vector(lmT[,"variant"]),getposvariant));

lmT[,"ps"] = as.numeric(laPos);

ldT = as.data.frame(lmT);
View(ldT)

################################################################################################################################

offsetlabel = 0.2;

sizeMutLabel = 2;

angleMutLabel = 0;

ybreaks = c(0,500000,1000000,1500000,2000000)

ylabels = c("0Mb","0.5Mb","1Mb","1.5Mb","2Mb")

######################################################################################################################################
###CREATING MANHATTAN PLOT USING GGPLOT	


print(paste("logpvalco used ", logpvalco,sep=""))

require(ggplot2)

p = ggplot(ldT, aes(x=ps, y=-log10(lrt.pvalue))) + 

geom_point(aes(size = af,colour = exp(beta))) +  

scale_colour_gradientn(colours=c("red","violet","blue")) +

geom_text(data = subset(ldT,-log10(lrt.pvalue)>=logpvalco), aes(x = ps, y = -log10(lrt.pvalue)+offsetlabel, label = variant), check_overlap = TRUE,size = sizeMutLabel,angle = angleMutLabel) +

labs(x = "Nucleotide position",y="-log10(p-val)",title = "Plot - variants identified through logistic LMM") +

theme(plot.title = element_text(family = "Calibri", face = "bold", size = (12)))

p


