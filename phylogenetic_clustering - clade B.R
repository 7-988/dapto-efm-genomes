####################################################################################
##### PHYLOGENETIC CLUSTERING SCRIPT - CLADE B ISOLATES ############################
####################################################################################

## This script, adapted from https://github.com/francesccoll/phylogenetic_clustering,
## groups isolates into phylogenetic clusters.

require(phangorn,quietly = TRUE)
require(ape,quietly = TRUE)

tree_file = "efm_dap_dataset_v1.cladeB.rmRCB.rmMGE.iqtree.contree"
tree = read.tree(tree_file)
cat(paste("\t",length(tree$tip.label)," taxa found from phylogenetic tree\n",sep=""))
# 	30 taxa found from phylogenetic tree

# rooting on outgroup ERR6615664
tree = root(tree, which(tree$tip.label=="ERR6615664"))

laBS = as.numeric(tree$node.label);
tree2 = tree;
tree$node.label = paste("node",seq(1,tree$Nnode,1),sep="")

snp_file = "efm_dap_dataset_v1.cladeB.rmRCB.rmMGE.pairsnp.csv"
dismat = read.delim(snp_file,sep=",",header=F)
dim(dismat)
# [1] 30 31

# Replacing hashes with underscores in sample names
laColNames = as.vector(dismat[,1])
laColNames = gsub("#", "_", laColNames)
dismat = dismat[,-1]
colnames(dismat)= laColNames
rownames(dismat)= laColNames

# Making sure sample names in the tree match those in the distance matrix
identical(sort(tree$tip.label), sort(colnames(dismat)))
# [1] TRUE

laBS = as.numeric(tree2$node.label);


######################################################################################################################################
######							CLUSTERING ANALYSIS							                                                                           ######
######################################################################################################################################


######################################################################################################################################
######							LOADING REQUIRED LIBRARIES						                                                                      ######
######################################################################################################################################

require(ape, quietly=TRUE)     ## ape and gieger to read the tree
require(geiger, quietly=TRUE)
require(igraph, quietly=TRUE)  # igraph for the search
require(corpcor, quietly=TRUE)
require(TreeSim, quietly=TRUE)
require(phylobase, quietly=TRUE)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/BMhyd/BMhyd_1.2-8.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
require(BMhyd)		# loaded to use 'GetAncestor' function

######################################################################################################################################
######################################################################################################################################



######################################################################################################################################
######							     FUNCTIONS								       ######
######################################################################################################################################

## Three helper functions

## Given
##   a node, tree, and distance matrix
## Return
##   median pairwise patristic distance (MPPD) of its leaves

# Example of how to use it: 
# get.node.leaf.MPPD.FCC(1740,mytree,dismat)
# [1] 8

get.node.leaf.MPPD.FCC <- function(node,tree,distmat)
{
	#nlist <- node.leaves(tree,node); # --> ânode.leavesâ is being deprecated: use âtipsâ instead
	# NOTE that node refers to the internal node index position
	nlist <- tips(tree,node)
	foo <- distmat[nlist,nlist] 
	dis = median(foo[upper.tri(foo,diag=FALSE)])
	print(paste("node ",node, " median distance ",dis,sep=""))
	return(dis)
}


## Given
##   a node, tree, and distance matrix
## Return
##   maximum pairwise patristic distance (MPPD) of its leaves

get.node.leaf.MaxPPD.FCC <- function(node,tree,distmat)
{
	#nlist <- node.leaves(tree,node); # --> ânode.leavesâ is being deprecated: use âtipsâ instead
	# NOTE that node refers to the internal node index position
	nlist <- tips(tree,node)
	foo <- distmat[nlist,nlist] 
	dis = max(foo[upper.tri(foo,diag=FALSE)])
	print(paste("node ",node, " maximum distance ",dis,sep=""))
	return(dis)
}

# NOTE: functions calculating SNP distances with internal nodes are not used, i.e. get.node.full.MPPD

## Given
##   a tree and a distance matrix
## Return
##   a vector giving the median pairwise
##   patristic distance of the subtree under
##   each internal node
## SLOW!! May be a good idea to save/cache results

pdist.median.clusttree <- function(tree,distmat)
{
  ntips<- Ntip(tree)
  nint <- tree$Nnode ## number of internal nodes
  return(sapply((ntips+1):(ntips+nint),get.node.leaf.MPPD.FCC,tree,distmat))
}

## Given
##   a tree and a distance matrix
## Return
##   a vector giving the median pairwise
##   patristic distance of the subtree under
##   each internal node
## SLOW!! May be a good idea to save/cache results

pdist.max.clusttree <- function(tree,distmat)
{
	ntips<- Ntip(tree)
	nint <- tree$Nnode ## number of internal nodes
	laDis=sapply((ntips+1):(ntips+nint),get.node.leaf.MaxPPD.FCC,tree,distmat)
	return(laDis)
}


######################################################################################################################################
######################################################################################################################################

## NOTE: if bootstrap values are NA, these must be converted to non-NA (e.g. 0)
laBS[which(is.na(laBS))] = 0

laIntNodes = tree$node.label
distmat = dismat
distvec = NULL
thresh = 300; # maximum 20 SNPs within the clade
rthresh = 70; # bootstrap value cut-off
reliabilityvec = laBS

lmClusNode = mat.or.vec(1,4); # table to store node-defining nodes for clusters

  retval="tips";

  # NOTE: change to pdist.median.clusttree if the median pairwise patristic distance is chosen
  if(is.null(distvec))
  {
    cat("Calculating MPPD for each node ...\n")
    distvec <- pdist.max.clusttree(tree, distmat)
  }


  if(is.null(rthresh) || is.null(reliabilityvec)){
    reliabilityvec=NULL
    cat("No reliability thresholding.\n")
  }

  ## set up clustering
  ntips<-Ntip(tree)
  cnum <- 0 ## cluster number
  assign <- rep(0,ntips+tree$Nnode) ## cluster assignment
  igraph.tree <- graph.edgelist(tree$edge) ## tree in igraph form
  dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',order=TRUE,dist=TRUE)

  ## travese the tree in depth first order
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    print(i)
    ## skip leaves
    if(node < ntips+1)
    {
	print(paste("Node ",tree$tip.label[node]," skipt as it is a leaf",sep=""))
 	next 
    } else
    {
	print(paste("Node ",tree$node.label[node-ntips]," to be investigated",sep=""))
    }
    
   
    ## If the node's subtree is below the threshold, mark it and
    ## its subtree as members of a new cluster
    if(distvec[node-ntips]<=thresh && assign[node]<=0)
    {
	 ## If unreliable node --> investigate ancestral nodes
   	 if(!is.null(reliabilityvec) && reliabilityvec[node-ntips] < rthresh)
    	 { 
		print(paste("Node ",tree$node.label[node-ntips]," skipt as it has low BS value of ",reliabilityvec[node-ntips],sep=""))
		
		# Select most recent node ancestor with a BS value greater than the cut-off
		lbStop = FALSE;
		node2 = node;
		while(lbStop==FALSE)
		{
			node2 = GetAncestor(tree,as.numeric(node2));
			if(!is.null(reliabilityvec) && reliabilityvec[node2-ntips] < rthresh){ next; }
			else { lbStop = TRUE; }
		}

		print(paste("Node ",tree$node.label[node2-ntips]," to be used instead",sep=""));
		cnum <- cnum+1
      		subtree <- graph.dfs(igraph.tree,node2,neimode='out',unreachable=FALSE)$order
      		subtree <- subtree[!is.na(subtree)]
      		assign[subtree] <- cnum
      		newrow = c(cnum,as.character(tree$node.label[node2-ntips]),as.character(reliabilityvec[node2-ntips]),distvec[node2-ntips])
      		lmClusNode = rbind(lmClusNode,newrow)
		
    	 } else
	 {
      		cnum <- cnum+1
      		subtree <- graph.dfs(igraph.tree,node,neimode='out',unreachable=FALSE)$order
      		subtree <- subtree[!is.na(subtree)]
      		assign[subtree] <- cnum
      		newrow = c(cnum,as.character(tree$node.label[node-ntips]),as.character(reliabilityvec[node-ntips]),distvec[node-ntips])
      		lmClusNode = rbind(lmClusNode,newrow)
	 }
    } # end of if
  } # end of for loop

  ans <- list(membership=assign,allcsize=table(assign),leafclustsize=table(assign[1:ntips]),ntips=ntips,threshold=thresh)
  if(retval=="tips"){ans$membership <- ans$membership[1:ntips]}
  class(ans) <- c(class(ans),'p.cluster')
  clustering <- ans$membership

# Saving clustering results

laSam = tree$tip.label
laSamOrg = vector()

for(s in 1:length(laSam))
{
	aa = strsplit(laSam[s],"_")
	bb = aa[[1]]
	lsSamOrg = paste(bb[1],paste(bb[2],bb[3],sep="#"),sep="_")
	laSamOrg = c(laSamOrg,lsSamOrg)
}

lmTableB = mat.or.vec(length(laSam),2)
lmTableB[,1] = laSamOrg
lmTableB[,2] = clustering
ii = which(lmTableB[,2]==0)
lmTableB[ii,2]="singleton"
ii = which(lmTableB[,2]=="singleton")
lmTableB[ii,2]=paste("singleton",seq(1,length(ii),1),sep="")

colnames(lmTableB) =c("Taxa","Cluster.DSF")

output = paste("output.DSF_",thresh,"SNPsMax.BS",rthresh,"_clusters.txt",sep="")

write.table(lmTableB,file=output,sep="\t",col.names=T,row.names=F,quote=F)

colnames(lmClusNode) =c("Cluster.DSF","node","BT","Max_distance")
lmClusNode = lmClusNode[2:nrow(lmClusNode),]
output = paste("output.DSF_",thresh,"SNPsMax.BS",rthresh,"_clusters_nodes.txt",sep="")
write.table(lmClusNode,file=output,sep="\t",col.names=T,row.names=F,quote=F)

View(lmTableB)
