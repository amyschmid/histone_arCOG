#arCOG functional enrichment -- Dulmage et al., "Growth phase specific modulation of cell morphology and gene expression by an archaeal histone protein"

#Load arCOG annotation file. 
load('~/NRC-1_arCOG_2012_final.Rdata')

ls() 

#Check to be sure an object called 'cogs' is listed this is your cog table

#Read in gene list file 

namelist = read.delim(file='~/namelist.txt', header=T, sep='\t')

#Be sure this list contains only the unique identifiers of your genes that matches those in the cog table (e.g. for H. salinarum, use VNG number identifiers)

#Calculate statistical information about the cogs represented in the list of genes using the function below.

# Cogtest: three inputs to this function: 
# 1) list of gene names from the cluster (namelist)
# 2) COG file
# 3) p-value cutoff.  

cogtest2=function(namelist,cogfile,pvalue) {
	cogs= subset(cogfile, is.element(cogfile$GeneName, namelist)==TRUE)
	clust= summary(cogs$funclass_name)
	resm <- matrix(0,31,3) 
	res <-data.frame(resm)
	rownames(res) <- names(clust)
	colnames(res) <- c("probability", "expect","count")
	all=summary(cogfile$funclass_name)
	totaln=sum(all)
	
		for (i in 1:length(clust)){
			
			res[i,1] <- phyper(clust[i], all[i], totaln-all[i], length(namelist),lower.tail=F)
			res[i,2] <- all[i]*(length(namelist)/totaln)
			res[i,3] <- clust[i]
				}
	subset(res,probability<=pvalue&count!=0)
		}
		
#Use the following function to look at the genes in your cluster associated with a particular COG.

#Note: the COGcategory name must be entered in quotes and add a space to the end of the category name of interest, e.g. 'transport '

cogset= function(namelist,cogfile, COGcategory) {
	subset(cogfile, is.element(cogfile$GeneName, namelist)&is.element(cogfile$funclass_name, COGcategory)==TRUE)
}

#Save results using write.table(), e.g. (x, file='clust296.txt',sep='\t')


