## -----------------------------------------------------------------
## -----------------------------------------------------------------
## Define important variables and run all analyses
## -----------------------------------------------------------------
## -----------------------------------------------------------------

options(width=60)

## Install needed packages
source("installNeededPackages.R")

##------------------------
##Files and other settings
##------------------------

## Settings for PCA.Rnw

## Important configuration:
output.dir <- "../output_allsamples"

##Input / output data files:
expression.infile <- "../input/gene_expression.tab"
metadata.infile <- "../input/gene expression metadata gold.txt"
scores.outfile <- paste(output.dir, "scores_PhenotypeISCORE_0.4.csv", sep="/")
clustering.outfile <- paste(output.dir, "/maaslin_clustering.pcl", sep="")
expr.all.outfile <- paste(output.dir, "expr_all.csv", sep="/")
bugs.all.outfile <- paste(output.dir, "bugs_all.csv", sep="/")

if( !file.exists( output.dir ) )
    dir.create( output.dir )

## runKnitr reads a file with extenstion .Rnw (e.g. myfile.Rnw), which
## should use argslist for any needed input variables, and create a PDF
## output.dir/myfile.pdf.
runKnitr <- function(rnw.file, argslist=NULL, output.dir=".", output.basename=NULL){
    library(knitr)
    input.basename <- sub("\\.[rR][nN][wW]", "", rnw.file)
    if(is.null(output.basename))
        output.basename <- input.basename
    opts_chunk$set(cache = TRUE, autodep=TRUE)
    opts_chunk$set(fig.path = paste("figures/", output.basename, sep=""))
    opts_chunk$set(cache.path = paste("cache/", output.basename, sep=""))
    knit(rnw.file)
    system(paste("pdflatex ", input.basename, ".tex"))
    system(paste("pdflatex ", input.basename, ".tex"))
    system(paste("mv ", input.basename, ".pdf ", output.dir, "/", output.basename, ".pdf", sep=""))
}

##---------------------------------------------
## Principal Components Analysis for Genes
##---------------------------------------------
pca.argslist <- list(
    expression.infile=expression.infile,
    metadata.infile=metadata.infile,
    repgenes.infile="../input/75.medoids_plus_PHs.txt",
    scores.outfile=scores.outfile,
    loadings.outfile=paste(output.dir, "loadings_PhenotypeISCORE_0.4.csv", sep="/"),
    keepPPIonly=FALSE,
    var.filter.quantile=0.4,
    scale.data=TRUE,
    outcome.regex="anova",
    ##GSEA settings:
    n.top.loadings=25,
    negative.loadings.in.brackets=TRUE,
    bonferroni.threshold=0.1,
    prop.var.explained=0.5,
    n.comp=NULL)
pca.argslist$anova.regex <- ifelse(pca.argslist$keepPPIonly, "OutcomeFAP|ISCORE", "OutcomeFAP|ISCORE|Location")
##
runKnitr("PCA.Rnw", argslist=pca.argslist, output.dir=output.dir)


##---------------------------------------------
##Prepare PCL files
##---------------------------------------------
preparePCLfiles.argslist <- list(
    expression.infile=expression.infile,
    scores.outfile=scores.outfile,
    metadata.infile=metadata.infile,
    repgenes.infile="../input/75.medoids_plus_PHs.txt",
    mod.abd.infile="../input/picrust_output.pcl",
    gene.loadings.infile=pca.argslist$loadings.outfile,
    expr.all.outfile=paste(output.dir, "expr_all.csv", sep="/"),
    bugs.all.outfile=paste(output.dir, "bugs_all.csv", sep="/"),
    terminal.outfile=paste(output.dir, "/maaslin_terminal.pcl", sep=""),
    genus.outfile=paste(output.dir, "/maaslin_genus.pcl", sep=""),
    clustering.outfile=paste(output.dir, "/maaslin_clustering.pcl", sep=""),
    rel.abd.infile=paste(output.dir, "/otu_table-mod.pcl", sep=""),
    discretize.labels=c("zero", "verylow", "low", "medium", "high"),  #categories for discretizing bugs
    discretize.cutpoints=c(0, 1e-100, 1e-4, 1e-1, 0.25, 1),
    set.cutree.options=list(h=0.5),   #h defines Pearson correlation for
                                    #cutting hierarchical clustering
                                    #tree, k specifies number of
                                    #clusters.
    funcBugFilt=function(x){mean(x) > 0.005 && max(x) > 0.05},  #if TRUE, all bugs will be kept
    funcGeneFilt=function(x){TRUE},   #if TRUE, all genes will be kept
    keepPPIonly=FALSE,
    prop.var.explained=0.5)
##
runKnitr("preparePCLfiles.Rnw", argslist=preparePCLfiles.argslist, output.dir=output.dir)

##---------------------------------------------
## sources of variation as p-value histograms
##---------------------------------------------
sourcesOfVariation.argslist <- list(
    clustering.outfile=clustering.outfile,
    expr.all.outfile=expr.all.outfile,
    bugs.all.outfile=bugs.all.outfile)
##
runKnitr("sourcesOfVariation.Rnw", argslist=sourcesOfVariation.argslist, output.dir=output.dir)

##---------------------------------------------
## mas-o-menos
##---------------------------------------------
masomenos.argslist <- list(
    clustering.outfile=clustering.outfile,
    scale.data=TRUE)
runKnitr("masomenos.Rnw", argslist=masomenos.argslist, output.dir=output.dir)

##---------------------------------------------------------------------
## LDA outcome prediction
##---------------------------------------------------------------------
ldaprediction.abxall <- list(
    subset.samples=matrix(c("Location", "PPI"), ncol=2, byrow=TRUE),
    clustering.outfile=clustering.outfile,
    scale.data=TRUE)
runKnitr("ldaprediction.Rnw", argslist=ldaprediction.abxall, output.dir=output.dir, output.basename="ldaprediction_abxall")
##
ldaprediction.abxyes <- list(
    subset.samples=matrix(c("Location", "PPI", "Antibiotics", "yes"), ncol=2, byrow=TRUE),
    clustering.outfile=clustering.outfile,
    scale.data=TRUE)
runKnitr("ldaprediction.Rnw", argslist=ldaprediction.abxyes, output.dir=output.dir, output.basename="ldaprediction_abxyes")
##
ldaprediction.abxno <- list(
    subset.samples=matrix(c("Location", "PPI", "Antibiotics", "no"), ncol=2, byrow=TRUE),
    clustering.outfile=clustering.outfile,
    scale.data=TRUE)
runKnitr("ldaprediction.Rnw", argslist=ldaprediction.abxno, output.dir=output.dir, output.basename="ldaprediction_abxno")


##---------------------------------------------------------------------
## Analysis of within and between-individual differences in microbiome
## and transcriptome profiles
## ---------------------------------------------------------------------
individual.argslist <- list(
    clustering.outfile=clustering.outfile,
    expr.all.outfile=expr.all.outfile,
    bugs.all.outfile=bugs.all.outfile)
##
runKnitr("individual.Rnw", argslist=individual.argslist, output.dir=output.dir)

    
runKnitr("corpower.Rnw", output.dir=output.dir)
