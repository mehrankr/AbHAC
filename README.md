# AbHAC R package 

## AbHAC: Aberration Hub analysis of Cancer


AbHAC is an R package for implementation of a simplistic approach for analysis of 
cancer genomics datasets in context of protein interaction networks. 
In AbHAC, each protein in the protein interaction network is considered as an individual subnetwork 
and based on abundance of molecules with aberrations at genomic
or transcriptomic levels among neighborhood of that protein as well as the whole interactome, a Fisher's exact test
p-value is calculated. The Fisher's exact test p-values are then corrected for multiple testing
by permutation of the protein interaction network. Details are available in the [paper](http://www.abhac.com/)/[thesis](http://www.abhac.com/) manuscript.



## Usage: _Required objects_

__snv__: a matrix/dataframe where column names represent name of samples
and rownames represent name of genes. Value of each cell can either be NA or a character (e,g. "Mutated")

__rna__: a numeric matrix/dataframe where column names are sample names similar and in same order with snv.
However, these names must be accompanied with __T__ at the end of their name. 
For example, if sample names in _snv_ are: a | b | c ..., in _rna_
they should be: aT | bT | cT ... . These must be followed with nontumour samples ending with N.
It is possible for the nontumour samples to have the same name (aN | bN ...) or something different (a2313N | a321bchN).


## Usage: _Optional objects_

__clinical__: A dataframe with first column having the same names as _snv_, and the second column providing information about samples.
These can be Metastasis/Primary, HighGrade/LowGrade or any other sets of strings describing patients subtypes.


__ppi.database__: By default, the package uses a protein interaction network built using [PSICQUIC](http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml) 
by querying for uniprot accession IDs obtained through [Uniprot.ws package](http://www.bioconductor.org/packages/release/bioc/html/UniProt.ws.html).
The databases used for generating this dataframe include:

DIP, InnateDB, IntAct, MatrixDB, MINT, I2D-IMEx, InnateDB-IMEx, MolCon and BindingDB


The AbHAC functions require the first 2 columns of this dataframe. The IDs must be uniprot accession.


__id.conversion.set__: A dataframe with the following column names:

| ENTREZ_GENE | UNIPROTKB | GENES   | ENSEMBL        | REFSEQ_PROTEIN    |
| ----------- |:---------:|:-------:|:--------------:| -----------------:|
| 7533        | Q04917    | YWHAH   |ENSG00000128245 | NP_003396         |



## Usage: _Examples_

Installing the package and all of its dependencies:

```R
install.packages(c("devtools", "foreach", "doMC", "iterators" ,"plyr"))
source("http://bioconductor.org/biocLite.R")
biocLite("EdgeR")
require(devtools)
install_github("AbHAC", username="mehrankr")
require("AbHAC")	
```
_abhac.brief_ is implemented to be used when a particular set of genes are of interest and we
want to investigate the proteins that might interact with a significant number of our set of genes.
These set of genes might be mutated (snv), upregulated (de.up) or downregulated(de.down).

Running _abhac.brief_ with vector of mutated/upregulated/downregulated genes:

```R
#Loading matrix of mutated genes and matrix of mRNA expression
data(snv)
data(rna)

#Randomly selecting the first 10/1000 genes
snv = sample(rownames(snv), 10)
de.up = sample(rownames(rna)[1:1000], 500)
de.down = sample(rownames(rna)[1001:2000], 500)

#Loading the default protein interaction data
data(ppi.database) 

#Loading dataframe used for converting IDs
data(id.conversion.set)

#Loading _fac_ which is a vector of all proteins existing inside _ppi.database_
data(fac) #vector of all proteins in ppi.database

abhac.brief.result = abhac.brief(de.up,de.down,fac=fac,snv=snv,
	enrichment.categories=c("snv.de","de.up"),
	ppi.database=ppi.database[,1:2],
	id.conversion.set=id.conversion.set)
```

If instead of particular selections of differentially expressed genes, we have 
an RPKM matrix of RNAseq or normalized mRNA expression values,
AbHAC can find differentially expressed genes using EdgeR/limma. This is through
the _set.abhac_ function which accepts _snv_ and _rna_ matrices as input.
The other important feature of this function is that you can provide
subtype / phenodata of patients in a two column object called _clinical_.
```R
#Loading example and default objects from the package
data(snv)
data(rna)
data(ppi.database) #2column whole human protein interaction database
data(id.conversion.set)
data(fac) #vector of all proteins in ppi.database


set.abhac.result = set.abhac(snv=snv,rna=rna,fac=fac,
   expression.method="Microarray",rna.paired=FALSE,
   fdr.cutoff=0.05,correction.method="BH",enrichment.categories=c("snv.de","de.up"),
   ppi.database=ppi.database[,1:2],id.conversion.set=id.conversion.set)

```

## Important parameters

_fisher.fdr_ : This parameter which is defaulted to using the permutation method described in the paper,
can be set to any of the parameters accepted by _p.adjust_.
The permutation based methods include _Permutation.FDR_ and _Permutation.FWER_.
if selecting any of these methods, other parameters described below would be important.


_fisher.fdr.cutoff_ : By default is set to 0.05.


_num.permuted.ppi_: Number of permuted protein interaction networks to generate for multiple testing correction.


_method.permuted.ppi_: There are three options: __AsPaper__, __ByDegree__ or __equal__.

* __AsPaper__ This method assumes all of the edge degrees with more than 4 proteins as a bin of proteins that proteins
will be permuted inside that bin. However for edge degrees with less than 4 proteins, they all will be groupd into k
bins defined by _bins.permuted.ppi_.

* __ByDegree__ Proteins will be grouped into k categories determined by _bins.permuted.ppi_ without breaking the edge degrees.
The bins created would have varying numbers; Some very low, some very high.

* __equal__ This would create bins with equal size of proteins. For bins that have thousand of proteins, it randomly distributes them to closest bins. It does this by ranking the proteins according to their edge degree and using the "random" method of ties.method in R rank function.


_bins.permuted.ppi_: Number of bins that proteins in the network are categorized into and then permuted within those bins. Read parameter specified by _method.permuted.ppi_ to understand more.


### Nomenclature:

In [old irish](http://en.wiktionary.org/wiki/abhac), abhac means a dwarf star. 


__Maintainer__: mehran dot karimzadeh at uhnresearch dot ca or mehran dot karimzadehreghbati at mail dot mcgil dot ca

