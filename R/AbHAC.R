###############################################
## Documenting available data in the package ##
###############################################

#' @name ppi.database
#' @title Whole Human Protein Interaction Network
#' @description This data set is obtained with bioconductor PSICQUIC package accessing following databases: DIP,InnateDB,IntAct,MatrixDB,MINT,I2D-IMEx,InnateDB-IMEx,MolCon,BindingDB
#' @docType data
#' @usage ppi.database
#' @format First two columns as proteins interacting with each other, and subsequenc columns representing information
#' @source PSICQUIC
#' @author Mehran Karimzadeh
NULL


#' @name id.conversion.set
#' @title Table for ID Conversions inside fucntions
#' @description This data set is obtained with UniProt.ws bioconductor package getting all entrez gene IDs and their associated uniprot IDs, ensembl gene IDs, Human Gene Symbols and Human Protein Refseq IDs
#' @docType data
#' @usage id.conversion.set
#' @format Columns are Entrez Gene ID, Uniprot Acession, Human Gene Symbols with Aliases seperated by space, Human Ensembl Gene IDs and Human Protein Refseq IDs
#' @source UniProt.ss
#' @author Mehran Karimzadeh
NULL

#' @name snv
#' @title TCGA BReast Cancer Mutation Matrix for 40 Patients
#' @description This data set is obtained From TCGA Breast Cancer Dataset and is incomplete and only used as an example to run functions
#' @docType data
#' @usage snv
#' @format rownames represent gene names as uniprot accession and each column is one patient. Values of NA or "Mutation" indicate status of each patient for each gene
#' @source TCGA
#' @author Mehran Karimzadeh
NULL

#' @name rna
#' @title TCGA BReast Cancer Expression Matrix for 40 Tumors and 22 Nontumor breast tissues
#' @description This data set is obtained From TCGA Breast Cancer Dataset and is incomplete and only used as an example to run functions
#' @docType data
#' @usage rna
#' @format rownames represent gene names as uniprot accession and each column is one patient. Values are obtained from normalized Agilient Microarray. Only top 2000 most variant genes are shown. Column names with T inicate Tumor and with N indicate Nontumor. First 40 patients are in same order as mutation matrix.
#' @source TCGA
#' @author Mehran Karimzadeh
NULL

#' @name clinical
#' @title TCGA BReast Cancer clinical information table for 40 patients
#' @description This data set is obtained From TCGA Breast Cancer Dataset and is incomplete and only used as an example to run functions
#' @docType data
#' @usage clinical
#' @format First column represent patients in same order as mutation matrix and expression matrix. Subsequent columns provide some clinical information.
#' @source TCGA
#' @author Mehran Karimzadeh
NULL

#' @name fac
#' @title Proteins Annotated in ppi.database data
#' @description This is simply a vector of all proteins that an interaction for them exists in ppi.database dataset
#' @docType data
#' @usage fac
#' @format Unipriot Accession
#' @source TCGA
#' @author Mehran Karimzadeh
NULL


###############################################
##   A bunch of scripts for ID conversion    ##
###############################################

#' Uniprot Accession to HGNC symbol Conversion
#' 
#' This function is used internally to convert uniprot Accession to HGNC symbol
#' @param uniprot is a vector of Uniprot accession IDs
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a HGNC is not found, it returns the ID itself instead of NA
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
#' @examples
#' uniprot.to.hgnc(ids.to.uniprot("VHL"))
uniprot.to.hgnc = function(uniprot, id.con.set=id.conversion.set){
    # If a uniprot ID is not found, it returns the ID itself instead of ""
    if(any(uniprot %in% id.con.set$UNIPROTKB)){
        temp_idconv = id.con.set[-which(is.na(id.con.set$GENES) | is.na(id.con.set$UNIPROTKB)),]
        temp_idconv = temp_idconv[-which(duplicated(temp_idconv$UNIPROTKB)), ]
        rownames(temp_idconv) = temp_idconv$UNIPROTKB
        hgnc = temp_idconv$GENES[uniprot]
    }else{
        warning(sprintf("Uniprot to HGNC conversion failed for %s %s %s ...\n",
                        uniprot[1], uniprot[2], uniprot[3]))
        hgnc = uniprot
    }
    if(any(is.na(hgnc))){
        hgnc[which(is.na(hgnc))] = uniprot[which(is.na(hgnc))]
    }
    return(hgnc)
    # Output is a vector of same length. Unmapped IDs as NA
}


#' HGNC Symbol to Uniprot Accession Conversion
#' 
#' Converts human gene symbols to uniprot swissprot acession
#' @param hgnc is a vector of HGNC symbols
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a Uniprot accession is not found, it returns the ID itself instead of ""
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
hgnc.to.uniprot = function(hgnc, id.con.set=id.conversion.set){
    if(any(hgnc %in% id.con.set$GENES)){
        temp_idconv = id.con.set[-which(is.na(id.con.set$GENES) | is.na(id.con.set$UNIPROTKB)),]
        temp_idconv = temp_idconv[-which(duplicated(temp_idconv$GENES)), ]
        rownames(temp_idconv) = temp_idconv$GENES
        uniprot = temp_idconv$UNIPROTKB[hgnc]
    }else{
        warning(sprintf("HGNC to Uniprot conversion failed for %s %s %s ...\n",
                        hgnc[1], hgnc[1], hgnc[3]))
        uniprot = hgnc
    }
    if(any(is.na(uniprot))){
        uniprot[which(is.na(uniprot))] = hgnc[which(is.na(uniprot))]
    }
    return(uniprot)
  # Output is a vector of same length. Unmapped IDs as NA
}


#' Uniprot Accession to Entrez Gene ID conversion
#' 
#' Converts Uniprot accession to Entrez Gene ID
#' @param uniprot is a vector of uniprot accession IDs
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a Entrez accession is not found, it returns the ID itself instead of "".
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
#' @examples
#' uniprot.to.entrez(ids.to.uniprot("VHL"))
uniprot.to.entrez = function(uniprot=NULL,id.con.set=id.conversion.set){
    if(any(uniprot %in% id.con.set$UNIPROTKB)){
        temp_idconv = id.con.set[-which(is.na(id.con.set$ENTREZ_GENE) | is.na(id.con.set$UNIPROTKB)),]
        temp_idconv = temp_idconv[-which(duplicated(temp_idconv$UNIPROTKB)), ]
        rownames(temp_idconv) = temp_idconv$UNIPROTKB
        entrez = temp_idconv$ENTREZ_GENE[uniprot]
    }else{
        warning(sprintf("Uniprot to Entrez conversion failed for %s %s %s ...\n",
                        uniprot[1], uniprot[2], uniprot[3]))
        entrez = uniprot
    }
    if(any(is.na(entrez))){
        entrez[which(is.na(entrez))] = uniprot[which(is.na(entrez))]
    }
    return(entrez)
  # Output is a vector of same length. Unmapped IDs returned as input
}


#' Entrez Gene ID to Uniprot accession
#' 
#' Converts Entrez Gene ID to uniprot accessions
#' @param entrez is vector of Entrez ID characters
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a Entrez accession is not found, the return will be the ID itself
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
#' @examples
#' entrez.to.uniprot(c("12", "5421", "65223"))
entrez.to.uniprot = function(entrez=NULL,id.con.set=id.conversion.set){
    if(any(entrez %in% id.con.set$ENTREZ_GENE)){
        temp_idconv = id.con.set[-which(is.na(id.con.set$ENTREZ_GENE) | is.na(id.con.set$UNIPROTKB)),]
        temp_idconv = temp_idconv[-which(duplicated(temp_idconv$ENTREZ_GENE)), ]
        rownames(temp_idconv) = temp_idconv$ENTREZ_GENE
        uniprot = temp_idconv$UNIPROTKB[entrez]
    }else{
        warning(sprintf("Uniprot to Entrez conversion failed for %s %s %s ...\n",
                        entrez[1], entrez[2], entrez[3]))
        uniprot = entrez
    }
    if(any(is.na(uniprot))){
        uniprot[which(is.na(uniprot))] = entrez[which(is.na(uniprot))]
    }
    return(uniprot)
  # Output is a vector of same length. Unmapped IDs as NA
}



#' Ensembl Gene ID to Uniprot Accession Conversion
#' 
#' Converts ENSEMBL Human Gene IDs to Uniprot
#' @param ensembl is a vector of ensembl gene IDs
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a Uniprot accession is not found, it returns the ID itself instead of ""
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
ensembl.to.uniprot = function(ensembl=NULL,id.con.set=id.conversion.set){
    if(any(ensembl %in% id.con.set$ENSEMBL)){
        temp_idconv = id.con.set[-which(is.na(id.con.set$ENSEMBL) | is.na(id.con.set$UNIPROTKB)),]
        temp_idconv = temp_idconv[-which(duplicated(temp_idconv$ENSEMBL)), ]
        rownames(temp_idconv) = temp_idconv$ENSEMBL
        uniprot = temp_idconv$UNIPROTKB[ensembl]
    }else{
        warning(sprintf("Ensembl to Uniprot conversion failed for %s %s %s ...\n",
                        ensembl[1], ensembl[2], ensembl[3]))
        uniprot = ensembl
    }
    if(any(is.na(uniprot))){
        uniprot[which(is.na(uniprot))] = ensembl[which(is.na(uniprot))]
    }
    return(uniprot)
  # Output is a vector of same length. Unmapped IDs as NA
}

#' Refseq Protein to Uniprot Conversion Function
#' 
#' Converts Human  Protein Refseq IDs to Uniprot
#' @param refseqp is one Human Protein Refseq ID accession string
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a Uniprot accession is not found, it returns the ID itself instead of ""
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
refseqp.to.uniprot = function(refseqp=NULL,
                               id.con.set=id.conversion.set){
    if(any(refseqp %in% id.con.set$REFSEQ_PROTEIN)){
        temp_idconv = id.con.set[-which(is.na(id.con.set$REFSEQ_PROTEIN) | is.na(id.con.set$UNIPROTKB)),]
        temp_idconv = temp_idconv[-which(duplicated(temp_idconv$REFSEQ_PROTEIN)), ]
        rownames(temp_idconv) = temp_idconv$REFSEQ_PROTEIN
        uniprot = temp_idconv$UNIPROTKB[refseqp]
    }else{
        warning(sprintf("Regseqp to Uniprot conversion failed for %s %s %s ...\n",
                        refseqp[1], refseqp[2], refseqp[3]))
        uniprot = refseqp
    }
    if(any(is.na(uniprot))){
        uniprot[which(is.na(uniprot))] = refseqp[which(is.na(uniprot))]
    }
    return(uniprot)
  # Output is a vector of same length. Unmapped IDs as NA
}

#' any ID to Uniprot Conversion
#' 
#' Converts a list of IDs (only human gene symbols, ensembl gene IDs or entrez gene IDs) to Uniprot IDs
#' @param ids is a vector of one type of ID that are either ensembl gene ID, Human gene symbol or Entez gene IDs
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return If a Uniprot accession is not found, it returns the ID itself instead of ""
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
#' @examples
#' ids.to.uniprot(c("VHL","PBRM1","BRCA1","TP53"))
#' ids.to.uniprot(c("ENSG00000175793","ENSG00000170027"))
ids.to.uniprot = function(ids=NULL, id.con.set=id.conversion.set){
    temp = ids
    if(any(ids%in%id.con.set[,2])){
        print(paste("Either all your IDs are in uniprot format, or they are in different formats"))
        temp = ids
    }else if(any(ids%in%id.con.set[,1])){
        print(paste("Converting IDs from Entrez to Uniprot"))
        temp = entrez.to.uniprot(ids)
    }else if(any(ids%in%id.con.set[,3])){
        print(paste("Converting IDs from HGNC to Uniprot"))
        temp = hgnc.to.uniprot(ids)
    }else if(any(ids%in%id.con.set[,4])){
        print(paste("Converting IDs from ENSEMBL to Uniprot"))
        temp = ensembl.to.uniprot(ids)
    }else if(any(ids%in%id.con.set[,5])){
        temp = refseqp.to.uniprot(ids)
    }
    if(all(temp==ids)){
        stop("Unsuccessful Attempt for ID conversion. Consult the documentation.")
    }
    return(temp)
}


#' Mutation Matrix ID Conversion
#' 
#' Converts non-uniprot rownames of mutation matrix to uniprot. merges rows that output rowname is same.
#' IDs should be in one of these formats: entrez gene IDs, HGNC or ENSEMBL gene IDs
#' @param snv for this function is a matrix of mutations with rownames as genes and colnames as samples
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return This function finds uniprot IDs, removes rows without uniprot ID or duplicate uniprot ID
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
snv.id.conversion = function(snv=NULL
                              ### a character matrix with row names as either uniprot IDs or gene symbols or entrez IDs or ensemble gene Ids.
                              ###If there is no mutation in an entry of matrix, it should be NA and if there is a mutation it can be any character.
                              ,id.con.set=id.conversion.set
                              ### a dataframe generated from acquire.ids()
){
    id.conversion.set = id.con.set
    snv.rn = ids.to.uniprot(rownames(snv)) ## converting snv rownames to uniprot
    if(any(is.na(snv.rn))){
        snv = snv[-which(is.na(snv.rn)),]
        snv.rn = snv.rn[-which(is.na(snv.rn))]
    }
    for(j in 1:dim(snv)[2]){
        snv[,j] = as.character(snv[,j])
    }
    if(any(duplicated(snv.rn))){
        snv[which(duplicated(snv.rn,fromLast=T)|duplicated(snv.rn,fromLast=F)),] =
            t(sapply(snv.rn[which(duplicated(snv.rn,fromLast=T)|duplicated(snv.rn,fromLast=F))],function(x){
                index = which(snv.rn%in%x)
                temp = as.character(apply(snv[index,],2,function(y){
                return(paste(y,collapse="|"))
                }))
            }))
        snv = snv[-which(duplicated(snv.rn)),]
        snv.rn = snv.rn[-which(duplicated(snv.rn))]
    }
    rownames(snv) = snv.rn
    return(snv)
    ### An appropriate matrix of mutations for abhac analysis
}


#' Expression Matrix ID Conversion
#' 
#' Converts non-uniprot rownames of expression matrix to uniprot. merges rows that output rowname is same.(by averaging)
#' IDs should be in one of these formats: entrez gene IDs, HGNC or ENSEMBL gene IDs
#' @param rna for this function is a matrix of expression with rownames as genes and colnames as samples
#' @param id.con.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @return This function finds uniprot IDs, removes rows without uniprot ID or duplicate uniprot ID while averaging values
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
rna.id.conversion = function(rna=NULL,
                              ###a numeric matrix with rownames either as gene symbols, ensemble genes IDs, entrez IDs or uniprot IDs.
                              ###Column names should be the same as mutation matrix accompanied with "T" or "N" specifying the tissue.
                              ###If there are "T"s or "N"s in the name of samples, they should be changed before using the function.
                              ###If the samples are not paired, Normal samples can have different names but must be accompanied with "N"
                              id.con.set=id.conversion.set
){
    id.conversion.set = id.con.set
    if(length(which(apply((rna==0),1,all)))>0){
        rna = rna[-which(apply((rna==0),1,all)),]
    }
    averiger = function(x,rna.rn,rna){
        i=which(rna.rn%in%x)
        return(as.numeric(apply(rna[i,],2,mean)))
    }
    rna.rn = ids.to.uniprot(rownames(rna))
    if(any(is.na(rna.rn))){
        rna = rna[-which(is.na(rna.rn)),]
        rna.rn = rna.rn[-which(is.na(rna.rn))]
    }
    if(any(duplicated(rna.rn))){
        rna[which(duplicated(rna.rn,fromLast=TRUE)|duplicated(rna.rn,fromLast=FALSE)),] =
            t(sapply(rna.rn[which(duplicated(rna.rn,fromLast=TRUE)|duplicated(rna.rn,fromLast=FALSE))],function(x){
            return(averiger(x=x,rna=rna,rna.rn=rna.rn))}))
        rna = rna[-which(duplicated(rna.rn)),]
        rna.rn = rna.rn[-which(duplicated(rna.rn))]
    }
    rownames(rna) = rna.rn
    return(rna)
    ### an appropriate matrix of expression for abhac analysis
}

##################################################
##  Two scripts for creating permuted networks  ##
##################################################

#' Making Protein Edge Reoccurence dataframe
#'
#' Used inside the fisher exact test motor function for creating subsequent permuted networks.
#' @param ppi.database 2 column undirected protein interaction network.
#' @param fac vector of unique proteins inside ppi.database
#' @return A dataframe describing frequency of proteins within each edge degree.
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
get_edgedegree_freq = function(ppi.database, fac){
    ppi_temp = ppi.database
    vec_freq = table(c(ppi_temp[,1], ppi_temp[,2]))
    df_freq = data.frame(Uniprot=names(vec_freq),
                         Num.Interactions=as.numeric(vec_freq))
    rownames(df_freq) = df_freq[,1]
    return(df_freq)
}


#' Create a permuted network
#'
#' Used inside the fisher exact test motor function, it will provide permuted networks for estimating FDR/FWER.
#' @param ppi.database 2 column undirected protein interaction network
#' @param df_pr_freq dataframe of first column as unique protein names of ppi.database, and 2nd column as number of times they have occured in the network. This is calculated by get_edgedegree_freq
#' @param method Specify "equal" for equal number of proteins in each bin sorted by edge degree, "ByDegree" to create bins solely based on edge degree and "AsPaper" to have one bin for each edge degree unless the number of proteins in the edge degree is lower than 3. For such proteins create k equal size bins.
#' @param k number of bins to divide the proteins and permute inside those bins
#' @return ppi_permuted 2 column undirected protein interaction network
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
permute_pdr = function(ppi.database, df_pr_freq, method="AsPaper", k=4, verbose=FALSE){
    # Create label of groups for each protein in the network
    ppi.permute = ppi.database[,1:2]
    freq_degrees_vec = table(df_pr_freq[,2])
    if(method == "AsPaper"){
        rare_indices = as.numeric(names(freq_degrees_vec[as.numeric(freq_degrees_vec) < 4]))
        df_pr_freq$Label.bins = as.character(df_pr_freq[,2]) # as.char BCZ some must be replaced
        index_rares = df_pr_freq[,2] %in% rare_indices
        rank_degrees = rank(df_pr_freq[index_rares, 2], ties.method = "random")
        df_pr_freq$Label.bins[index_rares] = as.character(cut(
            rank_degrees, breaks = seq(0, max(rank_degrees), length.out=k+1), include.lowest=TRUE))
        df_pr_freq$Label.bins = factor(df_pr_freq$Label.bins)
    }else if(method == "ByDegree"){
        rank_degrees = rank(df_pr_freq[,2], ties.method = "min")
        df_pr_freq$Label.bins = cut(rank_degrees, k, include.lowest=TRUE, right=TRUE)
    }else if(method == "equal"){
        rank_degrees = rank(df_pr_freq[,2], ties.method = "random", include.lowest=TRUE, right=TRUE)
        df_pr_freq$Label.bins = cut(rank_degrees, k)
    }
    if(verbose){
        cat("Printing the bins of protein edge degrees\nand number of proteins within each bin\n")
        print(table(df_pr_freq$Label.bins))
        if(any(table(df_pr_freq$Label.bins) > nrow(ppi.permute)/2)){
            warning("Be advised that one of the bins includes half the proteins!\nChange the method to 'equal' or 'AsPaper'\n")
        }
    }
    # Using the label created so far, permute proteins inside those labels
    df_pr_freq$swap = NA
    for(each_label in unique(df_pr_freq$Label.bins)){
        # 1. Adding the swap column to existing dataframe df_pr_freq
        select_proteins = df_pr_freq[df_pr_freq$Label.bins==each_label, 1]
        select_proteins = setdiff(select_proteins, df_pr_freq[-which(is.na(df_pr_freq$swap)),1])
        if(length(select_proteins) > 1){
            swap_prs = sample(select_proteins, length(select_proteins), replace=FALSE)
            df_pr_freq[select_proteins, "swap"] = swap_prs
            df_pr_freq[swap_prs, "swap"] = select_proteins
            # 2. Replacing entries in the network with swap column of df_pr_freq
            in_ppi_1 = ppi.permute[,1] %in% select_proteins
            ppi.permute[in_ppi_1, 1] = sapply(ppi.permute[in_ppi_1, 1], function(each_pr){
                                               return(df_pr_freq[each_pr,"swap"])})
            in_ppi_2 = ppi.permute[,2] %in% select_proteins
            ppi.permute[in_ppi_2, 2] = sapply(ppi.permute[in_ppi_2, 2], function(each_pr){
                                               return(df_pr_freq[each_pr, "swap"])})
        }
    }
    return(ppi.permute)
}

#################################################
##  AbHAC Fisher's exact test motor function   ##
#################################################


#' Keeps only one undirected edge between two proteins
#'
#' @param ppi.database 2 column undirected protein interaction network
#' @return A protein interaction 2 column dataframe of unique undirected interactions
#' @author Mehran Karimzadeh mehran dot karimzadehreghbati at mail dot mcgill dot ca 
remove_duplicate_interactions = function(ppi.database){
    ppi_temp = unique(ppi.database)
    ppi_temp$AB = paste(ppi_temp[,1], ppi_temp[,2], sep="_")
    ppi_temp$BA = paste(ppi_temp[,2], ppi_temp[,1], sep="_")
    # Removing reversed interactions
    if(any(ppi_temp$BA %in% ppi_temp$AB)){
        ppi_temp = ppi_temp[-which(ppi_temp$BA %in% ppi_temp$AB), 1:2]
    }
    return(ppi_temp)
}

#' Corrects p-values by the provided multiple testing correction method
#'
#' @param list.pvalues 2 column dataframes of protein names and associated p-values must be put in a list structure where the 1st list member corresponds to the p-values based on the true protein interaction network and other list members if they exist must be based on permuted interaction networks
#' @param fisher.fdr Multiple testing correction method from: "Permutation.FDR", "Permutation.FWER", or p.adjust methods.
#' @param fisher.fdr.cutoff the cutoff for FWER of FDR (numeric)
#' @author Mehran Karimzadeh mehran dot karimzadehreghbati at mail dot mcgill dot ca
multiple.testing.correction.handler = function(list.pvalues, fisher.fdr, fisher.fdr.cutoff){
    n.nets = length(list.pvalues)
    if(fisher.fdr == "Permutation.FDR"){
        FDR.lists = rep(NA, n.nets - 1)
        index.fdr.lists = 1
        ##For each of the Random networks, perform the calculation seperately
        for(ppi.net in list.pvalues[2:n.nets]){
            cutoff=0.000
            min.pval = min(list.pvalues[[1]][,2]) # Minimum p-value in real network
            fdr="Not Found"
            pvalues = list.pvalues[[1]][,2] # All p-values in real network 
            pvalues = sort(pvalues[pvalues<=0.2],TRUE)
            p.ind = 1
            while(fdr=="Not Found"){
                pv.cut=pvalues[p.ind]
                n.ran = length(which(ppi.net[,2]<=pv.cut))
                n.act = length(which(list.pvalues[[11]][,2]<=pv.cut))
                if(is.na(n.ran/n.act) | is.na(n.ran/n.act)){
                    fdr="Found"
                }else if( (n.ran/n.act) <= fisher.fdr.cutoff){
                    fdr="Found"
                    cutoff = pv.cut
                }
                if(min.pval >= pv.cut){
                    fdr="Found"
                }
                  p.ind=p.ind+1
            }
            FDR.lists[index.fdr.lists] = cutoff
            index.fdr.lists = index.fdr.lists + 1
        }
        df = list.pvalues[[1]]
        if(any(is.na(FDR.lists))){
            warning('Something was wrong in FDR permutation, NA was generated for at least one of the permuted networks, but excluded')
            FDR.lists = FDR.lists[-which(is.na(FDR.lists))]
        }
        cutoff=median(FDR.lists)
        df$FDR = ifelse(df[,2] <= cutoff, 0, 1)
    }else if(fisher.fdr == "Permutation.FWER"){
        min.pvalues = lapply(list.pvalues[2:n.nets], min)
        cutoff = quantile(min.pvalues, fisher.fdr.cutoff)
        df = list.pvalues[[1]]
        df$FDR = ifelse(df[,2] <= cutoff, 0, 1)
    }else{
        df = list.pvalues[[1]]
        df$FDR = p.adjust(df[,2], method=fisher.fdr)
        cutoff = max(which(df[,2]<=fisher.fdr.cutoff))
    }
    return(df)
}

#' AbHAC Internal Enrichment Calculator
#' 
#' Internal Function that performs abhac analysis inside other functions
#' @param ppi.database 2 column whole protein interaction network. Either loaded by data(ppi.database)(filtering is recommended based on types of interactions) or bu used.
#' @param list.categories (Internal) is list of proteins for each enrichment category with accurate names as of enrichment.categories
#' @param fac is all the proteins that exist in protein interaction network. If not using data(ppi.database), it is necessary to specify.
#' @param id.conversion.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @param fisher.fdr Can be either "Permutation.FDR", "Permutation.FWER" or any of the methods parsed into p.adjust. Type ?p.adjust for more details.
#' @param fisher.fdr.cutoff Cutoff to be used for fisher's exact test false discovery rates.
#' @param num.permuted.ppi If you have selected any of the two permutation based methods, the number of permuted networks to be used for multiple testing correction must be specified.
#' @param method.permuted.ppi If you have selected any of the two permutation based methods, the method for binning proteins by their edge degree for creting permuted networks for multiple testing correction must be specified. It should be one of ("AsPaper", "equal", "ByDegree").
#' @param bins.permuted.ppi If you have selected any of the two permutation based methods, specified the number of bins for proteins to be grouped into. If you have selected "AsPaper", you would better leave this as 4 (default). For the two other methods, we advise a number between 10-20.
#' @param num.cores Note that a parallel for loop using foreach package calculates the p-values for all the different permuted networks. The number of processors to be used for this foreach has to be set by you using the "registerDoMC(cores=4)". However this parameters determined the number of processors to be used for calculation of pvalues by an mclapply feature. So if you are using registerDoMC(cores=4) and you want to limit the analysis to 12 processors, you must specify num.cores = 3.
#' @return A dataframe with Benjamini-Hochberg FDR bases on fisher's one tail exact test for all enrichment categories for all proteins that interact with at least one of the input proteins
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
Integrator = function(ppi.database=NULL,
                       ### dataframe of whole human protein interaction network that can be obtained from protein.database.creator()
                       list.categories=NULL,
                       ### a list with 7 vectors  or lists of characters as described by enrichment categories
                       fac=NULL,
                       id.conversion.set=NULL,
                       fisher.fdr="Permutation.FDR",
                       fisher.fdr.cutoff=0.05,
                       ### Parameters to be specified if using permuted method 'Permutation.FDR' or 'Permutation.FWER'
                       num.permuted.ppi = 10,
                       method.permuted.ppi = "AsPaper",
                       bins.permuted.ppi = 4,
                       num.cores = 6)
{
    if(is.null(fac)){
        fac = AbHAC::fac
    }
    if(is.null(ppi.database)){
        ppi.database = AbHAC::ppi.database[,1:2]
    }
    if(is.null(id.conversion.set)){
        id.conversion.set = AbHAC::id.conversion.set
    }
    # 1. Removing duplicate interactions in any direction
    ppi.database = remove_duplicate_interactions(ppi.database)
    # 2. computations will be done on all the protein networks within ppi.lists and will be stored to list.pvalues
    ppi.lists = list()
    if(grepl("ermut", fisher.fdr)){
        ppi.lists = vector("list", (num.permuted.ppi + 1))
    }
    ppi.lists[[1]] = ppi.database
    names(ppi.lists)[1] = "MainNetwork"
    list.results = list()
    df_pr_freq = get_edgedegree_freq(ppi.database, fac)
    for(l in 1:length(list.categories)){
        ###Fina all the proteins that interact with at least one of the proteins in list.categories[[i]]
        list.pvalues = foreach(p=1:(num.permuted.ppi + 1)) %dopar% {
            ppi.dat = ppi.lists[[p]]
            if(is.null(ppi.dat)){
                ppi.dat = permute_pdr(ppi.database, df_pr_freq, method=method.permuted.ppi, k=bins.permuted.ppi)
            }
            case = list.categories[[l]]
            if(is.list(case)){
                for(z in 1:length(case)){
                    case[[z]] = intersect(fac,case[[z]])
                }
                all.case = unique(unlist(lapply(case,function(x){return(x)}))) #
                possible.prs = unique(union(ppi.dat[which(ppi.dat[,1]%in%case[[1]]),2],
                                     ppi.dat[which(ppi.dat[,2]%in%case[[1]]),1])) ##all proteins interacting with first list
                temp.ppi = ppi.dat[which(ppi.dat[,1]%in%possible.prs),] ##ppi dataframe of these prs
                temp.ppi.2 = ppi.dat[which(ppi.dat[,2]%in%possible.prs),] ##ppi dataframe of these prs
                temp.ppi.2 = temp.ppi.2[,2:1] 
                colnames(temp.ppi.2) = colnames(temp.ppi)
                temp.ppi = rbind(temp.ppi,temp.ppi.2)
                reference.proteins = union(
                    temp.ppi[which(temp.ppi[,1]%in%case[[2]]),2],
                        temp.ppi[which(temp.ppi[,2]%in%case[[2]]),1]) ##those proteins in previous ppi dataframe interacting with list 2
                if(length(case)>2){
                    temp.ppi.3 = temp.ppi[which(temp.ppi[,1]%in%reference.proteins),]
                    temp.ppi.2 = temp.ppi[which(temp.ppi[,2]%in%reference.proteins),]
                    temp.ppi.2 = temp.ppi.2[,2:1]
                    colnames(temp.ppi.2) = colnames(temp.ppi)
                    temp.ppi.3 = rbind(temp.ppi.3,temp.ppi.2)
                    reference.proteins = unique(union(temp.ppi.3[which(temp.ppi.3[,1]%in%case[[3]]),2],temp.ppi.3[which(temp.ppi.3[,2]%in%case[[3]]),1]))
                }
                case = all.case
            }else{
                case = intersect(case,fac)
                reference.proteins = unique(union(ppi.dat[which(ppi.dat[,1]%in%case),2],ppi.dat[which(ppi.dat[,2]%in%case),1]))
            }
            if(.Platform$OS.type!="windows"){
                df = unlist(mclapply(reference.proteins,function(protein){
                    interactors = unique(union(as.character(ppi.dat[which(ppi.dat[,1]%in%protein),2]),
                                          as.character(ppi.dat[which(ppi.dat[,2]%in%protein),1])))
                    A=length(intersect(interactors,case))##interactors of protein in category l
                    C=length(setdiff(case,interactors))##genes of category l not among interactors
                    B=length(setdiff(interactors,case))##interactors not in category i
                    D=length(unique(union(ppi.dat[,1],ppi.dat[,2]))) - (A+B+C)
                    p.value = fisher.test(matrix(c(A,C,B,D),nrow=2),alternative="greater")[[1]]
                    result = p.value
                    return(result)
                    },mc.cores=num.cores))
                df = data.frame(reference.proteins,as.numeric(df))
            }else{
                df = sapply(reference.proteins,function(protein){
                    interactors = unique(union(as.character(ppi.dat[which(ppi.dat[,1]%in%protein),2]),
                                          as.character(ppi.dat[which(ppi.dat[,2]%in%protein),1])))
                    A=length(intersect(interactors,case))##interactors of protein in category l
                    C=length(setdiff(case,interactors))##genes of category l not among interactors
                    B=length(setdiff(interactors,case))##interactors not in category i
                    D=length(fac) - (A+B+C)
                    p.value = fisher.test(matrix(c(A,C,B,D),nrow=2),alternative="greater")[[1]]
                    result = p.value
                    return(result)
                    })
                df = data.frame(reference.proteins,as.numeric(df))
            }
            colnames(df) = c("Protein",paste("P.Value",names(list.categories)[l],sep="_"))
            df = as.data.frame(df)
            df[,2] = as.numeric(df[,2])
            df
        }
        df = multiple.testing.correction.handler(list.pvalues, fisher.fdr, fisher.fdr.cutoff)
        list.results = c(list.results,list(df))
    }
    names(list.results) = names(list.categories)
    ###
    all.genes = unique(unlist(lapply(list.results,function(x){
        return(as.character(x[,1]))
    })))
    result = data.frame(all.genes)
    colnames(result) = "Protein"
    temp.tbl = matrix(1,nrow=length(all.genes),ncol=length(list.results))
    colnames(temp.tbl) = paste("FDR",names(list.results),sep="")
    for(j in 1:length(colnames(temp.tbl))){
        list.temp = list.results[[j]]
        ord.ind = numeric()
        for(J in 1:nrow(list.temp)){
            ord.ind = c(ord.ind,which(result[,1]%in%list.temp[J,1]))
        }
        temp.tbl[ord.ind,j] = list.temp[,3]
    }
    result = cbind(result,temp.tbl)
    temp.tbl = matrix(1,nrow=length(all.genes),ncol=length(list.results))
    colnames(temp.tbl) = paste("p.value",names(list.results),sep="")
    for(j in 1:length(colnames(temp.tbl))){
        list.temp = list.results[[j]]
        ord.ind = numeric()
        for(J in 1:nrow(list.temp)){
          ord.ind = c(ord.ind,which(result[,1]%in%list.temp[J,1]))
        }
        temp.tbl[ord.ind,j] = list.temp[,2]
    }
    result = cbind(result,temp.tbl)
    result$HGNC = uniprot.to.hgnc(result[,1])
    return(result)
    ###Output is is a dataframe listing FDRs for each protein in each enrichment category
}

#' Differential Expression calculator function
#' @param rna is expression matrix, method of determined expression, and index for samples to be included and pvalue adjustment method paired, a boolean, indicate if for each sample, we should expect adjacent normal tissue or not
#' @param i is index of columns to calculate differential expression agains columns having N
#' @param paired is a boolean showing if RNA samples are paired
#' @param expression.method is a string either "Microarray" or "RNAseq"
#' @param correction.method Any of the "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr" can be used.
#' @param fdr.cutoff 0.01 or 0.05 are suggested. All values in renge of 0 and 1 are accepted
#' @param fac is all the proteins that exist in protein interaction network. If not using data(ppi.database), it is necessary to specify.
#' @return a list with 3 string vectors: de.up, de.down and de
#' @export
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
deizer = function(rna=NULL,
                   ###a numeric matrix with rownames either as gene symbols, ensemble genes IDs, entrez IDs or uniprot IDs.
                   ###Column names should be the same as mutation matrix accompanied with "T" or "N" specifying the tissue.
                   ###If there are "T"s or "N"s in the name of samples, they should be changed before using the function.
                   ###If the samples are not paired, Normal samples can have different names but must be accompanied with "N"
                   expression.method=NULL,
                   ### "RNAseq" or "Microarray"
                   i=NULL,
                   ### index indicating which Tumor patients from columns of expression matrix shall be included in the analysis
                   paired=NULL,
                   ### If for each tumor exists a normal or nor (a boolean)
                   correction.method=NULL,
                   ### Any of the "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr" can be used.
                   fdr.cutoff=NULL,
                   ### 0.01 or 0.05 are suggested. All values in range of 0 and 1 are accepted
                   fac=NULL
){
  if(is.null(fac)){
    fac = AbHAC::fac
  }
  if(is.null(fdr.cutoff)){
    fdr.cutoff=0.01
  }
  if(is.null(correction.method)){
    correction.method="BH"
  }
  if(is.null(i)|is.null(paired)|is.null(rna)|is.null(expression.method)){
    stop("Please enter valid arguments for the function")
  }
  patients_raw = gsub("T", "", colnames(rna)[grep("T", colnames(rna))])
  rna = rna[which(rownames(rna)%in%fac),]
  if(expression.method=="Microarray"){
    if(paired) {
      index = 1:length(i)*2
      z=1
      for(test in 1:length(i)) {
        index[(z:(z+1))] = grep(patients_raw[test],colnames(rna))
        z=z+2
      }
      x = rna[,index]
      patient.stat = rep(1:(length(index)/2),2)
      status=rep(1,length(colnames(x)))
      status[grep("N",colnames(x))] = 2
      design = model.matrix(~patient.stat+status)
    } else  {
      index = which(gsub("T","",colnames(rna))%in% patients_raw[i])
      x = rna[ , c(index,grep("N",colnames(rna)))]
      status=rep(1, length(colnames(x)))
      status[grep("N",colnames(x))] = 2
      design = model.matrix(~ 0+factor(status))
    }        
    colnames(design) = c("Tumor","Normal")
    fit = lmFit(x,design)
    contrast.matrix = makeContrasts(Tumor-Normal,levels=design)
    fit2 = contrasts.fit(fit,contrast.matrix)
    fit2 = treat(fit2,lfc=1)
    deg = topTreat(fit2,number=Inf)
    deg$adj = p.adjust(deg$P.Value,method=correction.method)
    de.up = rownames(deg)[which(deg$adj<=fdr.cutoff & deg$logFC>0)]
    de.down = rownames(deg)[which(deg$adj<=fdr.cutoff & deg$logFC<0)]
    de = union(de.up,de.down)
    
  } else if (expression.method=="RNAseq") {
    if(paired) {
      index = 1:length(i)*2
      z=1
      for(test in 1:length(i)){
        index[(z:(z+1))] = grep(patients_raw[test],colnames(rna))
        z=z+2
      }
      x = rna[,index]
    }else {
      index = i
      x = rna[,c(index,grep("N",colnames(rna)))]
      ##BECAREFUL, NOT COMPLETE FOR PAIRED SAMPLES
    }
    group = rep(2,length(colnames(x)))
    group[grep("N",colnames(x))] = 1
    group = factor(group)
    y = DGEList(counts=x,group=group)
    y = calcNormFactors(y)
    y = estimateCommonDisp(y)
    y = estimateTagwiseDisp(y)
    et = exactTest(y)
    table = et$table
    table$FDR = p.adjust(table$PValue,method=correction.method)
    table = table[which(table$FDR<=fdr.cutoff),]
    de.up = rownames(table)[which(table[,1]>=1)]
    de.down = rownames(table)[which(table[,1]<=(-1))]
    de = union(de.up,de.down)
  }
  result=list(de.up,de.down,de)
  names(result) = c("de.up","de.down","de")
  print(paste("Differential Expression Analysis Successful",Sys.time()))
  return(result)
  ### a list with upregulated proteins, downregulated proteins, and all differentially expressed proteins
}

#' This Function performs integrative analysis on the whole dataset
#' @param snv a character matrix with row names as either of uniprot IDs, gene symbols, entrez IDs or ensemble gene Ids. When there is no mutation, NA is expected.
#' @param rna a numeric matrix with rownames either as gene symbols, ensemble genes IDs, entrez IDs or uniprot IDs and columns representing same patient IDs as mutation matrix but with T and N suffices standing for Tumor and Nontumor samples.
#' @param expression.method "RNAseq" or "Microarray"
#' @param correction.method for differential expression analysis "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr" 
#' @param ppi.database 2 column whole protein interaction network. Either loaded by data(ppi.database)(filtering is recommended based on types of interactions) or by user.
#' @param rna.paired boolean indicating if RNA samples are paired
#' @param fdr.cutoff values higher than 0.1 are not advised
#' @param clinical a dataframe of 2 columns with first column identical to snv dataframe or RNA dataframe (without "T"). Second column provides a particular categrical phenodata of patients that will be used as a factor to perform the analysis within the factor levels independently.
#' @param enrichment.categories can be all or any of the c("snv.de.up.de.down","de.up","de.down","de","snv","snv.de.up","snv.de.down","snv.de")
#' @param id.conversion.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @param fac is all the proteins that exist in protein interaction network. If not using data(ppi.database), it is necessary to specify.
#' @param fisher.fdr Can be either "Permutation.FDR", "Permutation.FWER" or any of the methods parsed into p.adjust. Type ?p.adjust for more details.
#' @param fisher.fdr.cutoff Cutoff used for false discovery rate cutoff in fisher's exact test. By default set to 0.2.
#' @param num.permuted.ppi If you have selected any of the two permutation based methods, the number of permuted networks to be used for multiple testing correction must be specified.
#' @param method.permuted.ppi If you have selected any of the two permutation based methods, the method for binning proteins by their edge degree for creting permuted networks for multiple testing correction must be specified. It should be one of ("AsPaper", "equal", "ByDegree").
#' @param bins.permuted.ppi If you have selected any of the two permutation based methods, specified the number of bins for proteins to be grouped into. If you have selected "AsPaper", you would better leave this as 4 (default). For the two other methods, we advise a number between 10-20.
#' @param num.cores Note that a parallel for loop using foreach package calculates the p-values for all the different permuted networks. The number of processors to be used for this foreach has to be set by you using the "registerDoMC(cores=4)". However this parameters determined the number of processors to be used for calculation of pvalues by an mclapply feature. So if you are using registerDoMC(cores=4) and you want to limit the analysis to 12 processors, you must specify num.cores = 3.
#' @return This function returns a data.frame with results of integrative network enrichment analysis
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
#' @examples
#' data(snv)
#' data(rna)
#' data(ppi.database) #2column whole human protein interaction database
#' data(id.conversion.set) 
#' data(fac) #vector of all proteins in ppi.database
#' set.abhac.result = set.abhac(snv=snv,rna=rna,fac=fac,expression.method="Microarray",rna.paired=FALSE,
#'    fdr.cutoff=0.05,correction.method="BH",enrichment.categories=c("snv.de","de.up"),
#'    ppi.database=ppi.database[,1:2],id.conversion.set=id.conversion.set)
set.abhac = function(ppi.database=NULL,
                     ### dataframe of whole human protein interaction network that can be obtained from protein.database.creator()
                     rna=NULL,
                     ###a numeric matrix with rownames either as gene symbols, ensemble genes IDs, entrez IDs or uniprot IDs.
                     ###Column names should be the same as mutation matrix accompanied with "T" or "N" specifying the tissue.
                     ###If there are "T"s or "N"s in the name of samples, they should be changed before using the function.
                     ###If the samples are not paired, Normal samples can have different names but must be accompanied with "N"
                     snv=NULL,
                     ### a character matrix with row names as either uniprot IDs or gene symbols or entrez IDs or ensemble gene Ids.
                     ###If there is no mutation in an entry of matrix, it should be NA and if there is a mutation it can be any character.
                     expression.method=NULL,
                     ### "RNAseq" or "Microarray"
                     rna.paired=FALSE,
                     ### a boolean indicating if adjacent normal tissue expression is available for each tumor tissue column.
                     correction.method="BH",
                     ### "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr"}
                     fdr.cutoff=0.05,
                     ### 0.01 or 0.05 are suggested. All values in range of 0 and 1 are accepted
                     enrichment.categories=NULL,
                     ### a string vector for any of these variables; "snv.in","de.up","de.down","de","snv","snv.de.up","snv.de.down"
                     id.conversion.set=NULL,
                     fac=NULL,
                     clinical=NULL,
                     fisher.fdr="Permutation.FDR",
                     fisher.fdr.cutoff=0.05,
                     num.permuted.ppi=10,
                     method.permuted.ppi="AsPaper",
                     bins.permuted.ppi=4,
                     num.cores=6)
{
  if(is.null(ppi.database)){
    ppi.database = AbHAC::ppi.database[,1:2]
  }
  if(is.null(fac)){
    fac = AbHAC::fac
  }
  if(is.null(id.conversion.set)){
    id.conversion.set = AbHAC::id.conversion.set
  }
  if(is.null(rna) & is.null(snv)){
    stop("Please make sure all the arguments of the function are entered correctly")
  }
  if(!is.null(rna) & !is.null(snv)){
    if(!all(colnames(snv)==gsub("T","",colnames(rna)[grep("T",colnames(rna))]))){
      stop("Make sure column names of snv and rna are the same and an additional T exists on each RNA columnname")
    }
  }
  if(!is.null(rna)){
    if(!any(grepl("N",colnames(rna)))){
      stop("Normal RNA samples should have a 'N' suffix")
    } 
    if(!any(rownames(rna)%in%id.conversion.set[,2])){
      print(paste("rna matrix IDs are not uniprot, trying for conversion, averiging duplicates, removing those without IDs"))
      rna = rna.id.conversion(rna,id.conversion.set)
    }
  }
  if(!is.null(snv)){
    if(!any(rownames(snv)%in%id.conversion.set[,2])){
      print(paste("snv matrix IDs are not uniprot, trying for conversion, averiging duplicates, removing those without IDs"))
      snv = snv.id.conversion(snv,id.conversion.set)
    }
  }
  if(is.null(clinical)){
      if(is.null(snv)){
        clinical = data.frame(Patient=gsub("T", "", colnames(rna)[grep("T", colnames(rna))]),
                               Status="Tumor")
      }else{
        clinical = data.frame(Patient=colnames(snv),Status="Tumor")
      }
  }
  #based on de method, find differentially expressed genes
  grs = unique(clinical[,2])
  list.abhacs = vector('list',length(grs))
  names(list.abhacs) = grs
  lcind=1
  for(gr in grs){
    cat(paste("Performing AbHAC analysis for:", gr, "\n"))
    i=which(clinical[,2]==gr)
    if(!is.null(rna)){
      DE = deizer(rna=rna,correction.method=correction.method,
		paired=rna.paired,expression.method=expression.method,
		i=i,fdr.cutoff=fdr.cutoff,
		fac=fac)
      de.up = DE$de.up
      de.down = DE$de.down
      de = DE$de
    }
    if(!is.null(snv)){
      SNV = rownames(snv)[apply(!is.na(snv[,i]),1,sum)>0]
    }
    list.categories = list()
    if(!is.null(snv) & !is.null(rna)){
      list.categories[[1]] = c(list(SNV),list(de.up),list(de.down))
      list.categories[[2]] = de.up
      list.categories[[3]] = de.down
      list.categories[[4]] = union(de.up,de.down)
      list.categories[[5]] = SNV
      list.categories[[6]] = c(list(SNV),list(de.up))
      list.categories[[7]] = c(list(SNV),list(de.down))
      list.categories[[8]] = c(list(SNV),list(de))
      names(list.categories) = c("snv.de.up.de.down","de.up","de.down",
                                "de","snv","snv.de.up","snv.de.down","snv.de")
    }else if(!is.null(snv)){
      list.categories = list(snv=SNV)
    }else{
      list.categories = list(de.up = de.up,
        de.down = de.down,
        de = de)
    }
    index.cat = numeric()
    for(z in 1:length(enrichment.categories)){
      index.cat = c(index.cat,which(names(list.categories)%in%enrichment.categories[z]))
    }
    list.categories = list.categories[index.cat] 
    list.categories = list.categories[which(names(list.categories)%in%enrichment.categories)]
    if(length(enrichment.categories)!=length(list.categories)){
      stop("Make sure proper enrichment categories are queried according to the provided datasets")
    }
    ## Local Function that gets snv.in, de.up, de.down, de and SNV and returns a dataframe of results
    abhac = Integrator(ppi.database=ppi.database, list.categories=list.categories,
	fac=fac, fisher.fdr=fisher.fdr,
	fisher.fdr.cutoff=fisher.fdr.cutoff,
        num.permuted.ppi=num.permuted.ppi,
        method.permuted.ppi=method.permuted.ppi,
        bins.permuted.ppi=bins.permuted.ppi,
        num.cores=num.cores)
    list.abhacs[[lcind]] = abhac
    lcind=lcind+1
  }
  return(list.abhacs)
  ### Returns a dataframe with each protein and enrichments in 3 different categories
}


#' AbHAC Based on Vector Inputs
#' This function performs AbHAC analysis on inputs of mutations, upregulated and downregulated proteins
#' @param de.up either ensembl gene IDs, HGNC, ENTREZ or uniprot IDs for a vector of upregulated genes 
#' @param de.down either ensembl gene IDs, HGNC, ENTREZ or uniprot IDs for a vector of downregulated genes
#' @param snv either ensembl gene IDs, HGNC, ENTREZ or uniprot IDs for a vector of mutated genes
#' @param ppi.database 2 column whole protein interaction network. Either loaded by data(ppi.database)(filtering is recommended based on types of interactions) or by user.
#' @param fac is all the proteins that exist in protein interaction network. If not using data(ppi.database), it is necessary to specify.
#' @param enrichment.categories can be all or any of the c("snv.de.up.de.down","de.up","de.down","de","snv","snv.de.up","snv.de.down","snv.de")
#' @param id.conversion.set A dataframe for ID conversions provided as global variable id.conversion.set. Columns represent Entrez gene ID, Uniprot Accession, Gene Symbol, Ensembl gene ID and refseq protein ID (all human)
#' @param fisher.fdr Can be either "Permutation.FDR", "Permutation.FWER" or any of the methods parsed into p.adjust. Type ?p.adjust for more details.
#' @param fisher.fdr.cutoff Cutoff used for false discovery rate cutoff in fisher's exact test. By default set to 0.2.
#' @param num.permuted.ppi If you have selected any of the two permutation based methods, the number of permuted networks to be used for multiple testing correction must be specified.
#' @param method.permuted.ppi If you have selected any of the two permutation based methods, the method for binning proteins by their edge degree for creting permuted networks for multiple testing correction must be specified. It should be one of ("AsPaper", "equal", "ByDegree").
#' @param bins.permuted.ppi If you have selected any of the two permutation based methods, specified the number of bins for proteins to be grouped into. If you have selected "AsPaper", you would better leave this as 4 (default). For the two other methods, we advise a number between 10-20.
#' @param num.cores Note that a parallel for loop using foreach package calculates the p-values for all the different permuted networks. The number of processors to be used for this foreach has to be set by you using the "registerDoMC(cores=4)". However this parameters determined the number of processors to be used for calculation of pvalues by an mclapply feature. So if you are using registerDoMC(cores=4) and you want to limit the analysis to 12 processors, you must specify num.cores = 3.
#' @return This function returns a data.frame with results of integrative network enrichment analysis
#' @author Mehran Karimzadeh mehran.karimzadehreghbati at mail dot mcgill dot ca
#' @export
#' @examples
#' data(snv)
#' data(rna)
#' snv = sample(rownames(snv),10)
#' de.up = sample(rownames(rna)[1:1000],500)
#' de.down = sample(rownames(rna)[1001:2000],500)
#' data(ppi.database) #2column whole human protein interaction database
#' data(id.conversion.set) 
#' data(fac) #vector of all proteins in ppi.database
#' abhac.brief.result = abhac.brief(de.up,de.down,fac=fac,snv=snv,enrichment.categories=c("snv.de","de.up"),ppi.database=ppi.database[,1:2],id.conversion.set=id.conversion.set)
abhac.brief = function(de.up=NULL,
                       ### Uniprot acession ID of upregulated genes
                       de.down=NULL,
                       ### Uniprot accession ID of downregulated proteins
                       snv=NULL,
                       ### Uniprot accession ID of mutated genes
                       ppi.database=NULL,
                       ### dataframe of whole human protein interaction network that can be obtained from protein.database.creator()
                       enrichment.categories=NULL,
                       fac=NULL,
                       fisher.fdr="Permutation.FDR",
                       fisher.fdr.cutoff=0.05,
                       id.conversion.set=NULL,
                       num.permuted.ppi=10,
                       method.permuted.ppi="AsPaper",
                       bins.permuted.ppi=4,
                       num.cores=6){
  if(is.null(de.up) & is.null(de.down) & is.null(snv)){
    stop("Please make sure the 3 main variables (de.up, de.down, snv) are entered correctly")
  }
  if(is.null(ppi.database)){
    ppi.database = AbHAC::ppi.database[,1:2]
  }
  if(is.null(id.conversion.set)){
    id.conversion.set = AbHAC::id.conversion.set
  }
  if(is.null(fac)){
    fac = AbHAC::fac
  }
  if(!is.null(de.up)){
    if(length(intersect(de.up, id.conversion.set[,2]))==0){
      de.up = ids.to.uniprot(de.up)
    }
  }
  if(!is.null(de.down)){
    if(length(intersect(de.down, id.conversion.set[,2]))==0){
      de.up = ids.to.uniprot(de.down)
    }
  }
  if(!is.null(snv)){
    if(length(intersect(snv, id.conversion.set[,2]))==0){
      snv = ids.to.uniprot(snv)
    }
  }

  id.conversion.set. = id.conversion.set
  list.categories = list()
  list.categories[[1]] = c(list(snv),list(de.up),list(de.down))
  list.categories[[2]] = de.up
  list.categories[[3]] = de.down
  list.categories[[4]] = union(de.up,de.down)
  list.categories[[5]] = snv
  list.categories[[6]] = c(list(snv),list(de.up))
  list.categories[[7]] = c(list(snv),list(de.down))
  list.categories[[8]] = c(list(snv),list(union(de.up,de.down)))
  names(list.categories) = c("snv.de.up.de.down","de.up","de.down",
                              "de","snv","snv.de.up","snv.de.down","snv.de")
  index.cat = numeric()
  for(z in 1:length(enrichment.categories)){
    index.cat = c(index.cat,which(names(list.categories)%in%enrichment.categories[z]))
  }
  list.categories = list.categories[index.cat] 
  list.categories = list.categories[which(names(list.categories)%in%enrichment.categories)]
  abhac = Integrator(ppi.database=ppi.database[,1:2],list.categories=list.categories,
                     fac=fac,id.conversion.set=id.conversion.set.,
                     fisher.fdr=fisher.fdr,fisher.fdr.cutoff=fisher.fdr.cutoff,
                     num.permuted.ppi=num.permuted.ppi,
                     method.permuted.ppi=method.permuted.ppi,
                     bins.permuted.ppi=bins.permuted.ppi,
                     num.cores=num.cores)
  return(abhac)
  ### Returns a dataframe with each protein, all its interactors, and enrichments in 3 different categories, before and after correction for multiple testing
}
