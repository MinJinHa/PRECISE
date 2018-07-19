getregcovDat <- function(Gmat,RPPAdat,mRNA=NULL,miRNA=NULL,Methylation=NULL,anno.miRNA=NULL,anno.methyl=NULL,integrative=T){
    #----------------------------------------
    # This function is to obtain a RPPA response vectors and covariate matrices from mRNA, miRNA, methylation data, and protein regulations from the estimated graph (G).
    # Input :
    #        - Gmat : the p by p adjacency matrix (row -> column) which includes weight
    #        - RPPAdat : n times p matrix
    #        - mRNA : a list including Data and Des
    #        - miRNA : a list including Data and Des
    #        - Methylation : a list including Data and Des
    #        - anno.miRNA : annotation matrix for miRNA including gene name in the first column and miRNA names in the second column
    #        - anno.methyl : annotation matrix for Methylation including gene name in the first column and Methylation names in the second column
    #        - integrative : TRUE if use miRNA, Methylation, and mRNA
    # Note: all samples should be matched. The rows of RPPAdata and columns of mRNA$Data, miRNA$Data, Methylation$Data should be matched.
    #----------------------------------------
    
    n = nrow(RPPAdat)
    p = ncol(RPPAdat)
    
    nodes = colnames(RPPAdat)
    genelist = strsplit(nodes,split=", ")
    genes = unlist(genelist)
    membership = rep(1:p,lapply(genelist,length))
    
    if (integrative) {
        # Indices for miRNA #
        if (!is.null(miRNA$Data)) {
            w.miRNA = grepVecAnno(x=genes,y=miRNA$Des[,1],anno=anno.miRNA)
        }
        if (!is.null(mRNA$Data)) {
            w.mRNA = match(genes,mRNA$Des[,1])
        }
    }
    
    ylist = Xlist  = vector("list",length=p)
    names(ylist) = names(Xlist) = nodes
    for (i in 1:p) {
        gg = genelist[[i]]
        ylist[[i]] = scale(RPPAdat[,i])
        
        X = numeric(0)
        typevar = numeric(0)
        # include regulator proteins in X
        pa = which(Gmat[,i] != 0)
        if (length(pa)>0) {
            X = RPPAdat[,pa]
            typevar=c(typevar,paste("pa_",nodes[pa],sep=""))
        }
        if (integrative) {
            if (!is.null(miRNA$Data)) {
                # include miRNA in X
                w = which(membership==i)
                for (j in 1:length(w)) {
                    if (length(w.miRNA[[w[j]]])>0) {
                        X = cbind(X,t(miRNA$Data)[,w.miRNA[[w[j]]]])
                        typevar = c(typevar,paste("miRNA_",miRNA$Des[w.miRNA[[w[j]]]],sep=""))
                    }
                }
            }
            # include Methylated gene expression and non-Methylated gene expression
            if (!is.null(mRNA$Data)|!is.null(Methylation$Data)) {
                for (j in 1:length(gg)) {
                    g = gg[j]
                    exp.val = mRNA$Data[which(mRNA$Des[,1] == g),]
                    cg = anno.methyl[which(anno.methyl[,1] == g),2]
                    if (!is.na(cg)) {
                        fit = lm(exp.val~Methylation$Data[which(Methylation$Des[,1] == cg),])
                        ME = NME = rep(NA,n)
                        ME[match(names(fit$fitted.values),names(exp.val))] = fit$fitted.values
                        NME[match(names(fit$residuals),names(exp.val))] = fit$residuals
                        X = cbind(X,ME,NME)
                        typevar = c(typevar,paste("ME_",g,sep=""),paste("NME_",g,sep=""))
                    }
                }
            }
        }else {
            others = setdiff(1:p,c(i,pa))
            if (length(others)>0) {
                X = cbind(X,RPPAdat[,others])
                typevar = c(typevar,paste("RPPA_",nodes[others],sep=""))
            }
        }
        colnames(X) = typevar
        Xlist[[i]] = scale(X) # missing values included
    }
    return(list(ylist=ylist,Xlist=Xlist))
}

