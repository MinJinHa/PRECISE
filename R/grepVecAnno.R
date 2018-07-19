grepVecAnno <- function(x,y,anno) {
    # ----------------------------------------
    # grep a vector x from a vector y using annotation file
    # Input :
    #     - x : a vector
    #     - y : a vector
    #     - anno : an annotation matrix which links x and y. The first and second column includes elements for x and y, respectively.
    # Output: a list which includes positions of x in y
    # Note : anno file is not necessarily one to one, multiple positions are possible for an element of x
    # ----------------------------------------
    stopifnot(is.matrix(anno),length(x)>0,length(y)>0)
    
    w.anno = grepVec(x,anno[,1],exact=T,Union=F)
    sapply(1:length(x),function(i) grepVec(anno[w.anno[[i]],2],y,exact=T,Union=T),simplify=F)
}
