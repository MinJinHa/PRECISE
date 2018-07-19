grepVec <- function(x,y,exact=F,Union=T,ignore.case=F) {
    # ----------------------------------------
    # grep a vector x from a vector y
    # Input :
    #     - x : a vector
    #     - y : a vector
    #     - exact : a logical value. If true, find exact matches. Otherwise, grep the pattern
    #     - union : a logical value. If true, find union of matches for y. Otherwise, find matches for each element of y.
    #     - ignore.case: when exact=F and ignore.case=T, case is ignored. Otherwise, case sensitive matching.
    # Output: position of x in y
    # ----------------------------------------
    
    if (Union) {
        if (!exact) {
            unique(unlist(sapply(1:length(x),function(i) grep(x[i],y,ignore.case=ignore.case),simplify=F)))
        }else {unique(unlist(sapply(1:length(x),function(i) which(y==x[i]),simplify=F)))}
    }else {
        if (!exact){
            sapply(1:length(x),function(i) grep(x[i],y,ignore.case=ignore.case),simplify=F)
        }else {sapply(1:length(x),function(i) which(y==x[i]),simplify=F)}
    }
}
