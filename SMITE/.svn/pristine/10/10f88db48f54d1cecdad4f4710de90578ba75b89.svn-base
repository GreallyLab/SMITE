setMethod (
    "print",
    signature="PvalueAnnotation",
    function(x,...){
        cat("*** Class P-value Annotation, method Print *** \n")
        if(nrow(exprs(slot(x,"expression"))) == 0){
            cat("* Expression data :"); print("0 entries")
        }
        else {
            cat("* Expression data :");
            print(Biobase::pData(slot(x,"expression")));
            print(paste("Expression has",
                        nrow(Biobase::pData(slot(x, "expression"))),
                        "entries"))
        }
        readline()
        if(length(unlist(slot(x,"modifications"))) == 0){
            cat("* Modification data :"); print("0 entries")
        }
        else{
            for(i in unique(unlist(slot(x, "modifications"))$type)){
                cat(paste("*", i, "data :"));
                print(subset(unlist(slot(x, "modifications")),
                             unlist(slot(x, "modifications"))$type==i));
                print(paste(i,"has", length(
                    subset(unlist(slot(x, "modifications")),
                           unlist(slot(x, "modifications"))$type==i)), "entries"))
            }
        }
        readline()
        if(nrow(slot(slot(x,"score_data"), "pval_data")) == 0){
            cat("* Score data :"); print("0 entries")
        }
        else {
            cat("* Score data :"); print(slot(slot(x,"score_data"),"scores"))
            print(paste("Score data has",
                        nrow(slot(slot(x, "score_data"), "scores")),
                        "entries"))
        }
        readline()
        cat("******* End Print (PvalueAnnotation) ******* \n")
    }
)

