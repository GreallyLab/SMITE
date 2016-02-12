
setMethod(
    "show",
    signature="PvalueAnnotation",
    function(object){
        cat("*** Class P-Value Annotation, method Show *** \n")
        lengthexpdata <- min(10,
                             nrow(Biobase::pData(slot(object, "expression"))))
        nrowShow <- min(10,nrow(unlist(slot(object,"modifications"))))
        if(lengthexpdata != 0){
            cat(" First 10 Expression data are shown: \n")
            print(Biobase::pData(slot(object,"expression"))[1:lengthexpdata,],
                        quote=FALSE)
        }
        else{}
        if(length(unlist(slot(object,"modifications"))) != 0){
            for(i in unique(unlist(slot(object,"modifications"))$type)){
                cat(paste( " First 10 rows of", i, " data are shown: \n"))
                print(subset(unlist(slot(object,"modifications")),
                             unlist(slot(object,"modifications"))$type ==
                                 i)[1:nrowShow,])
            }

        }
        else{}

        if(nrow(slot(slot(object,"score_data"),"scores")) != 0){
            cat(paste( " First 10 rows of scores  are shown: \n"))
            print(slot(slot(object,"score_data"),"scores")[1:nrowShow,])
        }
        else{}
        cat("******* End Show (pvalue.object) ******* \n")
    }
)

setMethod(
    "show",
    signature="PvalueObject",
    function(object){
        cat("*** Class P-Value Object, method Show *** \n")

        if(nrow(slot(object, "pval_data")) != 0){
            cat(" First 10 rows of pval_data data are shown: \n")
            temp <- slot(object, "pval_data")[1:10,]
            rownames(temp) <- slot(object, "genes")[1:10]
            print(temp)
        }
        else{}

        if(nrow(slot(object, "effect_data"))!=0){
            cat(" First 10 rows of effect_data data are shown: \n")
            temp <- slot(object, "effect_data")[1:10,]
            rownames(temp) <- slot(object, "genes")[1:10]
            print(temp)
        }
        else{}

        if(nrow(slot(object, "signs_index")) !=0 ){
            cat("Signs Index: \n")
            print(slot(object, "signs_index"))
        }
        else{}
        if(nrow(slot(object, "scores")) !=0 ){
            cat(" First 10 Scores are shown: \n")
            temp <- slot(object, "scores")[1:10,]
            names(temp) <- slot(object, "genes")[1:10]
            print(temp)
        }
        else{}

        if(length(slot(object,"module_output")) != 0){
            if(slot(object,"module_output")$moduleType == "spinglass"){
                cat("Spinglass data has been loaded \n")
            }
            if(slot(object,"module_output")$moduleType == "BioNet"){
                cat("BioNet data has been loaded \n")
            }
            if(length(slot(object,"module_output")$goseqOut) > 0){
                cat("GoSeq analysis has been loaded \n")
            }
        }
        cat("******* End Show (pvalue.object) ******* \n")
    }
)
