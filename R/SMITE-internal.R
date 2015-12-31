setMethod(
    f="annotationOutput",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation){
        metadata_sum <- slot(slot(pvalue_annotation, "modifications"), "metadata")$m_summary
        metadata_sum <- cbind(rownames(metadata_sum), metadata_sum)
        colnames(metadata_sum)[1] <- "genes"
        if(nrow(Biobase::pData(pvalue_annotation@expression)) > 0){

            exprs <- Biobase::pData(pvalue_annotation@expression)
            exprs <- cbind(rownames(exprs), exprs)
            metadata_sum <- merge(metadata_sum, exprs, by=1, all=TRUE)
        }
        metadata_sum
    }
)

setMethod(
    f="returnPvalueCol",
    signature="PvalueObject",
    definition=function(pval_object, col_name){

        return(slot(pval_object,
                    "pval_data")[, grep(col_name,
                                        colnames(slot(pval_object,
                                                      "pval_data")))])
    }
)
