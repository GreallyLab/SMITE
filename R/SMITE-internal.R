
setMethod(
    f="annotationOutput", 
    signature="PvalueAnnotation", 
    definition=function(pvalue_annotation){
        
<<<<<<< HEAD
        temp1 <- slot(slot(pvalue_annotation, "modifications"), "metadata")$m_summary
        temp1 <- cbind(rownames(temp1), temp1)
        colnames(temp1)[1] <- "genes"
        if(nrow(pData(pvalue_annotation@expression)) > 0){
            
            temp2 <- pData(pvalue_annotation@expression)
            temp2 <- cbind(rownames(temp2), temp2)
            temp1 <- merge(temp1, temp2, by=1, all=TRUE)
        }
        temp1
=======
        metadata_sum <- slot(slot(pvalue_annotation), "modifications"), "metadata")$m_summary
        metadata_sum <- cbind(rownames(metadata_sum), metadata_sum)
        colnames(metadata_sum)[1] <- "genes"
        if(nrow(pData(pvalue_annotation)@expression)) > 0){
            
            exprs <- pData(pvalue_annotation)@expression)
            exprs <- cbind(rownames(exprs), exprs)
            metadata_sum <- merge(metadata_sum, exprs, by=1, all=TRUE)
        }
        return(metadata_sum)
>>>>>>> ce95657b8041fa66d5bc5d2ebd3e9d315c4e9f27
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
