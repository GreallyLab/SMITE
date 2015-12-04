
setMethod(
    f="annotationOutput", 
    signature="PvalueAnnotation", 
    definition=function(pvalue_annotation)){
        
        temp1 <- slot(slot(pvalue_annotation), "modifications"), "metadata")$m_summary
        temp1 <- cbind(rownames(temp1), temp1)
        colnames(temp1)[1] <- "genes"
        if(nrow(pData(pvalue_annotation)@expression)) > 0){
            
            temp2 <- pData(pvalue_annotation)@expression)
            temp2 <- cbind(rownames(temp2), temp2)
            temp1 <- merge(temp1, temp2, by=1, all=TRUE)
        }
        return(temp1)
    }
)
	
setMethod(
    f="returnPvalueCol", 
    signature="PvalueObject", 
    definition=function(pval_object, c_name){        
        
        return(slot(pval_object, 
                    "pval_data")[, grep(c_name, 
                                        colnames(slot(pval_object, 
                                                      "pval_data")))])
    }
)
