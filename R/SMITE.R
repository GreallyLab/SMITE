##SMITE_1_0_2 10/9/2015


##internal function to perform a Stouffer's test for p value combination

setMethod(
    f="stoufferTest", 
    signature="vector", 
    definition=function(p, w)
    {
        if(is.null(w)){
            w <- rep(1, length(p))/length(p)
        } 
        
        Zi <- qnorm(1-p/2) 
        Z  <- sum(w*Zi)/sqrt(sum(w^2))
        p_val <- (1-pnorm(Z))*2
        p_val <- replace(p_val,p_val < .0000000000000001, .0000000000000001)
        p_val <- replace(p_val, p_val > .9999999999999999, .9999999999999999)
        p_val
    }
)

##SMITE Functions

setMethod(
    f="makePvalueAnnotation", 
    signature="data.frame", 
    definition=function(data, other_data=NULL, other_tss_distance=10000, 
                        promoter_upstream_distance=1000, promoter_downstream_distance=1000, 
                        strand_col=NULL, gene_name_col=NULL)
    {
        ##if the strand coloumn was not specified auto-detect
        if(is.null(strand_col)){strand_col<-which(data[1, ]%in%c("+", "-"))[1]}
        
        data_grange <- GRanges(seqnames=data[, 1], 
                             ranges=IRanges(start=data[, 2], end=data[, 3]), 
                             genenames=data[, gene_name_col], 
                             strand=data[, strand_col])
        
        if(any(duplicated(mcols(data_grange)$genenames)))
        {
            message("Genes are duplicated.  Removing duplicates")
            data_grange <-
                subset(data_grange,!duplicated(mcols(data_grange)$genenames))
        }
        
        mcols(data_grange)$feature <- "original"
        
        if(!is.null(other_data)){
            if(is.null(names(other_data))){
                otherdata_names <- paste("other", 1:length(other_data), sep="")
                names(other_data) <- otherdata_names
            }
            else{
                otherdata_names <- names(other_data)
            }
            
        }
        if(!is.null(other_data)){
            if(length(other_tss_distance) != length(other_data)){
                other_tss_distance <- rep(other_tss_distance[1], 
                                          length(other_data))
            }
        if(is.null(names(other_tss_distance))){
            names(other_tss_distance) <- otherdata_names}
        }
        
        tss <- shift(flank(data_grange, width=2), 1)
        mcols(tss)$feature <- "tss"
        
        if(!is.null(other_data)){
            other <- do.call(c, sapply(1:length(other_data), function(i){
                tempother<-c(GRanges(seqnames=other_data[[i]][, 1], 
                                     ranges=IRanges(start=other_data[[i]][, 2], 
                                                    end=other_data[[i]][, 3])))
                temp_other<-unique(temp_other)
                
                
                suppressWarnings(
                    overlap <- findOverlaps(flank(data_grange, 
                                                 other_tss_distance[otherdata_names[i]], 
                                                 start=TRUE), temp_other)
                )
                temp_other <- temp_other[as.numeric(subjectHits(overlap))]
                mcols(temp_other)$genenames <-
                    mcols(data_grange)$genenames[queryHits(overlap)]
                mcols(temp_other)$feature <- otherdata_names[i]
                temp_other
            })
            )
        }
        promoters_downstream <- flank(data_grange, -promoter_downstream_distance, 
                                    start=TRUE)
        promoters_upstream <- flank(data_grange, promoter_upstream_distance, start=TRUE)
        end(promoters_upstream) <- end(promoters_upstream)+1
        
        promoters <- punion(promoters_upstream, promoters_downstream)
        mcols(promoters)$genenames <- mcols(data_grange)$genenames
        mcols(promoters)$feature <- "promoter"
        
        body <- psetdiff(data_grange, promoters_downstream)
        mcols(body)$genenames <- mcols(data_grange)$genenames
        mcols(body)$feature <- "body"
        
        
        if(!is.null(other_data)){
            suppressWarnings(combined_data <- 
                                 c(data_grange, promoters, body, other, tss))
        }
        else{
            suppressWarnings(combined_data <- 
                                 c(data_grange, promoters, body, tss))
        }
        combined_data <- split(combined_data, mcols(combined_data)$genenames)
        slot(combined_data, "metadata")$params <-
            c(
                gene_name_col=gene_name_col, 
                promoter_upstream_distance_tss=promoter_upstream_distance, 
                promoter_downstream_distance_tss=promoter_downstream_distance, 
                other_annotation_distance_tss=other_tss_distance
            )
        slot(combined_data, "metadata")$sizes_summary <- lapply(
            split(unlist(combined_data), unlist(combined_data)$feature), 
            function(i){each_length <- width(i); 
                        c(summary(each_length), sds=sd(each_length))
            })
        new("PvalueAnnotation", annotation=combined_data)
    })


##Convert ids between types refseq, ensembleprot and uniprot to gene symbol

setMethod(
    f="convertGeneIds", 
    signature(a="character", a_type="character", b_type="character"), 
    definition=function(a, a_type, b_type, delim=NULL, verbose=FALSE)
    {
        if(any(duplicated(a))){stop(
            "Cannot convert duplicated ids. Please remove duplicates.")
        }
        
        if(!is.null(delim)){
            a <- do.call(rbind, strsplit(a, delim))[, 2]
        } 
        
        genes_old <- unique(as.character(a))
        a <- cbind(a, 1:length(a))
        
        if(a_type == "refseq"){
            genes_old <- subset(genes_old, genes_old %in% 
                                    (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egREFSEQ2EG)))
            
            if(b_type == "symbol"){
                
                eg <- unlist(AnnotationDbi::mget(genes_old,
                                                 org.Hs.eg.db::org.Hs.egREFSEQ2EG))
                symbol <- unlist(AnnotationDbi::mget(eg, 
                                                     org.Hs.eg.db::org.Hs.egSYMBOL))
                out <- cbind(names(eg), symbol)
            }
        }        
        
        else 
            if(a_type == "ensembleprot"){
                genes_old <- subset(genes_old, genes_old %in% 
                                    (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egENSEMBLPROT2EG)))
            
                if(b_type == "symbol"){
                eg <- unlist(AnnotationDbi::mget(genes_old, 
                                                 org.Hs.eg.db::org.Hs.egREFSEQ2EG))
                symbol <- unlist(AnnotationDbi::mget(eg, 
                                                     org.Hs.eg.db::org.Hs.egSYMBOL))
                out <- cbind(names(eg), symbol)
            }
        }
        else 
            if(a_type == "uniprot"){
                genes_old <- subset(genes_old, genes_old %in%
                                        (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egUNIPROT)))
                
                if(b_type == "symbol"){
                    eg <- unlist(AnnotationDbi::mget(genes_old, 
                                                     org.Hs.eg.db::org.Hs.egREFSEQ2EG))
                    symbol <- unlist(AnnotationDbi::mget(eg, 
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(eg), symbol)
                }
            }
        else 
            if(a_type == "ensemble"){
                genes_old <- subset(genes_old, genes_old %in%
                                        (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egENSEMBL2EG)))
                if(b_type == "symbol"){
                    eg <- unlist(AnnotationDbi::mget(genes_old, 
                                                     org.Hs.eg.db::org.Hs.egENSEMBL2EG))
                    symbol <- unlist(AnnotationDbi::mget(eg, 
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(eg), symbol)
                }
            }
        else
            if(a_type == "entrez"){
                if(b_type == "symbol"){
                    
                    symbol <- unlist(AnnotationDbi::mget(genes_old, 
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(symbol), symbol)
                }
            }
        else
            if(a_type == "symbol"){
                genes_old <- subset(genes_old, genes_old %in%
                                        (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egALIAS2EG)))
                if(b_type == "entrez"){
                    eg <- unlist(AnnotationDbi::mget(genes_old, 
                                                     org.Hs.eg.db::org.Hs.egALIAS2EG))
                    out <- cbind(names(eg), eg)
                }
                
            }
        
        out <- merge(a, out, by=1, all.x=TRUE)
        out <- out[order(as.numeric(out[, 2])), ]
        out <- subset(out,!duplicated(out[,1]))
        
        out[, 3]
    })

setMethod(
    f="annotateExpression", 
    signature="PvalueAnnotation", 
    definition=function(annotation, expdata, effect_col=NULL, pval_col=NULL){
        if(is.null(effect_col)){
            effect_col=grep("effect|odds|coeff|B", tolower(colnames(expdata)))
            if(length(effect_col) != 1){
                stop("Cannot determine effect column. Please specify with arg:effect_col")
            }
        }
        
        if(any(!c(-1,1) %in% unique(sign(expdata[,effect_col])))){
            message("WARNING: Effects should provide a direction, but these effects
            are all in the same direction.")
        }
        
        if(is.null(pval_col)){
            
            pval_col=grep("pval|p.val|p_val|sig", tolower(colnames(expdata)))
            if(length(effect_col) != 1){
                stop("Cannot determine p.value column. Please specify with arg:effect_col")}
        }
        if(any(expdata[,pval_col] < 0, expdata[,pval_col] > 1)){
            stop("P-values must be between 0 and 1")
        }
        
        if(any(expdata[, pval_col] < .0000001)){
            expdata[which(expdata[, pval_col] < .0000001), pval_col] <- .0000001
        }
        if(any(expdata[, pval_col] > .9999999)){
            expdata[which(expdata[, pval_col] > .9999999), pval_col] <- .9999999
        }
        
        slot(annotation, "expression") =
            ExpressionSet(as.matrix(expdata), featureNames=rownames(expdata))
        phenoData(slot(annotation, "expression")) =
            new("AnnotatedDataFrame", 
                data=as.data.frame(exprs(slot(annotation, 
                                              "expression"))[, 
                                                            c(effect_col, pval_col)]))
        varLabels(slot(annotation, "expression")) <- c("expression_effect", 
                                                      "expression_pvalue")
        
        annotation
        })



setMethod(
    f="annotateModification", 
    signature="PvalueAnnotation",
    definition=function(annotation, modData, weight_by=NULL, 
                        weight_by_method="Stouffer", modInclude=NULL, 
                        modCorr=TRUE, modType="methylation", verbose=FALSE){
        if(modType %in% names(slot(annotation,"modifications")@metadata$elements))
        {
            stop("Provided data set is already loaded as modType")
        } 
        
        ##no weights provided
        if(missing(weight_by)){
            weight_by <- rep("pval", length(unique(mcols(
                unlist(
                    slot(annotation, "annotation")))$feature)[!unique(mcols(
                        unlist(slot(annotation, "annotation")))$feature) %in% 
                            c("original", "tss")]))
        }
        
        ##no modInclude or weight names
        if(is.null(names(weight_by))){
            if(is.null(modInclude)){
                modInclude=unique(mcols(unlist(
                    slot(annotation, "annotation")))$feature)[!unique(mcols(
                        unlist(slot(annotation, "annotation")))$feature) %in%
                            c("original", "tss")]
            }
            names(weight_by) <- modInclude
        }
        
        ##weight names were provided
        if(!is.null(names(weight_by))){   
            modInclude <- names(weight_by)
            
            if(!all(modInclude%in%unique(mcols(unlist(
                slot(annotation, "annotation")))$feature))){
                
                stop("Provided weight names must match those in 
                     unique(mcols(unlist(annotation@annotation))$feature)")
            } 
        }
        
        if(any(!c(-1,1) %in% unique(sign(modData[,4])))){
            message("WARNING: Effects should provide a direction, 
                but these effects are all in the same direction.")
        }
        
        if(any(modData[,5]<0,modData[,5]>1)){
            stop("P-values must be between 0 and 1")
        }
        
        
        mod_grange <- GRanges(seqnames=modData[, 1], ranges=IRanges(
            start=modData[, 2], end=modData[, 3]), effect=modData[, 4], 
            pval=modData[, 5], type=modType)
        temp_annotation <- unlist(slot(annotation, "annotation"))
        overlap_sub <- findOverlaps(temp_annotation, mod_grange)
        out <- mod_grange[subjectHits(overlap_sub)]
        mcols(out) <- cbind(mcols(temp_annotation[as.numeric(
            queryHits(overlap_sub))]), mcols(out))
        out <- split(out, mcols(out)$genenames)
        temp_annotation <- split(temp_annotation, mcols(temp_annotation)$genenames)
        
        if(modCorr == TRUE){
            if(verbose == TRUE){ 
                message("Computing correlation matrices")
            }
            mod_grange_corr <- mod_grange[c(precede(mod_grange, mod_grange),
                                          follow(mod_grange, mod_grange))]
            mod_grange_corr$distance <- distance(c(mod_grange, mod_grange), 
                                               mod_grange_corr)
            mod_grange_corr$pval2 <- c(mod_grange, mod_grange)$pval
            cutpoints <- cut2(mod_grange_corr$distance, g=500, onlycuts=TRUE)
            cutpoints[length(cutpoints)] <- 250000000
            mod_grange_corr$cat <- cut(mod_grange_corr$distance, breaks=cutpoints)
            mod_grange_corr <- split(mod_grange_corr, mod_grange_corr$cat)
            mod_grange_corr2 <- lapply(mod_grange_corr, function(j)
            { 
                mean((sapply(1:100, function(i){
                    index <- sample(1:length(j), replace=TRUE);
                    cor(qnorm(1-j$pval[index]), qnorm(1-j$pval2[index]))
                })))
                
            })
            final_corr <- do.call(rbind, mod_grange_corr2)
            final_corr <- data.frame(final_corr, row.names(final_corr))
            final_corr <- rbind(c(.9, paste("(-1, ", cutpoints[1], "]", 
                                          sep="")), final_corr)
            rm(mod_grange_corr)
        }
        mylist <- sapply(modInclude, function(i){
            if(verbose == TRUE){
                message(paste("Combining p-values over:", i))
            }
            temp <- subset(unlist(out), mcols(unlist(out))$feature == i)
            
            
            ref_data <- slot(annotation, "annotation")[names(temp)]
            ref_data <-
                subset(unlist(slot(annotation, "annotation")),
                       mcols(unlist(
                           slot(annotation, "annotation")))$feature == "tss")
            ref_data <- ref_data[mcols(temp)$genenames]
            suppressWarnings(mcols(temp)$distance <- distance(ref_data, 
                                                            temp)+2)
            temp <- split(temp, mcols(temp)$genenames)
            forreturn <- lapply(temp, function(each){
                each_length <- length(each)
                each_effect <- each$effect[order(each$pval)]
                each_pval <- each$pval[order(each$pval)]
                distances <- each$distance[order(each$pval)]
                if(length(each_pval > 1)){
                    if(modCorr == TRUE){
                        corr_mat <- matrix(as.numeric(final_corr[match(cut(
                            as.matrix(dist(start(each)[order(each$pval)])), 
                            breaks=c(-1, cutpoints)), final_corr[, 2]), 1]), 
                            ncol=each_length)
                        diag(corr_mat) < -1
                        chol_d <- try(chol(corr_mat), silent=TRUE)
                        if(is(chol_d, "try-error"))
                        {
                            index<- as.numeric(strsplit(
                                strsplit(chol_d[1], 
                                         "the leading minor of order ")[[1]][2], 
                                " is not positive")[[1]][1])-1
                            chol_d <- chol(corr_mat[1:index, 1:index])
                            each_pval <- each_pval[1:index]
                            each_effect <- each_effect[1:index]
                            distances <- distances[1:index]
                            each_length <- index
                        }
                        
                        each_pval <- 1-pnorm(abs(solve(chol_d) %*% qnorm(
                            1-each_pval/2)))
                        if(any(each_pval == 0)){
                            each_pval[which(each_pval == 0)] <- 0.000000001
                        }
                        each_pval <- each_pval*2
                    } 
                    
                    if(weight_by_method == "Stouffer"){
                        
                        if(weight_by[i] == "distance"){
                            ##mean is weighted by distance
                            out_mean = weighted.mean(each_effect,
                                                    w=(1/log(distances)))
                            ##Stouffer test is weighted by distance
                            out_pval <- stoufferTest(each_pval, w=(1/log(distances)))
                        } 
                        else if(weight_by[i] %in%
                                      c("pval", "p.value", "pvalue", "p_val")){
                            ##mean is weight by pvalue
                            out_mean = weighted.mean(each_effect, w=-log(each_pval))
                            out_pval <- stoufferTest(each_pval, w=NULL)
                        } 
                        else {
                            ##mean is not weighted
                            out_mean = mean(each_effect, na.rm=TRUE)
                            out_pval <- stoufferTest(each_pval, w=NULL)    
                        }
                              
                    } 
                    else if(weight_by_method %in% c("minimum", "Sidak", "sidak")){
                        
                        index <- which(each_pval == min(each_pval))
                        if(length(index) > 1){
                            index <- index[which(
                                abs(each_effect[index]) == max(abs(each_effect[index])))][1]
                        }
                        out_mean = each_effect[index]
                        out_pval <- 1-(1-each_pval[index])^length(each_pval)
                        
                    } 
                    else if(weight_by_method == "binomial"){
                        
                        index <- which(each_pval == min(each_pval))
                        if(length(index) > 1){
                            index <- index[which(abs(each_effect[index]) == max(
                                abs(each_effect[index])))][1]
                        }
                        out_mean = each_effect[index]
                        out_pval <- (1-pbinom(q=length(which(each_pval<0.05)), 
                                            size=each_length, prob=0.05))
                    } else if(weight_by_method%in%
                                  c("Fisher", "fisher", "chisq", "chi")){
                        out_pval <- 1-pchisq(-2*sum(log(each_pval)), each_length*2)
                        out_mean <- mean(sign(each_effect), na.rm=TRUE)
                        
                    }
                }
                #if only one pval
                else{
                    out_mean <- each_effect
                    out_pval <- each_pval
                }
                
                c(out_mean, out_pval, each_length)
            })
            
            do.call(rbind, forreturn)
        })
        
        if(verbose == TRUE){ 
            message("Quantile permuting scores")
        }
        mylist <- lapply(mylist, function(each_feature){		
            cats <- data.frame(cats=as.numeric(cut2(each_feature[, 3], g=100)))
            cats_table <- data.frame(table(cats))
            
            trans_p <- cbind(trans=qnorm(1-each_feature[, 2]/2), 
                           join(cats, cats_table, by="cats"))
            if(any(is.infinite(trans_p[, 1]))){
                trans_p[which(is.infinite(trans_p[, 1])), 1] <- max(
                    subset(trans_p,!is.infinite(trans_p[, 1]))[,1])
            }
            
            num_list <- split(trans_p$trans, trans_p$cats)
            rand_list <- sapply(1:length(num_list), function(i){
                as.matrix(sapply(1:500, function(j){ 
                    sample(num_list[[as.numeric(i)]], replace=TRUE)
                }))
            })
            new_pval <- apply(trans_p, 1, function(i){ 
                length(which(
                    rand_list[[as.numeric(i[2])]] > as.numeric(i[1])))/(500*as.numeric(
                        i[3]))
            })
            new_pval <- replace(new_pval, new_pval == 0, 
                                min(subset(new_pval, new_pval != 0), na.rm=T))
            
            each_feature[,2] <- new_pval 
            each_feature <- as.data.frame(each_feature)
            each_feature
        })
        
        if(verbose == TRUE){ 
            message("Scores have been adjusted")
        }
        
        newmods <- c(unlist(slot(annotation, "modifications")), unlist(out))
        names(newmods) <- NULL
        newmods <- split(newmods, mcols(newmods)$genenames)
        
        final <- suppressWarnings(as.data.frame(c(list(names=names(out)), 
                                                lapply(mylist, function(x){
                                                    x[match(names(out), rownames(x)), 1:2]
                                                    }))))
        rownames(final) <- final[, 1]
        final <- final[, -1]
        colnames(final) <- paste(modType, apply(expand.grid(c("effect", "pvalue"),
                                                            modInclude),
                                                1, function(i){
                                                        paste(i[2], i[1], sep="_")
                                                    }), sep="_")
        
        newmetadata <- slot(slot(annotation, "modifications"), "metadata")
        if(is.null(newmetadata$m_summary)){
            newmetadata$m_summary <- final
        } 
        else{
            newmetadata$m_summary <- merge(newmetadata$m_summary, final, by=0, 
                                         all=TRUE)
            rownames(newmetadata$m_summary) <- newmetadata$m_summary[, 1]
            newmetadata$m_summary<-newmetadata$m_summary[, -1]
        }
        newmetadata[["elements"]][[modType]]$weight_by <- weight_by
        newmetadata[["elements"]][[modType]]$weight_by_method <- weight_by_method
        newmetadata$elementnames <- c(newmetadata$elementnames, paste(modType, 
                                                                    modInclude, sep="_"))
        slot(annotation, "modifications") <- newmods
        slot(slot(annotation, "modifications"), "metadata") <- newmetadata
        
        annotation
        })


setMethod(
    f="removeModification", 
    signature="PvalueAnnotation", 
    definition=function(annotation, modType="methylation"){
        
        if(!modType%in%names(slot(annotation,
                                  "modifications")@metadata$elements)){
            stop("Provided modType is not in the annotation")
        }
        
        temp_meta <- slot(slot(annotation,"modifications"),"metadata")
        temp <- unlist(slot(annotation,"modifications"))
        names(temp) <- NULL
        temp <- subset(temp,!(temp$type == modType))
        slot(temp,"metadata") <- temp_meta
        slot(temp,"metadata")$m_summary <-
            slot(temp,"metadata")$m_summary[, 
                                            -grep(modType, colnames(slot(temp,"metadata")$m_summary))]
        
        if(ncol(slot(temp,"metadata")$m_summary) == 0){
            slot(temp,"metadata")$m_summary <- NULL
        }
        
        slot(temp,"metadata")$elements[which(
            names(slot(temp,"metadata")$elements) == modType)] <- NULL
        
        
        slot(temp,"metadata")$elementnames <-      
            slot(temp,"metadata")$elementnames[
                -which(do.call(rbind, lapply(strsplit(
                    slot(temp,"metadata")$elementnames, "_"),
                    function(i)i[1]))[, 1] %in% modType)]
        temp_meta <- slot(temp,"metadata")
        slot(annotation,"modifications") <- split(temp, temp$genenames)
        slot(slot(annotation,"modifications"),"metadata") <- temp_meta
        annotation
    }
)

setMethod(
    f="makePvalueObject", 
    signature="PvalueAnnotation", 
    definition=function(object, effect_directions=NULL) {
        
        
        if(is.null(effect_directions)){
            effect_directions <- rep("bidirectional", 
                                     length(slot(
                                                slot(object, "modifications"),
                                                "metadata")$elementnames))
        }
        
        if(is.null(names(effect_directions))){
            names(effect_directions)<-slot(slot(object, "modifications"), 
                                           "metadata")$elementnames
        }
        
        if(any(!names(effect_directions) %in% slot(slot(object, "modifications"), 
                                                 "metadata")$elementnames)){
            stop("Effect name is invalid")
        }
        
        if(any(!(effect_directions %in% c("decrease", "increase", "bidirectional")
        ))){
            stop("Effect argument is invalid.")
        }
        
        
        exp_ind <- ifelse(nrow(pData(object@expression)) > 0, 1, 0)
        
        totalfactor <- length(effect_directions)
        
        signsindex <- merge(cbind(c("increase", "decrease", "bidirectional"),
                                c(1, -1, 2)), 
                          cbind(effect_directions, 1:totalfactor), by=1)
        signsindex <- signsindex[order(signsindex[, 3]), ]
        rownames(signsindex) <- NULL
        colnames(signsindex) <- c("expression_relationship", "B_coeff","name")
        
        temp1 <- annotationOutput(object)
        genenames <- temp1[, 1]
        data <- temp1[, -1]
        
        if(exp_ind == 1){
            data <- data[, c(as.numeric(sapply(names(effect_directions), 
                                             function(i){
                                                 grep(i, colnames(data))
                                             })),
                            grep("exp", colnames(data)))]
            
        }
        else {
            data <- data[, as.numeric(sapply(names(effect_directions), 
                                             function(i){
                                                 grep(i, colnames(data))
                                             }))]
        }
        signsindex[, 3] <- names(effect_directions)
        slot(object, "score_data") <- new(Class="PvalueObject", 
                                        pval_data=data[, grep("pval",
                                                              colnames(data))], 
                                        effect_data=data[, grep("effect", 
                                                                colnames(data))], 
                                        genes=genenames, 
                                        signsindex=signsindex)
        object
    }
)




setMethod(
    f="normalizePval", 
    signature="PvalueAnnotation", 
    definition=function(object, trans, ref="expression_pvalue", 
                        method="rescale"){
        temp_pval_data <- slot(slot(object,"score_data"),"pval_data")
        names_temp_pval_data <- colnames(temp_pval_data)
        
        if(nrow(temp_pval_data==0){stop(
            "Run makePvalueObject function first.")
        }
        if(!any(grepl(ref,names_temp_pval_data))){
               stop("paste(Reference is not one of the available:", 
                    names_temp_pval_data)
        }
        if(length(grep(ref,names_temp_pval_data)) > 1){
               stop("Reference was not specific enough.")
        }
           
        ref_index <- grep(ref,names_temp_pval_data)
        p_ref <- temp_pval_object[[ref_index]]
        logit_ref <- log(p_ref/(1-p_ref))
        
        if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
            dev.new(height=7,width=14)
        }
        
        par(mfrow=c(1, 2))
        plotDensityPval(object, ref=ref)
           
        if(method %in% c("Box-Cox", "box-cox", "boxcox", "Boxcox")){
            if(missing(trans)){
                message(paste("Auto-detecting best transformation"))
                expos <- c()
                for(x in subset(names_temp_pval_data,
                                   !names_temp_pval_data %in% ref)){
                       
                    p_temp<-temp_pval_object[[x]]
                    if(!all(is.na(p_temp))){
                        logit_temp <- log(p_temp/(1-p_temp))
                        each <- t(sapply(c(seq(.05, .95, .05), rev(1/seq(.05,
                                                                          .95, .05))),
                                          function(i){
                                              c(i, wilcox.test(logit_ref,
                                                               as.numeric(logit_temp)*i)$p_value)
                                              }))
                           each <- each[which(each[, 2] == max(each[, 2]))[1], 1]
                           p_temp <-
                               (exp(logit_temp)^(each))/(1+exp(logit_temp)^(each))
                    } 
                    else {
                        each<-1
                    }
                       expos <- c(expos, each)
                       names(expos)[length(expos)] <- x
                       temp_pval_object[[x]] <- p_temp
                       
                }
            }
            
            else {
                if(length(trans) != length(grep(paste(slot(slot(object,
                                                                   "score_data"), "signsindex")[, 3], collapse="|"), 
                                                   colnames(temp_pval_object)))){
                       stop("Length of p and transformations must equal!")
                }
                else {
                    expos <- trans 
                    names(expos) <- subset(names_temp_pval_data,!names_temp_pval_data  %in%ref)
                }
                for(x in names(expos)){
                    p_temp <- temp_pval_object[[x]]
                    if(!all(is.na(p_temp))){
                           p_temp <- p_temp^expos[x]
                    }
                    temp_pval_object[[x]] <- p_temp
                }
            }
        }
        else if(method%in%c("Rescale", "rescale")){
            for(i in slot(slot(object,"score_data"),"signsindex")[, 3]){
                if(!all(is.na(slot(
                    slot(object, "score_data"), "pval_data")[[grep(i,
                                                                   names_temp_pval_data)]])
                   )){
                    p_temp <- temp_pval_object[[grep(i, names_temp_pval_data)]]
                    logit_temp <- log(p_temp/(1-p_temp))
                    logit_temp <- rescale(logit_temp, to=range(logit_ref, 
                                                                    na.rm=TRUE))
                    temp_pval_object[[grep(i, 
                                           names_temp_pval_data
                                           )]] <- exp(logit_temp)/(1+exp(logit_temp))
                }
            }
        }
           
        slot(slot(object,"score_data"),"pval_data") <- temp_pval_data				
        plotDensityPval(object, ref=ref)   
        object
    }
)

setMethod(
    f="SMITEscorePval", 
    signature="PvalueAnnotation", 
    definition=function(object, weights){ 
        
        totalfactor <- ncol(slot(slot(object, "score_data"), "pval_data"))
        if(missing(weights)){
            weights <- rep((1/(totalfactor)), totalfactor)
            names(weights) = colnames(slot(slot(object, "score_data"), 
                                         "pval_data"))
        } 
        else {
            if(length(weights) != totalfactor){
                stop("ERROR: Number of factors(with expression data)
                     and weights must equal!")
            }
            else { 
                if(is.null(names(weights))){ 
                    names(weights)=colnames(slot(slot(object,
                                                      "score_data"), "pval_data"))
                }
                if(!all(sort(names(weights)) == sort(do.call(rbind, 
                                                            strsplit(colnames(slot(slot(object, "score_data"), 
                                                                                   "pval_data")), "_pvalue")))))
                {
                    stop(paste("Weight names must match the following:", 
                               paste(do.call(rbind, strsplit(colnames(slot(slot(object, 
                                                                                "score_data"), "pval_data")), "_pvalue")), collapse=", ")))
                }
            }
        }
        
        weights <- weights[match(do.call(rbind, strsplit(colnames(slot(
            slot(object, "score_data"), "pval_data")), "_pvalue")), names(weights))]
        
        message("The following weights are being applied")
        print(weights)
        
        scoringdata <- slot(slot(object, "score_data"), "pval_data")*sign(slot(
            slot(object, "score_data"), "effect_data"))
        
        weight_names_temp <- names(weights)
        
        if(any(grepl("exp", names(weights)))){
            weight_names_temp <- weight_names_temp[-grep("exp", names(weights))]
        }
        unidirectional <- as.numeric(slot(slot(object, "score_data"), 
                                        "signsindex")[match(slot(slot(object, "score_data"),
                                                                 "signsindex")[, 3], weight_names_temp), 2])
        unidirectional[which(unidirectional == 2)] <- 0
        bidirectional <- as.numeric(slot(slot(object, "score_data"), "signsindex")
                                  [match(slot(slot(object, "score_data"), "signsindex")[, 3], 
                                         weight_names_temp), 2])
        bidirectional[which(bidirectional != 2)] <- 0
        bidirectional[which(bidirectional == 2)] <- 1
        slot(slot(object, "score_data"), "scoring_vector") <- weights
        
        scoringdata <- qnorm(1-as.matrix(abs(scoringdata))/2)*sign(scoringdata)
        scoresout <- apply(scoringdata, 1, function(each){
            if(any(!is.na(each)))
            {
                if(any(grepl("exp", names(each)))){
                    
                    forreturn <- (sum(
                        abs(sum(c(as.numeric(each[[grep("exp", names(each))]]), 
                                  as.numeric(each[-grep("exp", names(each))]))*
                                    c(1, unidirectional)*weights[c(grep("exp",
                                                                        names(each)), which(!grepl("exp", 
                                                                                                   names(each))))], na.rm=TRUE)), 
                        (abs(as.numeric(each[-grep("exp",
                                                   names(each))])*bidirectional)*
                             weights[-grep("exp", names(each))]), 
                        na.rm=TRUE)/sum(weights^2)^.5)
                    
                } 
                else {
                    
                    forreturn <- (
                        sum(abs(sum(as.numeric(each)*unidirectional*weights,
                                    na.rm=TRUE)),(abs(as.numeric(each)*bidirectional)*weights),
                            na.rm=TRUE)/sum(weights^2)^.5)
                }
            } 
            else {
                forreturn <- (NA)
            }
            forreturn
        })
        
        scoresout <- (1-pnorm(as.numeric(scoresout)))*2
        if(any(scoresout == 0, na.rm=TRUE)){
            scoresout[which(scoresout == 0)] <- min(subset(scoresout, !(
                scoresout == 0)), na.rm=TRUE)
        }
        scoresout <- (-2*log(scoresout))
        rand_mat <- as.matrix(sapply(1:100, function(j){
            sample(scoresout, replace=TRUE)
            }))
        new_pval <- sapply(scoresout , function(i){ 
            length(which(rand_mat > as.numeric(i)))/(100*length(scoresout))
        })
        new_pval <- replace(new_pval, new_pval==0, 
                            min(subset(new_pval, new_pval!=0), na.rm=T))
        new_pval <- (-2)*log(new_pval[, 1])
        
        
        
        slot(slot(object, "score_data"), "scores") <- 
            data.frame(scores=as.numeric(new_pval), 
                       row.names=as.character(slot(slot(object,"score_data"), 
                                                   "genes")))
        object
    }
) 

setMethod(
    f="SMITErunSpinglass", 
    signature="PvalueAnnotation", 
    definition=function(object, network, random_alpha = 0.05, gam = 0.5, 
                        node_alpha = 0.05, maxsize = 500, minsize = 8,
                        niter = 1000, simplify=TRUE)
    {        
        
        if(inherits(network, what="graphNEL")){
            network <- graph_from_graphnel(network)
        }
        
        if(length(slot(slot(object,"score_data"),"module_output")) !=0 ){
            slot(slot(object,"score_data"),"module_output") <- list()
            message("Overwriting existing modules.")
        }
        
        if(simplify == TRUE){
            network <- igraph::simplify(network,remove.multiple=TRUE, 
                                        remove.loops=TRUE)
        }
        genes_in_network <- subset(slot(slot(object, "score_data"), "genes"), 
                                   slot(slot(object, "score_data"), "genes") %in% 
                                       V(network)$name)
        scores_in_network <- SMITEextractScores(object)[genes_in_network]
        
        ##should be FALSE, but just in case check for NAs 
        if(any(is.na(scores_in_network))){
            ##if there are NAs remove them
            scores_in_network <- subset(scores_in_network, !is.na(
                scores_in_network))
            genes_in_network <- names(scores_in_network)     
        }
        
        nodes_with_scores <- intersect(genes_in_network, V(network)$name)
        network <- induced_subgraph(network, nodes_with_scores)
        network_clusters <- clusters(network)
        
        ## choose largest connected nodes ## may want include > minsize		  
        maxclust <- which(network_clusters$csize == 
                              max(network_clusters$csize))[1]
        network <- induced_subgraph(network, 
                                    which(network_clusters$membership == maxclust)) 
        rm(network_clusters)     
        genes_in_network <- intersect(genes_in_network,V(network)$name)
        
        scores_in_network <- scores_in_network[genes_in_network]
        network.adj <-  as_adjacency_matrix(network)
        
        ## order of scores has to match order of adjacency rows 	   
        scores_in_network <- scores_in_network[rownames(network.adj)]
        genes_in_network <- names(scores_in_network)
        stat_scores <- as.numeric(scores_in_network)
        pval_scores <- exp(scores_in_network/(-2))
        
        temp1 <- apply(network.adj, 1, function(v) v*stat_scores)
        W <- (temp1 + t(temp1))
        rm(temp1)
        gc()
        W_vec <- (-2*log(1-pchisq(as.vector(W),4)))
        if(any(is.infinite(W_vec))){
            W_vec[which(is.infinite(W_vec))] <- 
                max(subset(W_vec[!is.infinite(W_vec)))
        }
        W <- matrix(W_vec, nrow=nrow(W))
        rm(W_vec)
        rownames(W) <- genes_in_network
        colnames(W) <- genes_in_network
        gc()
        final_network <- 
            graph_from_adjacency_matrix(W, mode = "undirected", weighted=TRUE)
        V(final_network)$weight <- stat_scores
        network <- final_network
        rm(final_network)
        gc() 
        
        ## For significant genes apply Spinglass algorithm 
        sig_genes <- subset(genes_in_network, pval_scores < node_alpha)
        sig_genes_counter <- 1:length(sig_genes)
        names(sig_genes_counter) <- sig_genes
        spin.glass.out <- lapply(sig_genes, function(j){
            message(paste("Computing modules: Vertex",sig_genes_counter[j],"of",
                          length(sig_genes), "significant genes is", j))
            genes_in_network[
                cluster_spinglass(network, weights = E(network)$weight, 
                                  vertex=j, gamma=gam)$community]
            })
        names(spin.glass.out) <- sig_genes 
        
        ## Select modules with size requirements
        spin.size <- do.call(c, lapply(spin.glass.out, length))
        spin.glass.out <- subset(spin.glass.out, spin.size >= minsize & spin.size 
                                <= maxsize)
        
        Modularity.edges = function(v, network)
        { 
            h <- induced_subgraph(network, v);
            c(sum( E(h)$weight ))
        }
        edge_sum <- do.call(c,lapply(spin.glass.out, function(j) { 
            Modularity.edges(j, network)
        }));
        
        nspin.glass.out <- length(spin.glass.out);
        random_edges <- lapply(1:nspin.glass.out, function(i) {
            j <- spin.glass.out[[i]]
            h <- induced_subgraph(network, j);
            B <- as_adjacency_matrix(h, sparse=FALSE);
            v <- sapply(1:niter, function(k){
                message(paste("Testing significance: module", i, "of", 
                              nspin.glass.out, "Randomization", k, "of", niter))
                atperm = sample(stat_scores, nrow(B) , replace=TRUE)
                temp1 = apply(B, 1, function(v) v*atperm)
                W <- (temp1 + t(temp1));
                W <- apply(W, 2, function(i){ i[which(i>0)]<-
                                                (-2*log(1-pchisq(subset(i,i>0),4))); i})
                sum(W)/2
            })
            v
        })
        names(random_edges) <- names(spin.glass.out);
        
        random_p <- lapply(1:nspin.glass.out, function(k){
            length(which(random_edges[[k]] >
                             edge_sum[k]))/niter})
        names(random_p) <- names(spin.glass.out)
        
        if(length(spin.glass.out[which(do.call(c,random_p) < random_alpha)]) == 0){
            stop("No modules found.  Please adjust the random_alpha and node_alpha 
                 parameters")
        }
        
        ## Assemble output data
        output <- list();
        output[[1]] <- subset(spin.glass.out, do.call(c,random_p) < random_alpha);
        output[[2]] <- pval_scores;
        output[[3]] <- stat_scores;
        output[[4]] <- induced_subgraph(network, as.character(
            unlist(output[[1]]))) ;
        output[[6]] <- "spinglass";
        
        
        names(output) <- c("modules",  "p_value", 
                          "statistic","network","moduleStats","moduleType");
        k <- 1
        while(k < length(output[[1]])){
            m <- k+1
            while(m <= length(output[[1]])){
                if(all(output[[1]][[m]] %in% output[[1]][[k]])){
                    names(output[[1]])[k] <- paste(names(output[[1]][k]),
                                                   names(output[[1]][m]), 
                                                   sep=":") 
                    output[[1]][[m]]<-NULL 
                    m <- m-1
                } 
                m <- m+1
            }
            k <- k+1
        }
        output[[5]] <- 
            lapply(output[[1]], function(i){
                stat.mod <- sum(abs(subset(stat_scores,genes_in_network %in% i)))
                pval.mod <- 1-pchisq(stat.mod,2*length(i))
                c(statMod=stat.mod, pvalMod=pval.mod, lenMod=length(i))
            })
        
        index <- order(do.call(rbind,output[[5]])[,2])
        output[[1]] <- output[[1]][index]
        output[[5]] <- output[[5]][index]

        slot(slot(object, "score_data"), "module_output") <- output
        object    
    }	
)


setMethod(
    f="SMITErunBioNet", 
    signature="PvalueAnnotation", 
    definition=function(object, network, alpha = 0.05)
    {
        if(any(alpha<0, alpha>1)){
            stop("Error: alpha must be between zero and one.")
        }
        if(length(slot(slot(object,"score_data"),"module_output")) !=0 ){
            slot(slot(object,"score_data"),"module_output") <- list()
            message("Overwriting existing modules.")
        }
        
        scores <- SMITEhighScores(object, alpha=alpha)
        pval.v <- exp(scores/(-2))
        g <- subNetwork(names(pval.v), network)
        g <- rmSelfLoops(g)
        if(inherits(g,what="igraph")){
            g <- as_graphnel(g)
        }
        
        heinzOut <- suppressWarnings(runFastHeinz(g, scores))
        
        output <- list();
        output[[1]] <- list(heinzOut@nodes);
        output[[2]] <- pval.v;
        output[[3]] <- scores;
        output[[4]] <- g;
        output[[6]] <- "BioNet";
        
        names(output) <- c("modules",  "p_value", 
                          "statistic","network","moduleStats","moduleType");
        sub_score <- subset(scores,names(scores)%in%output[[1]][[1]])
        stat.mod <- sum(abs(sub_score))
        pval.mod <- 1-pchisq(stat.mod,2*length(output[[1]][[1]]))
        output[[5]] <-
            list(c(statMod=stat.mod, pvalMod=pval.mod, 
                   lenMod=length(output[[1]][[1]]))
                 )
        
        
        names(output[[1]])[1] <- 
            names(sub_score[which(sub_score==max(sub_score))[1]]
                  )
        slot(slot(object,"score_data"),"module_output") <- output
        object
    }    
)

setMethod(
    f="SMITErunGOseq", 
    signature="PvalueAnnotation", 
    definition=function(object, p_thresh=0.05, coverage, type="reactome")
    {
        while(!names(dev.cur()) %in% c("pdf","null device")){
            dev.off()
        }
        names.eid <- names(slot(slot(object, "score_data"),
                              "module_output")$modules)
        eid <- slot(slot(object, "score_data"), "module_output")$modules
        pval <- slot(slot(object, "score_data"), "module_output")$p_value
        stat <- slot(slot(object, "score_data"), "module_output")$statistic
        genes <- names(slot(slot(object, "score_data"), "module_output")$modules)
        goseqdata <- slot(slot(object, "score_data"), "module_output")$modules
        
        sym2eg <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG)
        eg2sym <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
        
        if(type == "reactome"){
            
            eg2reactome <- 
                AnnotationDbi::as.list(reactome.db::reactomeEXTID2PATHID)
            temp <- sapply(names(eg2reactome), function(y){
                eg2sym[[as.character(y)]]})
            temp <- lapply(temp, function(y){if(is.null(y)){y<-NA};y})
            temp2 <- do.call("c", temp)
            sym2network <- eg2reactome
            names(sym2network) <- temp2
            PATHID2NAME <- AnnotationDbi::as.list(reactome.db::reactomePATHID2NAME)
            
        }
        
        if(type == "kegg"){
            
            eg2kegg <- AnnotationDbi::as.list(KEGG.db::KEGGEXTID2PATHID)
            temp <- sapply(names(eg2kegg), function(y){
                eg2sym[[as.character(y)]]})
            temp <- lapply(temp, function(y){if(is.null(y)){y<-NA};y})
            temp2 <- do.call("c", temp)
            sym2network <- eg2kegg
            names(sym2network) <- temp2
            PATHID2NAME <- AnnotationDbi::as.list(KEGG.db::KEGGPATHID2NAME)
            names(PATHID2NAME) <- paste("hsa", names(PATHID2NAME), sep="")
        }
        
        names(eid) <- paste("module_", names(eid), "_network", sep="")
        nl <- unlist(goseqdata)
        nl <- cbind(names(nl), nl)
        rownames(nl) <- NULL
        nl[, 1] <- gsub("[0-9]*$", "", nl[, 1])
        nl <- split(nl[, 1], nl[, 2])
        nl <- lapply(nl, function(i){unique(i)})
         
        if(inherits(coverage, what="data.frame")){
            a <- coverage
            a <- a[, 4:5]
        }
        if(inherits(coverage, what="character")){
            if(coverage == "refseq"){
                hg19.refGene.LENGTH <- NULL
                data(hg19.refGene.LENGTH,package="geneLenDataBase", 
                     envir=environment())
                a <- hg19.refGene.LENGTH[, c(1, 3)]}
            
            if(coverage == "symbol"){
                hg19.geneSymbol.LENGTH <- NULL
                data(hg19.geneSymbol.LENGTH,package="geneLenDataBase",
                     envir=environment())
                a <- hg19.geneSymbol.LENGTH[, c(1, 3)]
            }
        }
        if(any(duplicated(a[, 1]))){
            a <- split(a, a[, 1])
            a <- lapply(a, function(i){
                if(nrow(i) > 1){
                    i <- i[which(i[, 2] == max(i[, 2])), ]
                };
                i
            })
            a <- do.call(rbind, a)
        }
        
        slot(slot(object, "score_data"), "module_output")$goseqOut <- 
            lapply(slot(slot(object, "score_data"), 
                        "module_output")$modules, function(i){
                            b <- cbind(a[, 1], rep(0, nrow(a)))
                            b[which(b[, 1]%in%i), 2] <- 1
                            x <- as.vector(b[, 2])
                            x <- as.numeric(x)
                            names(x) <- b[, 1]
                            pwf <- nullp(x, 'h19', 'knownGene', bias.data=a[, 2])
                            path <- goseq(pwf, "hg19", "knownGene", 
                                          gene2cat=sym2network)
                            path <- cbind(path,(as.character(sapply(
                                path$category, function(i){PATHID2NAME[[i]]}))))
                            colnames(path)[6] <- "cat_name"
                            subset(path, path$over_represented_pvalue < p_thresh)
                        }
                   )
        object
    }
) 




setMethod(
    f="SMITEsearchGOseq", 
    signature="PvalueAnnotation", 
    definition=function(object, searchstring, wholeword=FALSE){
        options(stringsAsFactors=FALSE)
        searchstring<-tolower(searchstring)
        nums <- which(do.call("c", lapply(
            slot(slot(object, "score_data"), "module_output")$goseqOut,
            function(each){
                any(grepl(ifelse(wholeword == FALSE, searchstring, 
                                 paste("\\b", searchstring, "\\b", sep="")), 
                          tolower(each[, ncol(each)])))
                }
            ))
        )
        if(length(nums) > 0){ 
            out <- lapply(nums, function(i){
                pos <- grep(ifelse(wholeword == FALSE, searchstring, 
                                   paste("\\b", searchstring, "\\b", sep="")),
                            tolower(slot(slot(object, "score_data"), 
                                         "module_output")$goseqOut[[i]][, 6]))
                tot <- nrow(slot(slot(object, "score_data"), 
                                 "module_output")$goseqOut[[i]])
                outterm <- as.data.frame(as.character(slot(
                    slot(object, "score_data"), 
                    "module_output")$goseqOut[[i]][pos, 6]))
                cbind(outterm, as.data.frame(pos), 
                      as.data.frame(rep(tot,length(pos))))
            })
            out <- lapply(1:length(nums), function(i){
                old <- out[[i]]
                rownames(old) <- NULL
                old <- cbind(rep(names(nums)[i], nrow(old)), 
                        rep(as.numeric(nums[i]), nrow(old)), old)
                old
            })
            out <- do.call(rbind, out)       
            out[,2] <- sapply(out[,2], function(i){
                i <- as.numeric(i);
                i <- paste(i, "/", round(slot(
                    slot(object, "score_data"), 
                    "module_output")$moduleStats[[i]][2],4));
                i
            })
         
            colnames(out) <- c("epimod_name", "epimod_position_pval",
                               "term", "rank_of_term", "total_terms")
            rownames(out) <- NULL
            out
        }
        else { 
            message("Search term not found.")
        }
    }
)

setMethod(
    f="SMITEextractGOseq", 
    signature="PvalueAnnotation", 
    definition=function(object, which.network=NULL){  
        temp <- slot(slot(object, "score_data"), "module_output")$goseqOut
        if(!is.null(which.network)){temp <- temp[which.network]}
        temp
    }
)






setMethod(
    f="extractModification", 
    signature="PvalueAnnotation", 
    definition=function(annotation, modType=NULL){
        if(!is.null(modType)){
            if(!modType%in%names(slot(annotation,
                                      "modifications")@metadata$elements)){
                stop("Provided modType is not in the annotation")
            } 
        }
        
        temp_meta <- slot(slot(annotation,"modifications"),"metadata")
        temp <- unlist(slot(annotation,"modifications"))
        names(temp) <- NULL
        
        if(is.null(modType)){modType <- unique(temp$type)}
        temp<-subset(temp, temp$type %in% modType)
        slot(temp, "metadata") <- temp_meta
        temp
    }
)

setMethod(
    f="extractExpression", 
    signature="PvalueAnnotation", 
    definition=function(annotation){
        if(is.null(expression)){
            stop("No expression data loaded.")
        } 
        temp_exp <- slot(annotation,"expression")
        pData(temp_exp)
    }
)




setMethod(
    f="extractModSummary", 
    signature="PvalueAnnotation", 
    definition=function(annotation){ 
        temp_meta<-slot(slot(annotation,"modifications"),"metadata")
        temp_meta$m_summary
    }
)

setMethod(
    f="SMITEextractScores", 
    signature="PvalueAnnotation", 
    definition=function(annotation){
        if(nrow(slot(slot(annotation,"score_data"),"scores")) == 0){
            stop("Run SMITEscorePval function first.")
        }
        temp <- slot(slot(annotation,"score_data"),"scores")
        names_temp <- rownames(temp)
        temp <- as.vector(temp[,1])
        names(temp) <- names_temp
        temp
    }
)



setMethod(
    f="SMITEhighScores", 
    signature="PvalueAnnotation", 
    definition=function(annotation,alpha=0.05){
        if(nrow(slot(slot(annotation,"score_data"),"scores")) == 0){
            stop("Run SMITEscorePval function first.")
        }
        if(any(alpha < 0, alpha > 1)){
            stop("alpha must be between 0 and 1")
        }
        
        temp <- slot(slot(annotation,"score_data"),"scores")
        names_temp <- rownames(temp) 
        temp <- as.vector(temp[,1])
        names(temp) <- names_temp
        rand_mat <- as.matrix(sapply(1:100, function(j){ 
            sample(temp, replace=TRUE)
            }))
        new_pval <- sapply(temp , function(i){ 
            length(which(rand_mat > as.numeric(i)))/(100*length(temp)
            )
        })
        new_pval <- replace(new_pval, new_pval == 0, 
                            min(subset(new_pval, new_pval != 0), na.rm=T))
        temp[which(new_pval < alpha)]
    }
)

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
              labels, col=bg, ... )
    }
    text(xy$x, xy$y, labels, col=col, ... )
}

setMethod(
    f="addShadowText", 
    signature="ANY",
    definition=function(x, y=NULL, labels, col='white', bg='black',
                        theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ...) {
        
        xy <- xy.coords(x,y)
        xo <- r*strwidth('A')
        yo <- r*strheight('A')
        for (i in theta) {
            text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
                  labels, col=bg, ... )
        }
        text(xy$x, xy$y, labels, col=col, ... )        
    }
)

setMethod(
    f="SMITEextractModules", 
    signature="PvalueAnnotation", 
    definition=function(annotation,whichModule=NULL){
        if(length(slot(slot(annotation, "score_data"), 
                       "module_output")$modules) == 0) {
            stop("Spinglass or Bionet analysis has not been performed.")
        }
        temp <- slot(slot(annotation,"score_data"),"module_output")$modules
        if(!is.null(whichModule)){
            temp <- temp[whichModule]
        }
        temp
    }
)
