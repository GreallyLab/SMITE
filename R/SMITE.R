##SMITE_1_0_2 10/9/2015


##internal function to perform a Stouffer's test for p value combination

setMethod(
    f="stoufferTest",
    signature="vector",
    definition=function(pvalues, weights)
    {
        if(is.null(weights)){
            weights <- rep(1, length(pvalues))/length(pvalues)
        }

        Zi <- qnorm(1-pvalues/2)
        Z  <- sum(weights*Zi)/sqrt(sum(weights^2))
        new_pvalues <- (1-pnorm(Z))*2
        new_pvalues <- replace(new_pvalues,new_pvalues < .0000000000000001, .0000000000000001)
        new_pvalues <- replace(new_pvalues, new_pvalues > .9999999999999999, .9999999999999999)
        new_pvalues
    }
)

##SMITE Functions

setMethod(
    f="makePvalueAnnotation",
    signature="ANY",
    definition=function(data, other_data=NULL, other_tss_distance=10000,
                        promoter_upstream_distance=1000, promoter_downstream_distance=1000,
                        strand_col=NULL, gene_name_col=NULL)
    {
        ##create a Granges data object
        if(!inherits(data, "GRanges")){
            if(is.null(gene_name_col)){
                stop("Gene name column must be specified if GRanges not given")
            }
            ##if the strand column was not specified auto-detect
            if(is.null(strand_col)){
                strand_col <- which(data[1, ] %in% c("+", "-"))[1]
            }
            data_grange <- GenomicRanges::GRanges(seqnames=data[, 1],
                                   ranges=IRanges::IRanges(start=data[, 2], end=data[, 3]),
                                   name=data[, gene_name_col],
                                   strand=data[, strand_col])
        }
        else {
            data_grange <- data
            data_grange$score <- NULL
        }
        if(any(duplicated(data_grange$name)))
        {
            message("Genes are duplicated.  Removing duplicates")
            data_grange <-
                subset(data_grange,!duplicated(data_grange$name))
        }

        data_grange$feature <- "original"

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
            names(other_tss_distance) <- otherdata_names
            }
        }

        tss <- GenomicRanges::shift(GenomicRanges::flank(data_grange, width=2), 1)
        tss$feature <- "tss"

        if(!is.null(other_data)){
            other <- do.call(c, sapply(1:length(other_data), function(i){

                if(!inherits(other_data[[i]], "GRanges")){
                temp_other<-c(GenomicRanges::GRanges(
                    seqnames=other_data[[i]][, 1],
                    ranges=IRanges::IRanges(start=other_data[[i]][, 2],
                                   end=other_data[[i]][, 3])))
                }
                else {
                temp_other <- other_data[[i]]
                GenomicRanges::mcols(temp_other) <- NULL
                }
                temp_other <- unique(temp_other)


                suppressWarnings(
                    overlap <- findOverlaps(GenomicRanges::flank(data_grange,
                                                 other_tss_distance[otherdata_names[i]],
                                                 start=TRUE), temp_other)
                )
                temp_other <- temp_other[as.numeric(S4Vectors::subjectHits(overlap))]
                temp_other$name <-
                    data_grange$name[S4Vectors::queryHits(overlap)]
                temp_other$feature <- otherdata_names[i]
                temp_other
            })
            )
        }
        promoters_downstream <- GenomicRanges::flank(data_grange, -promoter_downstream_distance,
                                    start=TRUE)
        promoters_upstream <- GenomicRanges::flank(data_grange, 
            promoter_upstream_distance, start=TRUE)
        end(promoters_upstream) <- end(promoters_upstream)+1

        promoters <- GenomicRanges::punion(promoters_upstream, promoters_downstream)
        promoters$name <- data_grange$name
        promoters$feature <- "promoter"

        body <- GenomicRanges::psetdiff(data_grange, promoters_downstream)
        body$name <- data_grange$name
        body$feature <- "body"


        if(!is.null(other_data)){
            suppressWarnings(combined_data <-
                                 c(data_grange, promoters, body, other, tss))
        }
        else{
            suppressWarnings(combined_data <-
                                 c(data_grange, promoters, body, tss))
        }
        combined_data <- split(combined_data, combined_data$name)
        slot(combined_data, "metadata")$params <-
            c(
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
    signature(gene_IDs="character", ID_type="character", ID_convert_to="character"),
    definition=function(gene_IDs, ID_type, ID_convert_to, delim=NULL, verbose=FALSE)
    {
        if(any(duplicated(gene_IDs))){stop(
            "Cannot convert duplicated ids. Please remove duplicates.")
        }

        if(!is.null(delim)){
            gene_IDs <- do.call(rbind, strsplit(gene_IDs, delim))[, 2]
        }

        genes_old <- unique(as.character(gene_IDs))
        gene_IDs <- cbind(gene_IDs, 1:length(gene_IDs))

        if(ID_type == "refseq"){
            genes_old <- subset(genes_old, genes_old %in%
                                    (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egREFSEQ2EG)))

            if(ID_convert_to == "symbol"){

                eg <- unlist(AnnotationDbi::mget(genes_old,
                                                 org.Hs.eg.db::org.Hs.egREFSEQ2EG))
                symbol <- unlist(AnnotationDbi::mget(eg,
                                                     org.Hs.eg.db::org.Hs.egSYMBOL))
                out <- cbind(names(eg), symbol)
            }
        }

        else
            if(ID_type == "ensembleprot"){
                genes_old <- subset(genes_old, genes_old %in%
                            (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egENSEMBLPROT2EG)))

                if(ID_convert_to == "symbol"){
                    eg <- unlist(AnnotationDbi::mget(genes_old,
                                                     org.Hs.eg.db::org.Hs.egREFSEQ2EG))
                    symbol <- unlist(AnnotationDbi::mget(eg,
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(eg), symbol)
                }
            }
        else
            if(ID_type == "uniprot"){
                genes_old <- subset(genes_old, genes_old %in%
                                        (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egUNIPROT)))

                if(ID_convert_to == "symbol"){
                    eg <- unlist(AnnotationDbi::mget(genes_old,
                                                     org.Hs.eg.db::org.Hs.egREFSEQ2EG))
                    symbol <- unlist(AnnotationDbi::mget(eg,
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(eg), symbol)
                }
            }
        else
            if(ID_type == "ensemble"){
                genes_old <- subset(genes_old, genes_old %in%
                                        (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egENSEMBL2EG)))
                if(ID_convert_to == "symbol"){
                    eg <- unlist(AnnotationDbi::mget(genes_old,
                                                     org.Hs.eg.db::org.Hs.egENSEMBL2EG))
                    symbol <- unlist(AnnotationDbi::mget(eg,
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(eg), symbol)
                }
            }
        else
            if(ID_type == "entrez"){
                if(ID_convert_to == "symbol"){

                    symbol <- unlist(AnnotationDbi::mget(genes_old,
                                                         org.Hs.eg.db::org.Hs.egSYMBOL))
                    out <- cbind(names(symbol), symbol)
                }
            }
        else
            if(ID_type == "symbol"){
                genes_old <- subset(genes_old, genes_old %in%
                                        (AnnotationDbi::ls(org.Hs.eg.db::org.Hs.egALIAS2EG)))
                if(ID_convert_to == "entrez"){
                    eg <- unlist(AnnotationDbi::mget(genes_old,
                                                     org.Hs.eg.db::org.Hs.egALIAS2EG))
                    out <- cbind(names(eg), eg)
                }

            }

        out <- merge(gene_IDs, out, by=1, all.x=TRUE)
        out <- out[order(as.numeric(out[, 2])), ]
        out <- subset(out,!duplicated(out[,1]))

        out[, 3]
    })

setMethod(
    f="annotateExpression",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, expr_data, effect_col=NULL, pval_col=NULL){
        if(is.null(effect_col)){
            effect_col=grep("effect|odds|coeff|B", tolower(colnames(expr_data)))
            if(length(effect_col) != 1){
                stop("Cannot determine effect column. Please specify with arg:effect_col")
            }
        }

        if(any(!c(-1,1) %in% unique(sign(expr_data[, effect_col])))){
            message("WARNING: Effects should provide a direction, but these effects
            are all in the same direction.")
        }

        if(is.null(pval_col)){

            pval_col <- grep("pval|p.val|p_val|sig", tolower(colnames(expr_data)))
            if(length(effect_col) != 1){
                stop("Cannot determine p.value column. Please specify with arg:effect_col")
            }
        }
        if(any(expr_data[, pval_col] < 0, expr_data[, pval_col] > 1)){
            stop("P-values must be between 0 and 1")
        }
        temp_pval_col <- expr_data[, pval_col]
        temp_pval_col <- replace(temp_pval_col, temp_pval_col < .0000001,.0000001)
        temp_pval_col <- replace(temp_pval_col, temp_pval_col > .9999999,.9999999)
        expr_data[, pval_col] <- temp_pval_col

        expression_output <- ExpressionSet(as.matrix(expr_data), featureNames=rownames(expr_data))
        phenoData(expression_output) <-
            new("AnnotatedDataFrame",
                data=as.data.frame(exprs(expression_output)[,
                                                            c(effect_col, pval_col)]))
        varLabels(expression_output) <- c("expression_effect",
                                                      "expression_pvalue")
        slot(pvalue_annotation, "expression") <- expression_output
        pvalue_annotation
        })



setMethod(
    f="annotateModification",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, mod_data, weight_by=NULL,
                        weight_by_method="Stouffer", mod_included=NULL,
                        mod_corr=TRUE, mod_type="methylation", verbose=FALSE){
        if(mod_type %in% names(slot(slot(pvalue_annotation,"modifications"),"metadata")$elements))
        {
            stop("Provided data set is already loaded as mod_type")
        }
        
        unique_feature_names <- unique(
            unlist(slot(pvalue_annotation, "annotation"))$feature)
        
        ##no weights provided
        if(missing(weight_by)){
            weight_by <- rep("pval", length(unique_feature_names[!unique_feature_names %in%
                                                                     c("original", "tss")]))
        }
        
        ##no mod_included or weight names
        if(is.null(names(weight_by))){
            if(is.null(mod_included)){
                mod_included <- unique_feature_names[!unique_feature_names %in%
                                                         c("original", "tss")]
            }
            names(weight_by) <- mod_included
        }
        
        ##weight names were provided
        if(!is.null(names(weight_by))){
            mod_included <- names(weight_by)
            if(!all(mod_included %in% unique_feature_names)){
                stop("Provided weight names must match those in
                     unique(GenomicRanges::mcols(unlist(pvalue_annotation@annotation))$feature)")
            }
            }
        
        if(any(!c(-1, 1) %in% unique(sign(mod_data[, 4])))){
            message("WARNING: Effects should provide a direction,
                    but these effects are all in the same direction.")
        }
        
        if(any(mod_data[, 5] < 0,mod_data[, 5] > 1)){
            stop("P-values must be between 0 and 1")
        }
        
        
        mod_grange <- GenomicRanges::GRanges(seqnames=mod_data[, 1],
                                             ranges=IRanges::IRanges(start=mod_data[, 2], end=mod_data[, 3]),
                                             effect=mod_data[, 4], pval=mod_data[, 5], type=mod_type)
        
        temp_annotation <- unlist(slot(pvalue_annotation, "annotation"))
        
        overlap_mods <- GenomicRanges::findOverlaps(temp_annotation, mod_grange)
        mod_grange_overlaps <- mod_grange[S4Vectors::subjectHits(overlap_mods)]
        GenomicRanges::mcols(mod_grange_overlaps) <- cbind(GenomicRanges::mcols(
            temp_annotation[as.numeric(
                S4Vectors::queryHits(overlap_mods))]),
            GenomicRanges::mcols(mod_grange_overlaps))
        mod_grange_overlaps <- split(mod_grange_overlaps, mod_grange_overlaps$name)
        temp_annotation <- split(temp_annotation, temp_annotation$name)
        
        
        if(mod_corr == TRUE){
            if(verbose == TRUE){
                message("Computing correlation matrices")
            }
            temp_split_mod_grange <- split(mod_grange, GenomicRanges::seqnames(mod_grange))
            precede_follow_each_element <- lapply(temp_split_mod_grange,
                                                  function(chr){
                                                      temp_chr <- IRanges(start(chr), end(chr))
                                                      temp_precede <- precede(temp_chr)
                                                      temp_follow <- follow(temp_chr)
                                                      temp_precede[which(is.na(temp_precede))] <-
                                                          which(is.na(temp_precede))
                                                      temp_follow[which(is.na(temp_follow))] <-
                                                          which(is.na(temp_follow))
                                                      chr[c(temp_follow, temp_precede)]
                                                  })
            mod_grange_corr <- unlist(GRangesList(precede_follow_each_element))
            
            duplicate_each_chr <- lapply(temp_split_mod_grange, function(chr){
                c(chr,chr)
            })
            duplicate_each_chr <- unlist(GRangesList(duplicate_each_chr))
            mod_grange_corr$distance <- IRanges::distance(duplicate_each_chr,
                                                          mod_grange_corr)
            
            mod_grange_corr$pval2 <- duplicate_each_chr$pval
            mod_grange_corr <- mod_grange_corr[which(mod_grange_corr$pval2<0.05)]
            quantile_distances_mod_corr <- Hmisc::cut2(mod_grange_corr$distance,
                                                       g=500, onlycuts=TRUE)
            quantile_distances_mod_corr[length(quantile_distances_mod_corr)] <- 250000000
            mod_grange_corr$cat <- cut(mod_grange_corr$distance, breaks=quantile_distances_mod_corr)
            mod_grange_corr <- split(mod_grange_corr, mod_grange_corr$cat)
            mod_grange_corr2 <- lapply(mod_grange_corr, function(j){
                mean((sapply(1:500, function(i){
                    index <- sample(1:length(j), replace=TRUE);
                    cor(qnorm(1-j$pval[index]), qnorm(1-j$pval2[index]))
                })))
            })
            correlations <- as.data.frame(do.call(rbind, mod_grange_corr2))
            final_corr <- data.frame(correlations,
                                     as.character(names(mod_grange_corr2)),
                                     stringsAsFactors=FALSE)
            final_corr <- rbind(c(.9, paste("(-1, ", quantile_distances_mod_corr[1], "]",
                                            sep="")), final_corr)
            rm(mod_grange_corr)
        }
        combined_pvalues_list <- sapply(mod_included, function(i){
            if(verbose == TRUE){
                message(paste("Combining p-values over:", i))
            }
            temp <- subset(unlist(mod_grange_overlaps), unlist(mod_grange_overlaps)$feature == i)
            
            
            ref_data <- unlist(slot(pvalue_annotation, "annotation"))
            ref_data <-
                subset(ref_data,
                       ref_data$feature == "tss")
            ref_data <- ref_data[temp$name]
            suppressWarnings(temp$distance <- distance(ref_data,
                                                       temp)+2)
            temp <- split(temp, temp$name)
            forreturn <- lapply(temp, function(each){
                each_length <- length(each)
                each_effect <- each$effect[order(each$pval)]
                each_pval <- each$pval[order(each$pval)]
                distances <- each$distance[order(each$pval)]
                if(length(each_pval > 1)){
                    if(mod_corr == TRUE){
                        corr_mat <- matrix(as.numeric(final_corr[match(cut(
                            as.matrix(dist(start(each)[order(each$pval)])),
                            breaks=c(-1, quantile_distances_mod_corr)), final_corr[, 2]), 1]),
                            ncol=each_length)
                        diag(corr_mat) <- 1
                        corr_mat<-abs(corr_mat)
                        chol_d <- try(chol(corr_mat), silent=TRUE)
                        while(is(chol_d, "try-error"))
                        {
                            index<- as.numeric(strsplit(
                                strsplit(chol_d[1],
                                         "the leading minor of order ")[[1]][2],
                                " is not positive")[[1]][1])-1
                            chol_d <- try(chol(corr_mat[1:index, 1:index]),silent=TRUE)
                            each_pval <- each_pval[1:index]
                            each_effect <- each_effect[1:index]
                            distances <- distances[1:index]
                            each_length <- index
                        }
                        
                        each_pval <- 1-pnorm(abs(solve(t(chol_d)) %*% qnorm(
                            1-each_pval/2)))
                        each_pval<-replace(each_pval, each_pval == 0 , 0.000000001)
                        each_pval <- each_pval*2
                    }
                    
                    if(weight_by_method == "Stouffer"){
                        
                        if(weight_by[i] == "distance"){
                            ##mean is weighted by distance
                            out_mean <- weighted.mean(each_effect,
                                                      w=(1/log(distances)))
                            ##Stouffer test is weighted by distance
                            out_pval <- stoufferTest(each_pval, weights=(1/log(distances)))
                        }
                        else if(weight_by[i] %in%
                                c("pval", "p.value", "pvalue", "p_val")){
                            ##mean is weight by pvalue
                            out_mean <- weighted.mean(each_effect, w=-log(each_pval))
                            out_pval <- stoufferTest(each_pval, weights=NULL)
                        }
                        else {
                            ##mean is not weighted
                            out_mean <- mean(each_effect, na.rm=TRUE)
                            out_pval <- stoufferTest(each_pval, weights=NULL)
                        }
                        
                    }
                    else if(weight_by_method %in% c("minimum", "Sidak", "sidak")){
                        
                        index <- which(each_pval == min(each_pval))
                        if(length(index) > 1){
                            index <- index[which(
                                abs(each_effect[index]) == max(abs(each_effect[index])))][1]
                        }
                        out_mean <- each_effect[index]
                        out_pval <- 1-(1-each_pval[index])^length(each_pval)
                        
                    }
                    else if(weight_by_method == "binomial"){
                        
                        index <- which(each_pval == min(each_pval))
                        if(length(index) > 1){
                            index <- index[which(abs(each_effect[index]) == max(
                                abs(each_effect[index])))][1]
                        }
                        out_mean <- each_effect[index]
                        out_pval <- (1-pbinom(q=length(which(each_pval<0.05)),
                                              size=each_length, prob=0.05))
                    } else if(weight_by_method %in%
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
        combined_pvalues_list <- lapply(combined_pvalues_list, function(each_feature){
            categories <- data.frame(categories=as.numeric(
                Hmisc::cut2(each_feature[, 3], g=100)))
            categories_table <- data.frame(table(categories))
            
            trans_p <- cbind(trans=qnorm(1-each_feature[, 2]/2),
                             plyr::join(categories, categories_table, by="categories"))

            trans_p[, 1] <- replace(trans_p[, 1],is.infinite(trans_p[, 1]),
                                    max(subset(trans_p, !is.infinite(
                                        trans_p[, 1])), na.rm=TRUE))
            
            
            num_list <- split(trans_p$trans, trans_p$categories)
            rand_list <- sapply(1:length(num_list), function(i){
                as.matrix(sapply(1:500, function(j){
                    sample(num_list[[as.numeric(i)]], replace=TRUE)
                }))
            })
            new_pval <- apply(trans_p, 1, function(i){
                length(which(
                    rand_list[[as.numeric(i[2])]] > as.numeric(i[1])))/
                    (500*as.numeric(i[3]))
            })
            new_pval <- replace(new_pval, new_pval == 0,
                                min(subset(new_pval, new_pval != 0),
                                    na.rm=TRUE))
            
            each_feature[, 2] <- new_pval
            each_feature <- as.data.frame(each_feature)
            each_feature
        })
        
        if(verbose == TRUE){
            message("Scores have been adjusted")
        }
        
        newmods <- c(unlist(slot(pvalue_annotation, "modifications")),
                     unlist(mod_grange_overlaps))
        names(newmods) <- NULL
        newmods <- split(newmods, newmods$name)
        
        output_m_summary <-
            suppressWarnings(as.data.frame(c(list(names=names(mod_grange_overlaps)),
                                             lapply(combined_pvalues_list, function(x){
                                                 x[match(names(
                                                     mod_grange_overlaps),
                                                     rownames(x)), 1:2]
                                             }))))
        rownames(output_m_summary) <- output_m_summary[, 1]
        output_m_summary <- output_m_summary[, -1]
        colnames(output_m_summary) <- paste(mod_type,
                                            apply(expand.grid(c("effect", "pvalue"),
                                                              mod_included),1, function(i){
                                                                  paste(i[2], i[1],
                                                                        sep="_")
                                                              }), sep="_")
        
        newmetadata <- slot(slot(pvalue_annotation, "modifications"), "metadata")
        if(is.null(newmetadata$m_summary)){
            newmetadata$m_summary <- output_m_summary
        }
        else{
            newmetadata$m_summary <- merge(newmetadata$m_summary,
                                           output_m_summary, by=0, all=TRUE)
            rownames(newmetadata$m_summary) <- newmetadata$m_summary[, 1]
            newmetadata$m_summary <- newmetadata$m_summary[, -1]
        }
        newmetadata[["elements"]][[mod_type]]$weight_by <- weight_by
        newmetadata[["elements"]][[mod_type]]$weight_by_method <- weight_by_method
        newmetadata$elementnames <- c(newmetadata$elementnames,
                                      paste(mod_type, mod_included, sep="_"))
        
        slot(newmods, "metadata") <- newmetadata
        slot(pvalue_annotation, "modifications") <- newmods
        pvalue_annotation
        }
)

setMethod(
    f="removeModification",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, mod_type="methylation"){

        if(!mod_type%in%names(slot(pvalue_annotation,
                                  "modifications")@metadata$elements)){
            stop("Provided mod_type is not in the pvalue_annotation")
        }

        temp_meta <- slot(slot(pvalue_annotation,"modifications"),"metadata")
        temp <- unlist(slot(pvalue_annotation,"modifications"))
        names(temp) <- NULL
        temp <- subset(temp,!(temp$type == mod_type))
        slot(temp,"metadata") <- temp_meta
        temp_meta_colnams <- colnames(slot(temp,"metadata")$m_summary)
        slot(temp,"metadata")$m_summary <-
            slot(temp,"metadata")$m_summary[, -grep(mod_type, temp_meta_colnams)]

        if(ncol(slot(temp,"metadata")$m_summary) == 0){
            slot(temp,"metadata")$m_summary <- NULL
        }

        slot(temp,"metadata")$elements[which(
            names(slot(temp,"metadata")$elements) == mod_type)] <- NULL


        slot(temp,"metadata")$elementnames <-
            slot(temp,"metadata")$elementnames[
                -which(do.call(rbind, lapply(strsplit(
                    slot(temp,"metadata")$elementnames, "_"),
                    function(i)i[1]))[, 1] %in% mod_type)]

        temp_meta <- slot(temp,"metadata")
        slot(pvalue_annotation,"modifications") <- split(temp, temp$name)
        slot(slot(pvalue_annotation,"modifications"),"metadata") <- temp_meta
        pvalue_annotation
    }
)

setMethod(
    f="makePvalueObject",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, effect_directions=NULL) {


        if(is.null(effect_directions)){
            effect_directions <- rep("bidirectional",
                                     length(slot(
                                                slot(pvalue_annotation, "modifications"),
                                                "metadata")$elementnames))
        }

        if(is.null(names(effect_directions))){
            names(effect_directions)<-slot(slot(pvalue_annotation, "modifications"),
                                           "metadata")$elementnames
        }

        if(any(!names(effect_directions) %in% slot(slot(pvalue_annotation, "modifications"),
                                                 "metadata")$elementnames)){
            stop("Effect name is invalid")
        }

        if(any(!(effect_directions %in% c("decrease", "increase", "bidirectional")
        ))){
            stop("Effect argument is invalid.")
        }


        exp_ind <- ifelse(nrow(Biobase::pData(pvalue_annotation@expression)) >
                              0, 1, 0)

        total_num_factor <- length(effect_directions)

        signs_index <- merge(cbind(c("increase", "decrease", "bidirectional"),
                                c(1, -1, 2)),
                          cbind(effect_directions, 1:total_num_factor), by=1)
        signs_index <- signs_index[order(signs_index[, 3]), ]
        rownames(signs_index) <- NULL
        colnames(signs_index) <- c("expression_relationship", "B_coeff","name")

        temp1 <- annotationOutput(pvalue_annotation)
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
        signs_index[, 3] <- names(effect_directions)
        slot(pvalue_annotation, "score_data") <- new(Class="PvalueObject",
                                        pval_data=data[, grep("pval",
                                                              colnames(data))],
                                        effect_data=data[, grep("effect",
                                                                colnames(data))],
                                        genes=genenames,
                                        signs_index=signs_index)
        pvalue_annotation
    }
)




setMethod(
    f="normalizePval",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, trans, ref="expression_pvalue",
                        method="rescale"){

        temp_pval_data <- slot(slot(pvalue_annotation,"score_data"),"pval_data")

        names_temp_pval_data <- colnames(temp_pval_data)

        if(nrow(temp_pval_data) == 0){
            stop("Run makePvalueObject function first.")
        }
        if(!any(grepl(ref,names_temp_pval_data))){
               stop("paste(Reference is not one of the available:",
                    names_temp_pval_data)
        }
        if(length(grep(ref,names_temp_pval_data)) > 1){
               stop("Reference was not specific enough.")
        }

        ref_index <- grep(ref,names_temp_pval_data)
        p_ref <- temp_pval_data[[ref_index]]
        logit_ref <- log(p_ref/(1-p_ref))
        temp_signs_index<-slot(slot(pvalue_annotation,"score_data"),
                               "signs_index")[, 3]

        if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
            dev.new(height=7, width=14)
        }

        par(mfrow=c(1, 2))
        plotDensityPval(pvalue_annotation, ref=ref)

        if(method %in% c("Box-Cox", "box-cox", "boxcox", "Boxcox")){
            if(missing(trans)){
                message(paste("Auto-detecting best transformation"))
                optimal_boxcox_exponent <- c()
                for(x in names_temp_pval_data[-ref_index]){

                    p_temp <- temp_pval_data[[x]]
                    if(!all(is.na(p_temp))){
                        logit_temp <- log(p_temp/(1-p_temp))
                        nonparametric_comparison <- t(sapply(c(seq(.05, .95, .05),
                                           rev(1/seq(.05, .95, .05))),
                                           function(i){
                                               c(i,wilcox.test(logit_ref,
                                                        as.numeric(logit_temp)*i)$p_value)
                                               }))
                        nonparametric_comparison <-  subset(nonparametric_comparison, 
                        nonparametric_comparison[, 2] == max( 
                        nonparametric_comparison[, 2])[1])[, 1]
                        p_temp <-
                               (exp(logit_temp)^( nonparametric_comparison))/(1+exp(logit_temp)^( nonparametric_comparison))
                    }
                    else {
                        nonparametric_comparison <- 1
                    }
                       optimal_boxcox_exponent <- c(optimal_boxcox_exponent,  nonparametric_comparison)
                       names(optimal_boxcox_exponent)[length(optimal_boxcox_exponent)] <- x
                       temp_pval_data[[x]] <- p_temp

                }
            }

            else {



                if(length(trans) !=
                       length(grep(paste(temp_signs_index, collapse="|"),
                                   colnames(temp_pval_data)))){
                       stop("Length of p and transformations must equal!")
                }
                else {
                    optimal_boxcox_exponent <- trans
                    names(optimal_boxcox_exponent) <- subset(names_temp_pval_data,
                                                             !names_temp_pval_data %in% ref)
                }
                for(x in names(optimal_boxcox_exponent)){
                    p_temp <- temp_pval_data[[x]]
                    if(!all(is.na(p_temp))){
                           p_temp <- p_temp^optimal_boxcox_exponent[x]
                    }
                    temp_pval_data[[x]] <- p_temp
                }
            }
        }
        else if(method%in%c("Rescale", "rescale")){
            for(i in temp_signs_index){
                if(!all(is.na(slot(
                    slot(pvalue_annotation, "score_data"), "pval_data")[[grep(i,
                                                                   names_temp_pval_data)]])
                   )){
                    p_temp <- temp_pval_data[[grep(i, names_temp_pval_data)]]
                    logit_temp <- log(p_temp/(1-p_temp))
                    logit_temp <- scales::rescale(logit_temp, to=range(logit_ref,
                                                                    na.rm=TRUE))
                    temp_pval_data[[grep(i,
                                           names_temp_pval_data
                                           )]] <- exp(logit_temp)/(1+exp(logit_temp))
                }
            }
        }

        slot(slot(pvalue_annotation,"score_data"),"pval_data") <- temp_pval_data
        plotDensityPval(pvalue_annotation, ref=ref)
        pvalue_annotation
    }
)

setMethod(
    f="scorePval",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, weights){

        total_num_factor <- ncol(slot(slot(pvalue_annotation, "score_data"), "pval_data"))
        pval_score_colnames <- colnames(slot(slot(pvalue_annotation,"score_data"),
                                             "pval_data"))
        if(missing(weights)){
            weights <- rep((1/(total_num_factor)), total_num_factor)
            names(weights) <- pval_score_colnames
        }
        else {
            if(length(weights) != total_num_factor){
                stop("ERROR: Number of factors(with expression data)
                     and weights must equal!")
            }
            else {

                if(is.null(names(weights))){
                    names(weights) <- pval_score_colnames
                }
                sorted_pval_score_colnames <- sort(do.call(rbind, strsplit(pval_score_colnames,
                                             "_pvalue")))
                if(!all(sort(names(weights)) == sorted_pval_score_colnames))
                {
                    stop(paste("Weight names must match the following:",
                               paste(sorted_pval_score_colnames,
                                     collapse=", ")))
                }
            }
        }

        weights <- weights[match(
            do.call(rbind, strsplit(pval_score_colnames, "_pvalue")),
            names(weights))]

        message("The following weights are being applied")
        print(weights)
        temp_score_data <- slot(pvalue_annotation, "score_data")
        temp_pval_data <- slot(temp_score_data, "pval_data")
        temp_effect_data <- slot(temp_score_data, "effect_data")
        temp_signs_index <- slot(temp_score_data, "signs_index")
        scoringdata <- temp_pval_data*sign(temp_effect_data)

        weight_names_temp <- names(weights)

        if(any(grepl("exp", names(weights)))){
            weight_names_temp <- weight_names_temp[-grep("exp", names(weights))]
        }
        slot(temp_score_data, "scoring_vector") <- weights

        unidirectional <- as.numeric(temp_signs_index[match(temp_signs_index[, 3], weight_names_temp), 2])
        bidirectional <- unidirectional
        unidirectional[which(unidirectional == 2)] <- 0
        bidirectional[which(bidirectional != 2)] <- 0
        bidirectional[which(bidirectional == 2)] <- 1


        scoringdata <- qnorm(1-as.matrix(abs(scoringdata))/2)*sign(scoringdata)
        scoresout <- apply(scoringdata, 1, function(each){
            # for each gene
            if(any(!is.na(each)))
            {
              # not all missing
                if(any(grepl("exp", names(each)))){
                    # there is expression data
                    exp_index <- grep("exp", names(each))
                    forreturn <- (sum(
                        abs(sum(c(as.numeric(each[[exp_index]]),
                                  as.numeric(each[-exp_index]))*
                                    c(1, unidirectional)*
                                    weights[c(exp_index, 
                                              which(!grepl("exp", names(each))))],
                                na.rm=TRUE)),
                        (abs(as.numeric(each[-exp_index])*bidirectional)*
                             weights[-exp_index]),
                        na.rm=TRUE)/sum(weights^2)^.5)
                }
                else {

                    # there is no expression data

                    forreturn <- (
                        sum(abs(sum(as.numeric(each)*unidirectional*weights,
                                    na.rm=TRUE)), (abs(as.numeric(each)*
                                                           bidirectional)*
                                                       weights),
                            na.rm=TRUE)/sum(weights^2)^.5)
                }
            }
            else {
                forreturn <- (NA)
            }
            forreturn
        })

        scoresout <- (1-pnorm(as.numeric(scoresout)))*2
        iscoresout<-replace(scoresout, scoresout == 0 ,  min(subset(scoresout,
                 !(scoresout == 0)), na.rm=TRUE))
                scoresout <- (-2*log(scoresout))
        rand_mat <- as.matrix(sapply(1:100, function(j){
            sample(scoresout, replace=TRUE)
            }))
        new_pval <- sapply(scoresout , function(i){
            length(which(rand_mat > as.numeric(i)))/(100*length(scoresout))
        })
        new_pval <- replace(new_pval, new_pval == 0,
                            min(subset(new_pval, new_pval!=0), na.rm=TRUE))
        new_pval <- (-2)*log(new_pval)



        slot(temp_score_data, "scores") <-
            data.frame(scores=as.numeric(new_pval),
                       row.names=as.character(slot(temp_score_data,
                                                   "genes")))
        slot(pvalue_annotation, "score_data") <- temp_score_data
    pvalue_annotation
    }
)

setMethod(
    f="runSpinglass",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, network, random_alpha = 0.05, gam = 0.5,
                        node_alpha = 0.05, maxsize = 500, minsize = 8,
                        num_iterations = 1000, simplify=TRUE)
    {

        if(inherits(network, what="graphNEL")){
            network <- graph_from_graphnel(network)
        }

        if(length(slot(slot(pvalue_annotation,"score_data"),"module_output")) !=0 ){
            slot(slot(pvalue_annotation,"score_data"),"module_output") <- list()
            message("Overwriting existing modules.")
        }

        if(simplify == TRUE){
            network <- igraph::simplify(network, remove.multiple=TRUE,
                                        remove.loops=TRUE)
        }
        genes_in_network <- subset(slot(slot(pvalue_annotation, "score_data"),
                                        "genes"),
                                   slot(slot(pvalue_annotation, "score_data"),
                                        "genes") %in%
                                       igraph::V(network)$name)
        scores_in_network <- extractScores(pvalue_annotation)[genes_in_network]

        ##should be FALSE, but just in case check for NAs
        if(any(is.na(scores_in_network))){
            ##if there are NAs remove them
            scores_in_network <- subset(scores_in_network, !is.na(
                scores_in_network))
            genes_in_network <- names(scores_in_network)
        }

        nodes_with_scores <- base::intersect(genes_in_network, igraph::V(network)$name) 
        network <- igraph::induced_subgraph(network, nodes_with_scores)
        network_clusters <- igraph::clusters(network)

        ## choose largest connected nodes ## may want include > minsize
        maxclust <- which(network_clusters$csize ==
                              max(network_clusters$csize))[1]
        network <- igraph::induced_subgraph(network,
                                    which(network_clusters$membership == maxclust))
        rm(network_clusters)
        genes_in_network <- intersect(genes_in_network, igraph::V(network)$name) 
        
        scores_in_network <- scores_in_network[genes_in_network]
        network.adj <- igraph::as_adjacency_matrix(network)

        ## order of scores has to match order of adjacency rows
        scores_in_network <- scores_in_network[rownames(network.adj)]
        genes_in_network <- names(scores_in_network)
        stat_scores <- as.numeric(scores_in_network)
        pval_scores <- exp(scores_in_network/(-2))

        network_with_scores <- apply(network.adj, 1, function(v) v*stat_scores)
        W_for_spinglass <- (network_with_scores + t(network_with_scores))
        rm(network_with_scores)
        gc()
        W_vec <- (-2*log(1-pchisq(as.vector(W_for_spinglass),4)))

        W_vec <- replace(W_vec, is.infinite(W_vec),
                         max(subset(W_vec, !is.infinite(W_vec))) )

        W_for_spinglass <- matrix(W_vec, nrow=nrow(W_for_spinglass))
        rm(W_vec)
        rownames(W_for_spinglass) <- genes_in_network
        colnames(W_for_spinglass) <- genes_in_network
        gc()
        final_network <-
            igraph::graph_from_adjacency_matrix(W_for_spinglass,
                                                mode = "undirected",
                                                weighted=TRUE)
        igraph::V(final_network)$weight <- stat_scores
        network <- final_network
        rm(final_network)
        gc()

        ## For significant genes apply Spinglass algorithm
        sig_genes <- subset(genes_in_network, pval_scores < node_alpha)
        sig_genes_counter <- 1:length(sig_genes)
        names(sig_genes_counter) <- sig_genes
        spin_glass_out <- lapply(sig_genes, function(j){
            message(paste("Computing modules: Vertex",sig_genes_counter[j],"of",
                          length(sig_genes), "significant genes is", j))
            genes_in_network[
                cluster_spinglass(network, weights = igraph::E(network)$weight,
                                  vertex=j, gamma=gam)$community]
            })
        names(spin_glass_out) <- sig_genes

        ## Select modules with size requirements
        spin.size <- do.call(c, lapply(spin_glass_out, length))
        spin_glass_out <- subset(spin_glass_out, spin.size >= minsize & spin.size
                                <= maxsize)

        Modularity.edges = function(v, network)
        {
            h <- igraph::induced_subgraph(network, v);
            c(sum(igraph::E(h)$weight))
        }
        edge_sum <- do.call(c,lapply(spin_glass_out, function(j) {
            Modularity.edges(j, network)
        }));

        nspin_glass_out <- length(spin_glass_out);
        random_edges <- lapply(1:nspin_glass_out, function(i) {
            each_spin_result <- spin_glass_out[[i]]
            subnetwork <- igraph::induced_subgraph(network, each_spin_result);
            adjacency_of_subnetwork <- igraph::as_adjacency_matrix(subnetwork, sparse=FALSE);
            sapply(1:num_iterations, function(k){
                message(paste("Testing significance: module", i, "of",
                              nspin_glass_out, "Randomization", k, "of", num_iterations))
                random_sample_of_scores = sample(stat_scores, nrow(adjacency_of_subnetwork) , replace=TRUE)
                random_network_with_scores = apply(adjacency_of_subnetwork, 1, function(v) v*random_sample_of_scores)
                W_random <- (random_network_with_scores + t(random_network_with_scores));
                W_random <- apply(W_random, 2, function(i){
                    replace(i, i>0, (-2*log(1-pchisq(subset(i,i>0),4))))
                    })
                sum(W_random)/2
            })
        })
        names(random_edges) <- names(spin_glass_out);

        random_p <- lapply(1:nspin_glass_out, function(k){
            length(which(random_edges[[k]] >
                             edge_sum[k]))/num_iterations})
        names(random_p) <- names(spin_glass_out)

        if(length(spin_glass_out[which(do.call(c,random_p) < random_alpha)]) == 0){
            stop("No modules found.  Please adjust the random_alpha and node_alpha
                 parameters")
        }

        ## Assemble output data
        output <- list();
        output[[1]] <- subset(spin_glass_out, do.call(c,random_p) < random_alpha);
        output[[2]] <- pval_scores;
        output[[3]] <- stat_scores;
        output[[4]] <- igraph::induced_subgraph(network, as.character(
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

        slot(slot(pvalue_annotation, "score_data"), "module_output") <- output
        pvalue_annotation
    }
)


setMethod(
    f="runBioNet",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, network, alpha = 0.05)
    {
        if(any(alpha<0, alpha>1)){
            stop("Error: alpha must be between zero and one.")
        }
        if(length(slot(slot(pvalue_annotation,"score_data"),"module_output")) !=0 ){
            slot(slot(pvalue_annotation,"score_data"),"module_output") <- list()
            message("Overwriting existing modules.")
        }
        
        scores <- highScores(pvalue_annotation, alpha=alpha)
        pval.v <- exp(scores/(-2))
        g <- suppressWarnings(subNetwork(names(pval.v), network))
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
        slot(slot(pvalue_annotation,"score_data"),"module_output") <- output
        pvalue_annotation
    }
)

setMethod(
    f="runGOseq",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, p_thresh=0.05, coverage, type="reactome")
    {
        while(!names(dev.cur()) %in% c("pdf","null device")){
            dev.off()
        }
        names.eid <- names(slot(slot(pvalue_annotation, "score_data"),
                              "module_output")$modules)
        eid <- slot(slot(pvalue_annotation, "score_data"), "module_output")$modules
        pval <- slot(slot(pvalue_annotation, "score_data"), "module_output")$p_value
        stat <- slot(slot(pvalue_annotation, "score_data"), "module_output")$statistic
        genes <- names(slot(slot(pvalue_annotation, "score_data"), "module_output")$modules)
        goseqdata <- slot(slot(pvalue_annotation, "score_data"), "module_output")$modules

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

        slot(slot(pvalue_annotation, "score_data"), "module_output")$goseqOut <-
            lapply(slot(slot(pvalue_annotation, "score_data"),
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
        pvalue_annotation
    }
)




setMethod(
    f="searchGOseq",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, search_string, wholeword=FALSE){
        options(stringsAsFactors=FALSE)
        search_string<-tolower(search_string)
        nums <- which(do.call("c", lapply(
            slot(slot(pvalue_annotation, "score_data"), "module_output")$goseqOut,
            function(each){
                any(grepl(ifelse(wholeword == FALSE, search_string,
                                 paste("\\b", search_string, "\\b", sep="")),
                          tolower(each[, ncol(each)])))
                }
            ))
        )
        if(length(nums) > 0){
            out <- lapply(nums, function(i){
                pos <- grep(ifelse(wholeword == FALSE, search_string,
                                   paste("\\b", search_string, "\\b", sep="")),
                            tolower(slot(slot(pvalue_annotation, "score_data"),
                                         "module_output")$goseqOut[[i]][, 6]))
                tot <- nrow(slot(slot(pvalue_annotation, "score_data"),
                                 "module_output")$goseqOut[[i]])
                outterm <- as.data.frame(as.character(slot(
                    slot(pvalue_annotation, "score_data"),
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
                    slot(pvalue_annotation, "score_data"),
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
    f="extractGOseq",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, which_network=NULL){
        temp <- slot(slot(pvalue_annotation, "score_data"), "module_output")$goseqOut
        if(!is.null(which_network)){temp <- temp[which_network]}
        temp
    }
)






setMethod(
    f="extractModification",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, mod_type=NULL){
        if(!is.null(mod_type)){
            if(!mod_type%in%names(slot(pvalue_annotation,
                                      "modifications")@metadata$elements)){
                stop("Provided mod_type is not in the PvalueAnnotation object.")
            }
        }

        temp_meta <- slot(slot(pvalue_annotation,"modifications"),"metadata")
        temp <- unlist(slot(pvalue_annotation,"modifications"))
        names(temp) <- NULL

        if(is.null(mod_type)){mod_type <- unique(temp$type)}
        temp<-subset(temp, temp$type %in% mod_type)
        slot(temp, "metadata") <- temp_meta
        temp
    }
)

setMethod(
    f="extractExpression",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation){
        if(is.null(expression)){
            stop("No expression data loaded.")
        }
        temp_exp <- slot(pvalue_annotation,"expression")
        Biobase::pData(temp_exp)
    }
)




setMethod(
    f="extractModSummary",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation){
        temp_meta<-slot(slot(pvalue_annotation,"modifications"),"metadata")
        temp_meta$m_summary
    }
)

setMethod(
    f="extractScores",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation){
        if(nrow(slot(slot(pvalue_annotation,"score_data"),"scores")) == 0){
            stop("Run scorePval function first.")
        }
        temp <- slot(slot(pvalue_annotation,"score_data"),"scores")
        names_temp <- rownames(temp)
        temp <- as.vector(temp[,1])
        names(temp) <- names_temp
        temp
    }
)

setMethod(
    f="highScores",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, alpha=0.05){
        if(nrow(slot(slot(pvalue_annotation,"score_data"),"scores")) == 0){
            stop("Run scorePval function first.")
        }
        if(any(alpha < 0, alpha > 1)){
            stop("alpha must be between 0 and 1")
        }

        temp <- slot(slot(pvalue_annotation,"score_data"),"scores")
        names_temp <- rownames(temp)
        temp <- as.vector(temp[,1])
        names(temp) <- names_temp
        rand_mat <- as.matrix(sapply(1:100, function(j){
            sample(temp, replace=TRUE)
            }))
        new_pval <- sapply(temp , function(i){
            length(which(rand_mat > as.numeric(i)))/(100*length(temp))
        })
        new_pval <- replace(new_pval, new_pval == 0,
                            min(subset(new_pval, new_pval != 0), na.rm=TRUE))
        temp[which(new_pval < alpha)]
    }
)

setMethod(
    f="addShadowText",
    signature="ANY",
    definition=function(x, y=NULL, labels, col='white', bg='black',
                        theta=seq(pi/4, 2*pi, length.out=8), r=0.1, ...) {

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
    f="extractModules",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, which_module=NULL){
        if(length(slot(slot(pvalue_annotation, "score_data"),
                       "module_output")$modules) == 0) {
            stop("Spinglass or Bionet analysis has not been performed.")
        }
        temp <- slot(slot(pvalue_annotation,"score_data"),"module_output")$modules
        if(!is.null(which_module)){
            temp <- temp[which_module]
        }
        temp
    }
)
