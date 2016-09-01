setMethod(
    f="plotCompareScores",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, x_name, y_name, ...)
    {
        pval_data <- slot(slot(pvalue_annotation, "score_data"), "pval_data")
        pval_col_names <- colnames(pval_data)
        effect_data <- slot(slot(pvalue_annotation, "score_data"), "effect_data")
        effect_col_names <- colnames(effect_data)

        if(all(!grepl(x_name, pval_col_names)))
        {
            stop(paste("Provided x_name was not one of the following:",
                    paste(pval_col_names, collapse=",")))
        }

        if(all(!grepl(y_name, pval_col_names)))
        {
            stop(paste("Provided y_name was not one of the following:",
                    paste(pval_col_names, collapse=",")))
        }

        if(length(grep(x_name, pval_col_names)) > 1)
        {
            stop(paste("Provided x_name was not unique"))
        }

        if(length(grep(y_name, colnames(pval_col_names))) > 1)
        {
            stop(paste("Provided y_name was not unique"))
        }

        x <- (-log(abs(as.numeric(pval_data[[grep(x_name, pval_col_names)]])))*
                  sign(effect_data[[grep(x_name, effect_col_names)]]))

        y <- (-log(abs(as.numeric(pval_data[[grep(y_name, pval_col_names)]])))*
                  sign(effect_data[[grep(y_name, effect_col_names)]]))

        ggplot(data=data.frame(x=x, y=y), aes(x=x, y=y))+
            stat_binhex(bins=50)+
            geom_hline(yintercept=0, colour="red", linetype = "longdash")+
            geom_vline(xintercept=0, colour="red", linetype = "longdash")+
            labs(title = paste(y_name, "vs", x_name,
                               "p-values comparing effect direction"),
            x=paste(x_name, "Score * Effect Direction"), y=paste(y_name,
            "Score * Effect Direction"))

    }
)


setMethod(
    f="plotModule",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, p_thresh=0.05, which_network=1, goseq=FALSE,
                        layout="fr", legend=TRUE, namestyle="symbol",
                        suppress_details=FALSE, meth_hi_col="blue",
                        meth_low_col="yellow1", meth_mid_col="gray90",
                        exp_hi_col="red1", exp_low_col="chartreuse1",
                        exp_mid_col="gray90", label_scale=TRUE,
                        compare_plot=FALSE, pdf_out=NULL)
    {
        
        module_output <- slot(slot(pvalue_annotation, "score_data"), "module_output")
        
        if(!is.null(pdf_out)){
            unlink(pdf_out)
            if(compare_plot==TRUE){
                pdf(pdf_out, width=24, height=8.5)
            }
            else {
                pdf(pdf_out, width=16, height=8.5)
            }
        }
        
        ## if two side by side plots are not needed
        if(compare_plot == FALSE){
            ## no pdf
            par(mfrow=c(1,1))
            if(is.null(pdf_out)){
                ##no goseq
                if(goseq == FALSE){
                    ##yes legend
                    if(legend == TRUE){
                        if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
                            dev.new(height=10, width=12)
                        }
                    }
                    ##no legend
                    else{
                        if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
                            dev.new(height=8, width=8)
                        }
                    }
                }
                ##yes goseq and legend
                else {
                    if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
                        dev.new(height=10, width=16)
                    }
                }
            }
        }
        ##compare plot is TRUE
        else{
            legend <- FALSE
            goseq <- FALSE
            suppress_details <- TRUE
            if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
                dev.new(height=10, width=20)
            }
            par(mfrow=c(1,2))
        }
        
        for(n_plot in which_network){ # This is the big loop, with n_plot being the chosen GoSeq networks
            if(goseq == TRUE){
                if(length(module_output$goseqOut) == 0){
                    stop("Goseq analysis has not been performed.")
                }
            }
            name.eid <- names(module_output$modules)[n_plot]
            eid <- module_output$modules[[n_plot]]
            network <- module_output$network
            stat <- extractScores(pvalue_annotation)
            pval <- exp(stat/(-2))
            
            if(class(network) == "graphNEL"){
                network <- igraph::graph_from_graphnel(network)
                adj_mat_network <- igraph::as_adjacency_matrix(network)
                stat.v <- stat[which(names(stat) %in%
                                         rownames(adj_mat_network))]
                stat.v <- stat.v[order(names(stat.v))]
                adj_mat_network <- adj_mat_network[order(
                    rownames(adj_mat_network)),
                    order(colnames(adj_mat_network))]
                temp_vstat <- apply(adj_mat_network, 1, function(v) return(v*stat.v)) 
                W <- (temp_vstat + t(temp_vstat))/2
                Graph_adj_mat <- igraph::graph_from_adjacency_matrix(
                    W, mode ="undirected", weighted=TRUE)
                igraph::V(Graph_adj_mat)$weight <- stat.v
                network <- Graph_adj_mat
            }
            
            vect2color <- function(v, palettev, breaks) {
                w <- v;
                for (i in 1:length(palettev)) {
                    w[which( v >= breaks[i] & v < breaks[i+1] )] = palettev[i]
                }
                return(w);
            }
            
            ## Vertex palette:
            vertexPalette.v <-
                colorRampPalette(c("white","gray85","gray65","salmon"))(50)
            vertexBreaks.v <- rev(-2*log(Hmisc::cut2(pval, g=50,
                                                     onlycuts=TRUE)))
            vertexBreaks.v[51] <- vertexBreaks.v[51]+0.001
            ## Edge palette: grey to black
            edgePalette.v <-
                colorRampPalette(c("white","gray85","gray65","salmon"))(50);
            edgeBreaks.v <- Hmisc::cut2(igraph::E(network)$weight, g=50,
                                        onlycuts=TRUE)
            
            ## Compute iGraph object
            h <- igraph::induced_subgraph(network, eid);
            stat.v <- stat[igraph::V(h)$name];
            pval.v <- pval[igraph::V(h)$name];
            pval.v[which(pval.v == 0)] <- 0.000000001
            
            par(mar=c(4, 0, 2, 0))
            ## Color edges between grey and red according to significance
            igraph::E(h)$color <- vect2color(igraph::E(h)$weight, edgePalette.v, edgeBreaks.v);
            ## Color nodes blue to yellow according to hyper/hypo-methylation
            igraph::V(h)$color <- vect2color(stat.v, vertexPalette.v, vertexBreaks.v);
            #### this is where I could Modify to allow other annotation names
            ## vl = unlist(entrez2symbol[V(h)$name])
            vl <- igraph::V(h)$name
            igraph::V(h)$color[which(1-pchisq(igraph::V(h)$weight,2) < p_thresh)] <- "red"
            igraph::E(h)$color[which(1-pchisq(igraph::E(h)$weight,4) < p_thresh)] <- "red"
            
            if(layout == "circle"){layout1 <- igraph::layout_in_circle(h)}
            if(layout == "fr"){layout1 <- igraph::layout_with_fr(h)}
            if(layout == "dh"){layout1 <- igraph::layout_with_dh(h)}
            if(layout == "kk"){layout1 <- igraph::layout_with_kk(h)}
            
            layout1_scaled <- cbind(scales::rescale(layout1[, 1], to=c(-1, 1)),
                                    scales::rescale(layout1[, 2], to=c(-1, 1)))
            
            counter <- 0
            while(counter < 2){
                plot(h, layout=layout1, ## same layout each time
                     vertex.label="",
                     vertex.frame.color="black",
                     vertex.label.dist=.1,
                     vertex.label.font=3, vertex.label.color="black",
                     vertex.size =
                         if(length(igraph::V(h)) < 50) { 15 }
                     else {15*13/length(igraph::V(h))},
                     edge.width = if(length(igraph::V(h))< 50) { 2 } else { 1 },
                     ylim=c(-1, 1.5), xlim=c(-1, 1)
                )
                
                halfCircle <- function(x, y, r, r2=.75, quarter=FALSE, start=0,
                                       end=pi, nsteps=30, col=NULL, lwd=1,
                                       border=NULL){
                    if(isTRUE(quarter))
                    {
                        
                        rs <- seq(start, end, len=nsteps)
                        xc <- x+r*cos(rs)
                        yc <- y+r*sin(rs)
                        xc2 <- x+r*r2*cos(rs)
                        yc2 <- y+r*r2*sin(rs)
                        polygon(c(xc, rev(xc2)), c(yc, rev(yc2)), col=col,
                                lwd=lwd, border=border)
                        
                        
                    }
                    if(!isTRUE(quarter))
                    {
                        rs <- seq(start, end, len=nsteps)
                        xc <- x+r*cos(rs)
                        yc <- y+r*sin(rs)
                        polygon(xc, yc, col=col, lwd=lwd, border=border)
                    }
                }
                
                arctext <- function(x, y, r, start, end, what, cex=1){
                    delta=(end-start)/4
                    rs <- seq(start, end, len=3)
                    xc <- x+r*cos(rs[2])
                    yc <- y+r*sin(rs[2])
                    text(xc, yc, what, srt=180+atan2(((y+r*sin(rs[3]))-
                                                          (y+r*sin(rs[1]))),((x+r*cos(rs[3]))-(x+r*cos(rs[1]))
                                                          ))*180/pi, cex=cex )
                }
                
                methcol <- c(meth_low_col, meth_mid_col, meth_hi_col, "white")
                names(methcol) <- c("Low", "Med", "High", "NoData")
                expcol <- c(exp_low_col, exp_mid_col, exp_hi_col, "white")
                names(expcol) <- c("Low", "Med", "High", "NoData")
                
                pval_data <- slot(slot(pvalue_annotation, "score_data"), "pval_data")
                genes_score <- slot(pvalue_annotation, "score_data")@genes
                effect_data <- slot(slot(pvalue_annotation, "score_data"), "effect_data")
                signs_idx <- slot(slot(pvalue_annotation, "score_data"), "signs_index")
                if(any(suppress_details == FALSE, counter == 1)){
                    for(i in 1:nrow(layout1_scaled)){
                        
                        halfCircle(x=layout1_scaled[i, 1], y=layout1_scaled[i, 2],
                                   r=ifelse(length(igraph::V(h)) < 50, 0.075, 0.025), start=pi/2,
                                   end=2*pi/2, quarter=TRUE, lwd=1,
                                   col=ifelse(!is.na(pval_data$expression_pvalue[which(
                                       genes_score %in% igraph::V(h)$name[i])]),
                                       expcol[ifelse(abs(pval_data$expression_pvalue[
                                           which(genes_score %in% igraph::V(h)$name[i])
                                           ]) < p_thresh, ifelse(sign(
                                               effect_data$expression_effect[which(
                                                   genes_score %in%
                                                       igraph::V(h)$name[i])]) ==
                                                   1, 3, 1), 2)], expcol[4])
                        )
                        
                        start <- pi
                        delta <- (3*pi/2)/nrow(signs_idx)
                        
                        for(j in signs_idx[, 3]){
                            
                            score_graph_col<-returnPvalueCol(slot(
                                pvalue_annotation, "score_data"),
                                j)[which(genes_score %in% igraph::V(h)$name[i])]
                            
                            halfCircle(x=layout1_scaled[i, 1],
                                       y=layout1_scaled[i, 2],
                                       r=ifelse(length(igraph::V(h)) <
                                                    50, 0.075, 0.025),
                                       start=start, end=start+delta, quarter=TRUE,
                                       col=ifelse(!is.na(score_graph_col),
                                                  methcol[ifelse(abs(score_graph_col) <
                                                                     p_thresh,
                                                                 ifelse(effect_data[, grep(j, colnames(
                                                                     effect_data))][which(genes_score %in%
                                                                                              igraph::V(h)$name[i])] == 1,
                                                                     3, 1), 2)],methcol[4])
                            )
                            start <- start+delta
                        }
                    }
                }
                
                if(any(legend == TRUE, counter == 1)){
                    if(any(suppress_details == FALSE, counter == 1)){
                        num_factors <- nrow(signs_idx)
                        halfCircle(x=-1.25, y=1.25, r=.4, start=pi/2, end=pi, quarter=TRUE)
                        start <- pi/2
                        delta <- pi/8
                        
                        for(g in 1:4){
                            halfCircle(x=-1.25, y=1.25, r=.37, r2=.89, start=start+delta*(g-1),
                                       end=start+delta*g, col=expcol[g], quarter=TRUE)
                            arctext(x=-1.25, y=1.25, r=.35, start=start+delta*(g-1),
                                    end=start+delta*g, names(expcol)[g], cex=.5)
                        }
                        text(-1.25+1.2*cos(seq(start, start+delta*2, len=30)[15]),
                             1.25+.4*.75*sin(seq(start,start+delta*2, len=30)[15]),
                             "expression", cex=.6)
                        
                        start <- pi
                        delta <- (3*pi/2)/num_factors
                        
                        for(j in signs_idx[, 3]){
                            halfCircle(x=-1.25, y=1.25, r=.4, start=start,
                                       end=start+delta, quarter=TRUE)
                            
                            for(g in 1:4){
                                start <- start
                                end <- start+delta
                                delta2 <- (end-start)/4
                                halfCircle(x=-1.25, y=1.25, r=.37, r2=.89,
                                           start=start+delta2*(g-1),
                                           end=start+delta2*g,
                                           col=methcol[g], quarter=TRUE)
                                
                                arctext(x=-1.25, y=1.25, r=.35,
                                        start=start+delta2*(g-1),
                                        end=start+delta2*g,
                                        ifelse(num_factors<=4,
                                               names(methcol)[g],
                                               substring(names(methcol)[g],
                                                         1, 1)),
                                        cex=.5)
                            }
                            
                            text(-1.25+.57*cos(seq(start, start+delta, len=30)[15]),
                                 1.25+.57*.75*sin(seq(start, start+delta, len=30)[15]),
                                 paste(strsplit(j, "_")[[1]], collapse="\n"), cex=.6)
                            start <- start+delta
                        }
                    }
                    halfCircle(x=-1.25, y=1.25, r=.3, end=2*pi, col="white")
                    
                    vertexPalette.v[which(
                        vertexBreaks.v >= min(igraph::V(h)$weight[which(
                            igraph::V(h)$weight >= qchisq(1-p_thresh,2))])) - 1] <- "red"
                    
                    points(seq(-1.45, -1.05, length.out=50), rep(1.35, 50),
                           col=vertexPalette.v, pch=15, cex=2.5)
                    
                    text(-1.20, 1.44, expression('node             ',Chi[2]^2))
                    text(-1.45, 1.295, round(vertexBreaks.v[1], 2))
                    text(-1.25, 1.295, round(vertexBreaks.v[26], 2))
                    text(-1.07, 1.295, paste(">",round(qchisq(1-p_thresh,2), 2),sep=""))
                    
                    if(length(which(
                        igraph::E(h)$weight >= qchisq(1-p_thresh, 4))) > 0){
                        
                        edgePalette.v[which(
                            edgeBreaks.v >= min(igraph::E(h)$weight[which(
                                igraph::E(h)$weight >= qchisq(1-p_thresh, 4))]))-1] <- "red"
                    }
                    points(seq(-1.45, -1.05, length.out=50), rep(1.15, 50),
                           col=edgePalette.v, pch=15, cex=2.5)
                    text(-1.20, 1.24, expression('edge             ',Chi[4]^2))
                    text(-1.45, 1.095, round(edgeBreaks.v[1], 2))
                    text(-1.25, 1.095, round(edgeBreaks.v[26], 2))
                    text(-1.07, 1.095, paste(">", round(qchisq(1-p_thresh,4), 2), sep=""))
                }
                
                addShadowText(layout1_scaled[, 1], layout1_scaled[, 2], vl, font=2,
                              cex=if(label_scale == TRUE){
                                  scales::rescale(stat.v, to=(c(.5, 2)))
                              }
                              else{.5},
                              bg="white", col="black")
                
                if(namestyle == "refseq"){
                    ref2eg <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egREFSEQ2EG)
                    eg2sym <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
                    
                    text(layout1_scaled[, 1], layout1_scaled[, 2]-.05,
                         sapply(vl,
                                function(k){
                                    ifelse(is.null(ref2eg[[k]]), return(NA),
                                           return(eg2sym[[ref2eg[[k]]]]))
                                })
                    )
                    
                    text(0, 1.65, paste("Network built around", name.eid,
                                        ifelse(is.null(ref2eg[[name.eid]]), NA,
                                               eg2sym[[ref2eg[[name.eid]]]])))
                }
                if(namestyle == "symbol"){
                    text(0, 1.7, paste(
                        "Network built around", name.eid, "\nChi-square P-value=",
                        round(module_output$moduleStats[[n_plot]][2],4)))
                }
                
                if(goseq == TRUE){
                    if(nrow(module_output$goseqOut[[n_plot]]) > 0){
                        
                        text("Num\nGenes", x=1.35, y=1.6, font=2)
                        text("Enriched\nPathway/Term", x=2, y=1.6, font=2)
                        for(i in 1:nrow(module_output$goseqOut[[n_plot]])){
                            text(module_output$goseqOut[[n_plot]][i, 4], x=1.3,
                                 y=seq(1.4, -1, length.out=nrow(module_output$goseqOut[[n_plot]]))[i],
                                 adj = c(0, 0))
                            
                            text(module_output$goseqOut[[n_plot]][i, 6], x=1.5,
                                 y=seq(1.4, -1, length.out=nrow(module_output$goseqOut[[n_plot]]))[i],
                                 adj= c(0, 0))
                        }
                        
                    }
                    else {
                        text("No enriched terms from\nGoseq", x=1.5, y=1.6, font=2)
                        
                    }
                }
                
                if(counter == 1){break}
                counter <- 2
            }
            
            if(is.null(pdf_out)){
                if(n_plot != which_network[length(which_network)]){
                    message("Press key to go to next plot")
                    readline()
                }
            }
        }
        if(!is.null(pdf_out)){dev.off()}
    }
)



setMethod(
    f="plotDensityPval",
    signature="PvalueAnnotation",
    definition=function(pvalue_annotation, ref="expression_pvalue", ...)
    {
        palette(c("red", "green", "blue", "gold", "orange",
                        "purple", "magenta"))
        pval_data <- slot(slot(pvalue_annotation, "score_data"), "pval_data")
        pval_col_names <- colnames(pval_data)

        if(!any(grepl(ref, pval_col_names))){
            stop("paste(Reference is not one of the available:",
                 pval_col_names)
        }
        if(length(grep(ref, pval_col_names)) > 1){
            stop("Reference was not specific enough.")
        }

        ref <- pval_col_names[grep(ref, pval_col_names)]
        nonref <- pval_col_names[which(!grepl(ref, pval_col_names))]
        dens_range <- lapply(pval_data,
                             function(i){
                                dens <- density(as.numeric(i), na.rm=TRUE)
                                x <- dens$x
                                y <- dens$y
                                return(as.data.frame(cbind(x,y)))
                             }
                     )
        dens_range <- do.call(rbind, dens_range)
        x_range <- range(dens_range[,1])
        y_range <- range(dens_range[,2])
        y_range[2] <- y_range[2]+diff(y_range)*.1

        if(!names(dev.cur()) %in% c("RStudioGD")){
            par(mar=c(4,4,4,10), xpd=FALSE)
        }
        else{
            par(mar=c(4,4,4,4))
        }
        plot.new()
        plot.window(xlim=x_range, ylim=y_range)
        axis(1)
        axis(2)
        box()
        title("Density of P-values/Scores", xlab="P-values",
              ylab="Density")

        score_sign_idx <- slot(slot(pvalue_annotation, "score_data"), "signs_index")
        sapply(score_sign_idx[, 3], function(i){
            message(paste("Plotting: ", i))
            if(!all(is.na(pval_data[[grep(i, pval_col_names)]]))){
                lines(density(abs(as.numeric(pval_data[[grep(i, pval_col_names)]])),
                              na.rm=TRUE), col= which(score_sign_idx[, 3] == i))
            }
        })

        lines(density(abs(
            returnPvalueCol(slot(pvalue_annotation, "score_data"), ref)
            ), na.rm=TRUE), col="black")

        legend(x_range[2]+.01, y_range[2]+.01, gsub("_", "\n", c(paste(
            strsplit(ref, "_pvalue")[[1]], "(REF)"), score_sign_idx[, 3])),
            fill=c("black", 1:length(score_sign_idx[, 3])), xpd=TRUE,bty="n")
    }
)
