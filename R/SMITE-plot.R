
setMethod(
    f="plotCompareScores", 
    signature="PvalueAnnotation", 
    definition=function(object, x_name, y_name, ...){
        if(all(!grepl(x_name,colnames(slot(slot(object, "score_data")
                                          , "pval_data")))))
        {
            stop(paste("Provided x_name was not one of the following:",
                    paste(colnames(slot(slot(object, "score_data")
                                        , "pval_data")), collapse=",")))
        }
              
        if(all(!grepl(y_name,colnames(slot(slot(object, "score_data")
                                          , "pval_data")))))
        {
            stop(paste("Provided y_name was not one of the following:",
                    paste(colnames(slot(slot(object, "score_data")
                                        , "pval_data")),collapse=",")))
        }
            
        if(length(grep(x_name,colnames(slot(slot(object, "score_data")
                                                , "pval_data")))) > 1)
        {
            stop(paste("Provided x_name was not unique"))
        }
            
        if(length(grep(y_name,colnames(slot(slot(object, "score_data")
                                                , "pval_data")))) > 1)
        {
            stop(paste("Provided y_name was not unique"))
        }   
        
        x <- (-log(abs(as.numeric(
            slot(
                slot(object, "score_data"),"pval_data")[[grep(x_name, colnames(
                    slot(
                        slot(object, "score_data"), "pval_data")))]])))*sign(
                            slot(
                                slot(object, "score_data"), "effect_data")[[grep(
                                    x_name, colnames(slot(slot(object, "score_data"), 
                                                          "effect_data")))]]
                            )
            )
        y <- (-log(abs(as.numeric(slot(slot(object, "score_data"),
                    "pval_data")[[grep(y_name, colnames(slot(slot(object, 
                    "score_data"), "pval_data")))]]))
                    )*sign(slot(slot(object, "score_data"), 
                    "effect_data")[[grep(y_name, colnames(slot(slot(object, 
                                            "score_data"), "effect_data")))]])
              )
        ggplot(data=data.frame(x=x, y=y), 
            aes(x=x, y=y))+stat_binhex(bins=50)+geom_hline(yintercept=0, 
            colour="red", linetype = "longdash")+geom_vline(xintercept=0,
            colour="red", linetype = "longdash") +labs(title = paste(y_name, 
            "vs", x_name, "p-values comparing effect direction"), 
            x=paste(x_name, "Score * Effect Direction"), y=paste( y_name, 
            "Score * Effect Direction")
        )
        
    }
)
		  
		  
setMethod(
    f="SMITEplotModule", 
    signature="PvalueAnnotation", 
    definition=function(object, p_thresh=0.05, which.network=1, goseq=FALSE, 
                        layout="fr", legend=TRUE, namestyle="symbol",
                        suppressDetails=FALSE, meth_hi_col="blue", 
                        meth_low_col="yellow1", meth_mid_col="gray90", 
                        exp_hi_col="red1", exp_low_col="chartreuse1", 
                        exp_mid_col="gray90", label_scale=TRUE,
                        comparePlot=FALSE,pdfOut=NULL)
    {
        
        if(!is.null(pdfOut)){
            unlink(pdfOut)
            if(comparePlot==TRUE){
                pdf(pdfOut, width=24, height=8.5)
            }
            else {
                pdf(pdfOut, width=16, height=8.5)
            }
        }
        
        for(n_plot in which.network){
            if(goseq == TRUE){
                if(length(slot(slot(object, "score_data"), 
                           "module_output")$goseqOut) == 0){ 
                stop("Goseq analysis has not been performed.")
                }
            }
            name.eid <- 
                names(slot(slot(object, "score_data"), 
                       "module_output")$modules)[n_plot]
            eid <- slot(slot(object, "score_data"), 
                 "module_output")$modules[[n_plot]]
            g <- slot(slot(object, "score_data"), "module_output")$network
            stat <- SMITEextractScores(object)
            pval <- exp(stat/(-2))
        
            if(class(g)=="graphNEL"){
                g <- graph_from_graphnel(g)
                A <- as_adjacency_matrix(g)
                stat.v <- stat[which(names(stat) %in% rownames(A))]
                stat.v <- stat.v[order(names(stat.v))]
                A <- A[order(rownames(A)),order(colnames(A))]
                temp1 <- apply(A, 1, function(v) return(v*stat.v))
                W <- (temp1 + t(temp1))/2
                G <- graph_from_adjacency_matrix(W, mode = "undirected", weighted=TRUE)
                V(G)$weight <- stat.v
                g <- G
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
            vertexBreaks.v <- rev(-2*log(cut2(pval, g=50, onlycuts=TRUE)))
            vertexBreaks.v[51] <- vertexBreaks.v[51]+0.001
            ## Edge palette: grey to black 
            edgePalette.v <- 
                colorRampPalette(c("white","gray85","gray65","salmon"))(50);
            edgeBreaks.v <- cut2(E(g)$weight, g=50, onlycuts=TRUE)         
            
            ## Compute iGraph object
            h <- induced_subgraph(g, eid);
            stat.v <- stat[V(h)$name];
            pval.v <- pval[V(h)$name];
            pval.v[which(pval.v == 0)] <- 0.000000001
            
            par(mar=c(4, 0, 2, 0))
            ## Color edges between grey and red according to significance
            E(h)$color <- vect2color(E(h)$weight, edgePalette.v, edgeBreaks.v);
            ## Color nodes blue to yellow according to hyper/hypo-methylation
            V(h)$color <- vect2color(stat.v, vertexPalette.v, vertexBreaks.v);
            #### this is where I could Modify to allow other annotation names
            ## vl = unlist(entrez2symbol[V(h)$name])
            vl <- V(h)$name
            V(h)$color[which(1-pchisq(V(h)$weight,2) < p_thresh)] <- "red"
            E(h)$color[which(1-pchisq(E(h)$weight,4) < p_thresh)] <- "red"
            
            if(layout == "circle"){layout1 <- layout_in_circle(h)}
            if(layout == "fr"){layout1 <- layout_with_fr(h)}
            ## if two side by side plots are not needed
            if(comparePlot == FALSE){
                ## no pdf
                if(is.null(pdfOut)){
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
                suppressDetails <- TRUE
                if(!names(dev.cur()) %in% c("RStudioGD","pdf")){
                    dev.new(height=10, width=20)
                }
                par(mfrow=c(1,2))
            }
            
            layout1_scaled <- cbind(rescale(layout1[, 1], to=c(-1, 1)),
                                    rescale(layout1[, 2], to=c(-1, 1)))
            
            counter <- 0
            while(counter < 2){
                plot(h, layout=layout1, ## same layout each time
                     vertex.label="", 
                     vertex.frame.color="black", 
                     vertex.label.dist=.1, 
                     vertex.label.font=3, vertex.label.color="black", 
                     ##vertex.size = 15*13/length(V(h)), 
                     ##edge.width = 160/length(V(h))
                     vertex.size = 
                         if(length(V(h)) < 50) { 15 } else {15*13/length(V(h))}, 
                     edge.width = if(length(V(h))< 50) { 2 } else { 1 }, 
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
                
                if(any(suppressDetails == FALSE, counter == 1)){
                    for(i in 1:nrow(layout1_scaled)){
                        
                        halfCircle(x=layout1_scaled[i, 1], y=layout1_scaled[i, 2], 
                            r=ifelse(length(V(h)) < 50, 0.075, 0.025), start=pi/2, 
                            end=2*pi/2, quarter=TRUE, lwd=1, 
                            col=ifelse(!is.na(slot(slot(object, "score_data"),
                            "pval_data")$expression_pvalue[which(slot(object,
                            "score_data")@genes%in%V(h)$name[i])]), 
                            expcol[ifelse(abs(slot(slot(object, "score_data"), 
                            "pval_data")$expression_pvalue[which(slot(object, 
                            "score_data")@genes%in%V(h)$name[i])]) < p_thresh, 
                            ifelse(sign(slot(slot(object, "score_data"), 
                            "effect_data")$expression_effect[
                            which(slot(object,"score_data")@genes %in% V(h)$name[i])])
                            == 1, 3, 1), 2)], expcol[4])
                        )
                        
                        start <- pi
                        delta <- (3*pi/2)/nrow(slot(slot(object, "score_data"), 
                                                 "signsindex"))
                        
                        for(j in slot(slot(object, "score_data"), "signsindex")[, 3]){
                            
                            halfCircle(x=layout1_scaled[i, 1], y=layout1_scaled[i, 2], 
                                       r=ifelse(length(V(h))< 50, 0.075, 0.025), 
                                       start=start, end=start+delta, quarter=TRUE, 
                                       
                                      col=ifelse(!is.na(returnPvalueCol(slot(object, 
                                      "score_data"), j)[which(slot(object, 
                                        "score_data")@genes%in%V(h)$name[i])]), 
                                {
                                methcol[ ifelse(abs(returnPvalueCol(slot(object, "score_data"), j)[
                                which(slot(object, "score_data")@genes%in%V(h)$name[i])])<p_thresh, 
                                ifelse(sign(slot(slot(object, "score_data"), 
                                "effect_data")[, grep(j, colnames(slot(slot(object, 
                                "score_data"), "effect_data")))][which(slot(object, 
                                "score_data")@genes%in%V(h)$name[i])]) == 1, 3, 1), 2)]
                                }, methcol[4])
                            )
                            start <- start+delta
                        }    
                    }	
                }
    
                if(any(legend == TRUE, counter == 1)){
                    if(any(suppressDetails == FALSE, counter == 1)){
                        num_factors <- nrow(slot(slot(object, "score_data"), "signsindex"))
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
                
                        for(j in slot(slot(object, "score_data"), "signsindex")[, 3]){
                            halfCircle(x=-1.25, y=1.25, r=.4, start=start, end=start+delta,
                                       quarter=TRUE)
                            for(g in 1:4){
                                start <- start
                                end <- start+delta
                                delta2 <- (end-start)/4
                                halfCircle(x=-1.25, y=1.25, r=.37, r2=.89, 
                                           start=start+delta2*(g-1), end=start+delta2*g, 
                                           col=methcol[g], quarter=TRUE)
                        
                                arctext(x=-1.25, y=1.25, r=.35, start=start+delta2*(g-1), 
                                        end=start+delta2*g, ifelse(num_factors<=4, 
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
                        vertexBreaks.v >= min(V(h)$weight[which(
                            V(h)$weight >= qchisq(1-p_thresh,2))]))-1] <- "red"
                    
                    points(seq(-1.45, -1.05, length.out=50), rep(1.35, 50), 
                           col=vertexPalette.v, pch=15, cex=2.5)
        
                    text(-1.20, 1.44, expression('node             ',Chi[2]^2))
                    text(-1.45, 1.295, round(vertexBreaks.v[1], 2))
                    text(-1.25, 1.295, round(vertexBreaks.v[26], 2))
                    text(-1.07, 1.295, paste(">",round(qchisq(1-p_thresh,2), 2),sep=""))
        
                    if(length(which(
                        E(h)$weight >= qchisq(1-p_thresh,4))) > 0){
                        
                        edgePalette.v[which(
                            edgeBreaks.v >= min(E(h)$weight[which(
                                E(h)$weight >= qchisq(1-p_thresh,4))]))-1] <- "red"
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
                                  rescale(stat.v, to=(c(.5, 2)))
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
                        round(slot(slot(object, "score_data"), 
                                   "module_output")$moduleStats[[n_plot]][2],4)))
                }
    
                if(goseq == TRUE){
                    if(nrow(slot(slot(object, "score_data"),
                                 "module_output")$goseqOut[[n_plot]]) > 0){
            
                        text("Num\nGenes", x=1.35, y=1.6, font=2)
                        text("Enriched\nPathway/Term", x=2, y=1.6, font=2)
                        for(i in 1:nrow(slot(slot(object, "score_data"), 
                                             "module_output")$goseqOut[[n_plot]])){
                            text(slot(slot(object, "score_data"), 
                                      "module_output")$goseqOut[[n_plot]][i, 4], x=1.3,
                                 y=seq(1.4, -1, length.out= nrow(slot(
                                     slot(object, "score_data"),
                                     "module_output")$goseqOut[[n_plot]]))[i], 
                                 adj = c(0, 0))
                        
                            text(slot(slot(object, "score_data"),
                                      "module_output")$goseqOut[[n_plot]][i, 6], x=1.5, 
                                 y=seq(1.4, -1, length.out= nrow(slot(
                                     slot(object, "score_data"),
                                     "module_output")$goseqOut[[n_plot]]))[i], 
                                 adj= c(0, 0))
                        }
            
                    } 
                    else {
                        text("No enriched terms from\nGoseq", x=1.5, y=1.6, font=2)
            
                    }
                }
                
                if(counter == 1){break}
                counter <- 2
                if(comparePlot == TRUE){counter <- 1}
            }
            
            if(is.null(pdfOut)){
                if(n_plot != which.network[length(which.network)]){
                    message("Press key to go to next plot")
                    readline()
                }
            }
        }
        if(!is.null(pdfOut)){dev.off()}
    }
)
	
	
	
setMethod(
    f="plotDensityPval", 
    signature="PvalueAnnotation", 
    definition=function(x, ref="expression_pvalue", ...){
        palette(c("red", "green", "blue", "gold", "orange", 
                        "purple", "magenta"))
           
        if(!any(grepl(ref,colnames(slot(slot(x, "score_data"), 
                                        "pval_data"))))){
            stop("paste(Reference is not one of the available:", 
                 colnames(slot(slot(x, "score_data"), "pval_data")))
        }
        if(length(grep(ref,colnames(slot(slot(x, "score_data"), 
                                        "pval_data"))))>1){
            stop("Reference was not specific enough.")
        }
              
        ref <- colnames(slot(slot(x, "score_data"), "pval_data"))[
            grep(ref,colnames(slot(slot(x, "score_data"), "pval_data")))]
        nonref <- colnames(slot(slot(x, "score_data"), "pval_data"))[
            which(!grepl(ref,colnames(slot(slot(x, "score_data"), "pval_data"))))]
        dens_range <- lapply(slot(slot(x,"score_data"), "pval_data"), 
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
              
        sapply(slot(slot(x, "score_data"), "signsindex")[, 3], 
                    function(i){
                        message(paste("Plotting: ",i))
                        if(!all(is.na(slot(slot(x, "score_data"), "pval_data")[[
                            grep(i, colnames(slot(slot(x, "score_data"), "pval_data")
                                             ))]]
                        ))){ 
                            
                            lines(density(abs(as.numeric(slot(
                                slot(x, "score_data"), "pval_data")
                                [[grep(i, colnames(slot(
                                    slot(x, "score_data"),"pval_data")))]])),
                                na.rm=TRUE), col= which(slot(
                                    slot(x, "score_data"),
                                    "signsindex")[, 3] == i))
                        }
                    }
               )
        
        lines(density(abs(
            returnPvalueCol(slot(x, "score_data"), ref)
            ), na.rm=TRUE), col="black")
              
        legend(x_range[2]+.01,y_range[2]+.01, gsub("_","\n",c(paste(
                strsplit(ref,"_pvalue")[[1]], "(REF)"),slot(slot(x, 
                "score_data"), "signsindex")[, 3])), fill=c("black", 
                1:length(slot(slot(x, "score_data"), "signsindex")[, 3])),
                xpd=TRUE,bty="n")         
    }
)
