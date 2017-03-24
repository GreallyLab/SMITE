setGeneric(
    name="stoufferTest",
    def=function(pvalues, weights)
    {
        standardGeneric("stoufferTest")
    }
)

setGeneric(
    name="annotationOutput",
    def=function(pvalue_annotation)
    {
        standardGeneric("annotationOutput")
    }
)

setGeneric(
    name="returnPvalueCol",
    def=function(pval_object, col_name)
    {
        standardGeneric("returnPvalueCol")
    }
)

setGeneric(
    name="makePvalueAnnotation",
    def=function(data, other_data=NULL, other_tss_distance=10000,
                 promoter_upstream_distance=1000, promoter_downstream_distance=1000,
                 strand_col=NULL, gene_name_col=NULL)
    {
        standardGeneric("makePvalueAnnotation")
    }
)

setGeneric(
    name="convertGeneIds",
    def=function(gene_IDs, ID_type, ID_convert_to, delim=NULL, verbose=FALSE)
    {
        standardGeneric("convertGeneIds")
    }
)


setGeneric(
    name="annotateExpression",
    def=function(pvalue_annotation, expr_data, effect_col=NULL, pval_col=NULL)
    {
        standardGeneric("annotateExpression")
    }
)


setGeneric(
    name="annotateModification",
    def=function(pvalue_annotation, mod_data, weight_by=NULL,
                 weight_by_method="Stouffer", mod_included=NULL, mod_corr=TRUE,
                 mod_type="methylation", verbose=FALSE)
    {
        standardGeneric("annotateModification")
    }
)


setGeneric(
    name="removeModification",
    def=function(pvalue_annotation, mod_type="methylation")
    {
        standardGeneric("removeModification")
    }
)


setGeneric(
    name="makePvalueObject",
    def=function(pvalue_annotation, effect_directions=NULL)
    {
        standardGeneric("makePvalueObject")
    }
)


setGeneric(
    name="plotDensityPval",
    def=function(pvalue_annotation, ref="expression_pvalue", ...)
    {
        standardGeneric("plotDensityPval")
    }
)

setGeneric(
    name="normalizePval",
    def=function(pvalue_annotation, trans, ref="expression_pvalue", method="rescale")
    {
        standardGeneric("normalizePval")
    }
)


setGeneric(
    name="scorePval",
    def=function(pvalue_annotation, weights)
    {
        standardGeneric("scorePval")
    }
)

setGeneric(
    name="runSpinglass",
    def=function(pvalue_annotation, network, random_alpha = 0.05, gam = 0.5,
                 node_alpha = 0.05, maxsize = 500, minsize = 8, num_iterations = 1000,
                 simplify=TRUE)
    {
        standardGeneric("runSpinglass")
    }
)



setGeneric(
    name="runBioNet",
    def=function(pvalue_annotation, network, alpha = 0.05)
    {
        standardGeneric("runBioNet")
    }
)

setGeneric(
    name="plotModule",
    def=function(pvalue_annotation, p_thresh=0.05, which_network=1, goseq=FALSE,
                 layout="fr",legend=TRUE, namestyle="symbol",
                 suppress_details=FALSE,
                 meth_hi_col="blue", meth_low_col="yellow1",
                 meth_mid_col="gray90", exp_hi_col="red1",
                 exp_low_col="chartreuse1", exp_mid_col="gray90",
                 label_scale=TRUE,compare_plot=FALSE, pdf_out=NULL)
    {
        standardGeneric("plotModule")
    }
)


setGeneric(
    name="runGOseq",
    def=function(pvalue_annotation, p_thresh=0.05, coverage, type="reactome")
    {
        standardGeneric("runGOseq")
    }
)


setGeneric(
    name="searchGOseq",
    def=function(pvalue_annotation, search_string, wholeword=FALSE)
    {
        standardGeneric("searchGOseq")
    }
)



setGeneric(
    name="extractGOseq",
    def=function(pvalue_annotation, which_network=NULL)
    {
        standardGeneric("extractGOseq")
    }
)


setGeneric(
    name="plotCompareScores",
    def=function(pvalue_annotation, x_name, y_name, ...)
    {
        standardGeneric("plotCompareScores")
    }
)


setGeneric(
    name="extractModification",
    def=function(pvalue_annotation, mod_type="methylation")
    {
        standardGeneric("extractModification")
    }
)

setGeneric(
    name="extractExpression",
    def=function(pvalue_annotation)
    {
        standardGeneric("extractExpression")
    }
)

setGeneric(
    name="extractModSummary",
    def=function(pvalue_annotation)
    {
        standardGeneric("extractModSummary")
    }
)

setGeneric(
    name="extractScores",
    def=function(pvalue_annotation)
    {
        standardGeneric("extractScores")
    }
)

setGeneric(
    name="addShadowText",
    def=function(x, y=NULL, labels, col='white', bg='black',
                 theta=seq(pi/4, 2*pi, length.out=8), r=0.1, ...)
    {
        standardGeneric("addShadowText")
    }
)


setGeneric(
    name="extractModules",
    def=function(pvalue_annotation, which_module=NULL)
    {
        standardGeneric("extractModules")
    }
)

setGeneric(
    name="highScores",
    def=function(pvalue_annotation, alpha=0.05)
    {
        standardGeneric("highScores")
    }
)



