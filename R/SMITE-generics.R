setGeneric(
    name="stoufferTest", 
    def=function(p, w)
    {
        standardGeneric("stoufferTest")
    }
)

setGeneric(
    name="annotationOutput", 
    def=function(object){standardGeneric("annotationOutput")}
)

setGeneric(
    name="returnPvalueCol", 
    def=function(pval_object, c_name){standardGeneric("returnPvalueCol")}
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
    def=function(a, a_type, b_type, delim=NULL, verbose=FALSE)
    {standardGeneric("convertGeneIds")}
)


setGeneric(
    name="annotateExpression", 
    def=function(annotation, expdata, effect_col=NULL, pval_col=NULL){
        standardGeneric("annotateExpression")}
)


setGeneric(
    name="annotateModification", 
    def=function(annotation, modData, weight_by=NULL, 
                 weight_by_method="Stouffer", modInclude=NULL, modCorr=TRUE, 
                 modType="methylation", verbose=FALSE){
        standardGeneric("annotateModification")}
)


setGeneric(
    name="removeModification", 
    def=function(annotation, modType="methylation"){
        standardGeneric("removeModification")}
)


setGeneric(
    name="makePvalueObject", 
    def=function(object, effect_directions=NULL){
        standardGeneric("makePvalueObject")}
)


setGeneric(
    name="plotDensityPval", 
    def=function(x, ref="expression_pvalue", ...){standardGeneric(
        "plotDensityPval")}
)

setGeneric(
    name="normalizePval", 
    def=function(object, trans, ref="expression_pvalue", method="rescale"){
        standardGeneric("normalizePval")}
)


setGeneric(
    name="SMITEscorePval", 
    def=function(object, weights){standardGeneric("SMITEscorePval")}
)

setGeneric(
    name="SMITErunSpinglass", 
    def=function(object, network, random_alpha = 0.05, gam = 0.5, 
                 node_alpha = 0.05, maxsize = 500, minsize = 8, niter = 1000,
                 simplify=TRUE){
        standardGeneric("SMITErunSpinglass")
    }
)



setGeneric(
    name="SMITErunBioNet", 
    def=function(object, network, alpha = 0.05){
        standardGeneric("SMITErunBioNet")
    }
)

setGeneric(
    name="SMITEplotModule", 
    def=function(object, p_thresh=0.05, which.network=1, goseq=FALSE, 
                 layout="fr",legend=TRUE, namestyle="symbol", 
                 suppressDetails=FALSE, 
                 meth_hi_col="blue", meth_low_col="yellow1", 
                 meth_mid_col="gray90", exp_hi_col="red1", 
                 exp_low_col="chartreuse1", exp_mid_col="gray90", 
                 label_scale=TRUE,comparePlot=FALSE, pdfOut=NULL
         ){
        standardGeneric("SMITEplotModule")}
)


setGeneric(
    name="SMITErunGOseq", 
    def=function(object, p_thresh=0.05, coverage, type="reactome"){
        standardGeneric("SMITErunGOseq")}
)


setGeneric(
    name="SMITEsearchGOseq", 
    def=function(object, searchstring, wholeword=FALSE){
        standardGeneric("SMITEsearchGOseq")}
)



setGeneric(
    name="SMITEextractGOseq", 
    def=function(object, which.network=NULL){
        standardGeneric("SMITEextractGOseq")}
)


setGeneric(
    name="plotCompareScores", 
    def=function(object, x_name, y_name, ...){
        standardGeneric("plotCompareScores")}
)


setGeneric(
    name="extractModification", 
    def=function(annotation, modType="methylation"){
        standardGeneric("extractModification")}
)

setGeneric(
    name="extractExpression", 
    def=function(annotation){
        standardGeneric("extractExpression")}
)

setGeneric(
    name="extractModSummary", 
    def=function(annotation){
        standardGeneric("extractModSummary")}
)

setGeneric(
    name="SMITEextractScores", 
    def=function(annotation){
        standardGeneric("SMITEextractScores")}
)

setGeneric(
    name="addShadowText", 
    def=function(x, y=NULL, labels, col='white', bg='black',
                 theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ...){
        standardGeneric("addShadowText")}
)


setGeneric(
    name="SMITEextractModules", 
    def=function(annotation,whichModule=NULL){
        standardGeneric("SMITEextractModules")}
)

setGeneric(
    name="SMITEhighScores", 
    def=function(annotation,alpha=0.05){
        standardGeneric("SMITEhighScores")}
)



