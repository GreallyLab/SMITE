##Classes

setClass(
    Class="PvalueObject",
    representation=representation(
        pval_data="data.frame",
        effect_data="data.frame",
        genes="character",
        signs_index="data.frame",
        scores="data.frame",
        trans="numeric",
        scoring_vector="numeric",
        module_output="list")
)

setClass(
    Class="PvalueAnnotation",
    representation=representation(
        annotation="GRangesList",
        modifications="GRangesList",
        expression="ExpressionSet",
        score_data="PvalueObject"))

