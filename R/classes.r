combineSpectraParam = setClass("combineSpectraParam",
                               slots=c(
                                 fcol="character",
                                 mzd="numeric",
                                 intensityFun="ANY"
                               ),
                               prototype=prototype(
                                 fcol="peak_id",
                                 mzd=0.005,
                                 intensityFun=mean
                               ))
