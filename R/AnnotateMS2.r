#' MS/MS Spectra Extraction
#'
#' @description 
#' @param
#' @export
#' @seealso 
#' 
Dulce_MS2SpectraExtract = function(data, features, annotations){
  
  if(any(table(data@featureData@data$msLevel, data@featureData@data$fileIdx)[2,]==0)){
    stop("Dulce error: Some of the files do not have MS level 2 data.")
  }
  
  data = as(data, "XCMSnExp")
  
  chromPeaks(data) = chromPeaks(features)
  featureDefinition = featureDefinitions(features)
  featureDefinition_filtered = featureDefinition[paste0(featureDefinition$mzmed,"_", featureDefinition$rtmed) %in%
                                                   paste0(annotations$i.mz,"_",annotations$rt),]
  
  featureDefinitions(data) = featureDefinition_filtered
  
  MS2Spectra = featureSpectra(data, msLevel = 2, expandMz = 0.005)
  
  return(MS2Spectra)
}

# This should be in a different script file.
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

#' MS/MS Spectra Combination and Normalization
#'
#' @description 
#' @param
#' @export
#' @seealso 
#' 
Dulce_MS2CombNorm = function(data, csp=NULL, norm.method="max"){
  
  if (class(data)!="MSpectra"){stop("Dulce error: 'data' object not from 'MSpectra' class.")}
  if (is.null(csp)){csp = combineSpectraParam()}
  
  data@listData = lapply(data, clean, all=T)
  data = combineSpectra(data, fcol=csp@fcol, mzd=csp@mzd, intensityFun=csp@intensityFun)
  data@listData = lapply(data, MSnbase::normalise, method=norm.method)
  
  return(data)
}

#' MS/MS Precursor Matching
#'
#' @description 
#' @param
#' @export
#' @seealso 
#' 
Dulce_MS2PrecursorMatch = function(data, annotations){
  
  precursorsMz = precursorMz(data)
  precursorsRt = rtime(data)
  
  precursorsData = data.frame(peakId = names(precursorsMz),
                              precursorMz = precursorsMz,
                              rtime = precursorsRt,
                              mzmin = precursorsMz - 0.0025,
                              mzmax = precursorsMz + 0.0025)
  
  setDT(precursorsData)
  setDT(annotations)
  setkey(precursorsData, mzmin, mzmax)
  
  data = foverlaps(annotations, precursorsData) %>% drop_na()
  
  return(data)
}

#' MS/MS Distances Calculation
#'
#' @description 
#' @param
#' @export
#' @seealso 
#' 
Dulce_calculateMZDistances = function(data, annotations, MSpectra){
  
  distances = lapply(MSpectra@listData, function(x) c(dist(x@mz)))
  distances = lapply(distances, function(x) data.frame("mzDifference"=x)) %>% bind_rows(.id="peakId")
  
  # Discuss this with Margot (should be right_join, if left we are adding a 
  # lot of not annotated differences.)
  distances = dplyr::right_join(distances, annotations, by="peakId") %>%
    dplyr::select(peakId, mzDifference, precursorMz, rtime, name, ion) %>%
    dplyr::distinct()
  
  distances$sample = distances$peakId %>%
    gsub(".*\\.F|\\.S\\d+$", "", .) %>% as.numeric()
  
  distances$sample = pData(data)$name[distances$sample]
  
  return(distances)
}

#' MS/MS Distances Filtration and Highlighting
#'
#' @description 
#' @param
#' @export
#' @seealso 
#' 
Dulce_filterMZDistances = function(data, sigdiff=NULL, n=NULL, mzdiff=NULL){
  
  sigdiff.options = c(73.089, 116.131, 162.053, 132.042)
  names(sigdiff.options) = c("procA-R1", "procA-R2", "hexose", "pentose")
  
  if(is.null(sigdiff)){
    sigdiff.used = sigdiff.options
  } else{
    
    if(!any(sigdiff %in% names(sigdiff.options))){
      stop("Dulce error: significant differences you requested are not options!")}
    
    sigdiff.used = sigdiff.options[names(sigdiff.options) %in% sigdiff]
  }
  
  if(is.null(n)){
    n = length(sigdiff.used)
  }
  if(n > length(sigdiff.used)){
    stop("Dulce error: cannot ask for higher n than possible options!")
  }
  if(is.null(mzdiff)){
    mzdiff = 0.001
  }
  
  data.filt.list = list()
  
  for(i in 1:length(sigdiff.used)){
    data.filt.list[[i]] = data %>%
      filter(between(mzDifference,
                     sigdiff.used[i]-mzdiff, sigdiff.used[i]+mzdiff))
  }
  
  names(data.filt.list) = names(sigdiff.used)
  data.filtered = rbindlist(data.filt.list, idcol = "MSMS.fragment")
  
  
  #summarise number of times each sig diff found in each spectrum
  data.filtered.sum <- data.filtered %>%
    dplyr::group_by(sample, peakId, MSMS.fragment) %>%
    dplyr::summarise(n = n())
  
  #pivot wider
  data.filtered.sum.w <- data.filtered.sum %>%
    pivot_wider(names_from = MSMS.fragment, values_from = n)
  
  #count how many different sig diff per spectrum
  data.filtered.sum = data.filtered.sum %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(sample, peakId) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= !!n) %>% 
    dplyr::left_join(data.filtered.sum.w, by=c("sample","peakId"))
  
  #filter the original data to only contain these spectra
  spectra.filtered = data %>% 
    dplyr::select(-mzDifference) %>%
    dplyr::distinct(peakId, .keep_all = T) %>% 
    dplyr::right_join(data.filtered.sum, by=c("sample","peakId"))
  
  return(list(data.filtered, spectra.filtered))
}





