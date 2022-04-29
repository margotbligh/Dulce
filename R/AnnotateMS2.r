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
                                                   paste0(annotations$mz,"_",annotations$rt),]
  
  featureDefinitions(data) = featureDefinition_filtered
  
  MS2Spectra = featureSpectra(data, msLevel = 2, expandMz = 0.005)
  
  return(MS2Spectra)
}


#' MS/MS Spectra Combination and Normalization
#' @include classes.r
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
  
  data = foverlaps(annotations, precursorsData) %>% 
    tidyr::drop_na(peakId) %>% 
    filter(rtime > rtmin - 15 & rtime < rtmax + 15)
  
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
  
  distances = lapply(data_MS2CombNorm@listData, calc_distances) %>% bind_rows(.id="peakId")
  
  distances = dplyr::left_join(distances, annotations, by="peakId") %>%
    dplyr::select(peakId, f1, f2, mzDifference, precursorMz, rtime, name, ion) %>%
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
  data.filtered = as.data.frame(rbindlist(data.filt.list, idcol = "MSMS.fragment"))
  
  
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
    dplyr::select(-mzDifference, -f1, -f2) %>%
    dplyr::distinct(peakId, .keep_all = T) %>% 
    dplyr::right_join(data.filtered.sum, by=c("sample","peakId"))
  
  return(list(data.filtered, spectra.filtered))
}

#' MS/MS Feature plotting 
#'
#' @description 
#' @param
#' @export
#' @seealso 
#' 
Dulce_MS2FeaturePlot = function(MSpectra, distances, peakId){
  
  ymax = 1.1
  
  df = distances[[1]] %>% dplyr::filter(peakId==!!peakId)
  df2 = data.frame(mz = MSpectra[[df$peakId[1]]]@mz,
                   intensity = MSpectra[[df$peakId[1]]]@intensity)
  
  ggplot() + 
    ggplot2::geom_segment(data=df, aes(yend=seq(ymax, 0.5, length.out=length(f1)),
                                    xend=f1,
                                    x=f1,
                                    y=0,
                                    color=MSMS.fragment), 
                       linetype="dashed") +
    ggplot2::geom_segment(data=df, aes(yend=seq(ymax, 0.5, length.out=length(f1)),
                                       xend=f2,
                                       x=f2,
                                       y=0,
                                       color=MSMS.fragment), 
                          linetype="dashed") +
    ggplot2::geom_segment(data=df2, aes(x=mz, xend=mz, y=0, yend=intensity)) +
    ggplot2::geom_segment(data=df, aes(x=f1, 
                              xend=f2, 
                              y=seq(ymax, 0.5, length.out=length(f1)), 
                              yend=seq(ymax, 0.5, length.out=length(f1)),
                              color=MSMS.fragment),
                 size=1.3) + 
    ggplot2::labs(x = expression(italic(m/z)), 
         y="Relative intesity",
         title = paste0(df$sample[1]," - ",
                        df$name[1],":",
                        df$ion[1],", ",
                        round(df$rtime[1]/60, 1)," min"),
         color="Loss", fill="Loss") + 
    ggplot2::geom_point(data=df, aes(x=precursorMz, y = ymax),
               shape = 25, 
               size = 3,
               fill = "black", 
               color = "black") +
    ggplot2::geom_point(aes(x = 0, y = ymax), 
               shape = 21, 
               fill = "black",
               colour = "black", 
               size = 3) +
    ggplot2::geom_segment(aes(x = 0, y = ymax, xend = 0, yend = ymax*0.9),
                 arrow = arrow(length = unit(0.1, "cm")),
                 colour = "black") +
    ggplot2::geom_point(aes(x=0, y=ymax*0.89), 
               shape=21, 
               fill="white",
               color="black", 
               size=3) +
    ggplot2::geom_text(data=df, aes(x = 20, y = ymax, label = round(precursorMz,4)),
              hjust = 0) +
    ggplot2::geom_text(aes(x = 20, y = ymax*0.95), label = "CID",
              hjust = 0) +
    ggplot2::geom_text(data = df2 %>% dplyr::filter(mz %in% c(df$f1, df$f2)), aes(x = mz, y = intensity,
                              label = as.character(round(mz,4))),
              angle = 90, vjust = 0.5, nudge_y=ymax/20,
              hjust = 0)+
    ggplot2::scale_color_manual(values=c("#F2B84B", "#E28139", "#734460", "#D1D1FE")) + 
    ggplot2::scale_fill_manual(values=c("#F2B84B", "#E28139", "#734460", "#D1D1FE")) +
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5))
}


