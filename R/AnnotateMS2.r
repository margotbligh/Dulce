#' MS/MS Spectra Extraction
#'
#' @description This function looks into an `MSnBase` object and extracts all MS/MS spectra related to MS1 features defined in the "annotations" object. It uses the `featureSpectra()` function.
#' @param data `MSnBase` original data object. It must have msLevel=2 spectra.
#' @param features `XCMSnExp` object with msLevel=1 and featured defined.
#' @param annotations features annotated `data.frame` with mz values column.
#'
#' @return This function returns a MSpectra object with all msLevel=2 spectra associated with the features annotated `data.frame`.
#'
#' @export
#' @seealso
#'
Dulce_MS2SpectraExtract = function(data, features, annotations){

  # Check if "data" has msLevel = 2 spectra.
  if(any(table(data@featureData@data$msLevel, data@featureData@data$fileIdx)[2,]==0)){
    stop("Dulce error: Some of the files do not have MS level 2 data.")
  }

  # Change "data" into "XCMSnExp"
  data = as(data, "XCMSnExp")

  # Overwrite "data" peaks and features for the ones in "features"
  chromPeaks(data) = chromPeaks(features)
  featureDefinition = featureDefinitions(features)
  featureDefinition_filtered = featureDefinition[paste0(featureDefinition$mzmed,"_", featureDefinition$rtmed) %in%
                                                   paste0(annotations$mz,"_",annotations$rt),]

  featureDefinitions(data) = featureDefinition_filtered

  # Get msLevel=2 spectra that is related to defined features.
  MS2Spectra = featureSpectra(data, msLevel = 2, expandMz = 0.005)

  return(MS2Spectra)
}


#' MS/MS Spectra Combination and Normalization
#' @include classes.r
#' @description This function cleans, combines and normalizes MS/MS spectra.
#' @param data `MSpectra` object with msLevel=2 spectra.
#' @param csp  `combineSpectraParam()` object.
#' @param `norm.method` normalization method for `MSnbase::normalise()` function.
#'
#' @return This function returns normalized and combined MS/MS spectra as `MSpectra` object.
#'
#' @export
#' @seealso
#'
Dulce_MS2CombNorm = function(data, csp=NULL, norm.method="max"){

  # Check if "data" is an "MSpectra" object.
  if (class(data)!="MSpectra"){stop("Dulce error: 'data' object not from 'MSpectra' class.")}
  # If no parameters are defined, default "combineSpectraParam" parameters are used.
  if (is.null(csp)){csp = combineSpectraParam()}

  # Clean, combine and normalize spectra from same precursor
  data@listData = lapply(data, clean, all=T)
  data = combineSpectra(data, fcol=csp@fcol, mzd=csp@mzd, intensityFun=csp@intensityFun)
  data@listData = lapply(data, MSnbase::normalise, method=norm.method)

  return(data)
}

#' MS/MS Precursor Matching
#'
#' @description This function generates an annotated `data.frame` of the MS/MS
#' spectra information with the precursor annotations. It uses the
#' function `data.table::foverlaps()` on the mzmin and mzmax values of both
#' groups of spectra.
#'
#' @param data `MSpectra` object.
#' @param annotations `data.frame` object with annotated MS features.
#'
#' @return This function returns a `data.frame` with annotated MS/MS spectra.
#'
#' @export
#' @seealso
#'
Dulce_MS2PrecursorMatch = function(data, annotations){

  # Gather spectra's precursor mz value and retention time.
  precursorsMz = precursorMz(data)
  precursorsRt = rtime(data)

  # Create data frame with precursor information
  precursorsData = data.frame(peakId = names(precursorsMz),
                              precursorMz = precursorsMz,
                              rtime = precursorsRt,
                              mzmin = precursorsMz - 0.0025,
                              mzmax = precursorsMz + 0.0025)

  setDT(precursorsData)
  setDT(annotations)
  setkey(precursorsData, mzmin, mzmax)

  # Overlap annotations and precursors
  data = foverlaps(annotations, precursorsData) %>%
    tidyr::drop_na(peakId) %>%
    filter(rtime > rtmin - 15 & rtime < rtmax + 15)

  return(data)
}

#' MS/MS Distances Calculation
#'
#' @description This function calculate mz differences between MS/MS fragments.
#' @param data `MSnBase` object
#' @param annotations `data.frame` with annotated MS/MS spectra.
#' @param MSpectra `MSpectra` object with MS/MS spectra.
#'
#' @return This function returns a `data.frame` object with mzDifferences,
#' their respective precursor, and the fragments used to calculate that distance.
#'
#' @export
#' @seealso
#'
Dulce_calculateMZDistances = function(data, annotations, MSpectra){

  # Calculate differences
  distances = lapply(MSpectra@listData, calc_distances) %>% bind_rows(.id="peakId")

  # Join differences with annotations
  distances = dplyr::left_join(distances, annotations, by="peakId") %>%
    dplyr::select(peakId, f1, f2, mzDifference, precursorMz, rtime, name, ion) %>%
    dplyr::distinct() %>%
    dplyr::mutate(name_ion = paste0(name, ":", ion)) %>%
    dplyr::group_by(across(-c(name, ion, name_ion))) %>%
    dplyr::summarise(name_ion=sapply(list(name_ion), paste, collapse=",")) %>%
    dplyr::ungroup()

  distances$sample = distances$peakId %>%
    gsub(".*\\.F|\\.S\\d+$", "", .) %>% as.numeric()

  # Join differences with samples
  distances$sample = pData(data)$name[distances$sample]

  # Return data frame with annotated differnces between mz values.
  return(distances)
}

#' MS/MS Distances Filtration and Highlighting
#'
#' @description This function filters the MS/MS spectra by the amount of
#' different meaningful differences.
#' @param data `data.frame` object with calculated distances.
#' @param sigdiff numeric vector with significant m/z differences values.
#' @param n numeric value. Indicates the amount of different significant
#' differences that a spectra should have in order to be included (not filtered out).
#' If `NULL`, then `n=length(sigdiff)`
#' @param mzdiff numeric value. Indicates the allowed range of mz values by which
#' differences should be filtered.
#'
#' @return This function returns a `data.frame` with the filtered spectra based
#' on the significant differences that were present.
#'
#' @export
#' @seealso
#'
Dulce_filterMZDistances = function(data, sigdiff=NULL, n=NULL, mzdiff=NULL){

  # Define significant differences
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

  # Identify significant differences
  for(i in 1:length(sigdiff.used)){
    data.filt.list[[i]] = data %>%
      filter(between(mzDifference,
                     sigdiff.used[i]-mzdiff, sigdiff.used[i]+mzdiff))
  }

  names(data.filt.list) = names(sigdiff.used)
  data.filtered = as.data.frame(rbindlist(data.filt.list, idcol = "MSMS.fragment"))


  # Count how many times each significant diff found in each spectrum
  data.filtered.sum <- data.filtered %>%
    dplyr::group_by(sample, peakId, MSMS.fragment) %>%
    dplyr::summarise(n = n())

  data.filtered.sum.w <- data.filtered.sum %>%
    pivot_wider(names_from = MSMS.fragment, values_from = n)

  # Count how many different significant differences there are per spectrum
  data.filtered.sum = data.filtered.sum %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, peakId) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n >= !!n) %>%
    dplyr::left_join(data.filtered.sum.w, by=c("sample","peakId"))

  # Filter the original data to only contain these spectra
  spectra.filtered = data %>%
    dplyr::select(-mzDifference, -f1, -f2) %>%
    dplyr::distinct(peakId, .keep_all = T) %>%
    dplyr::right_join(data.filtered.sum, by=c("sample","peakId"))

  return(list(data.filtered, spectra.filtered))
}

#' MS/MS Feature plotting
#'
#' @description This function plots each MS/MS spectra, highlighting the meaningful differences.
#' @param MSpectra `MSpectra` object with MS/MS spectra.
#' @param distances `data.frame` with significant differences annotated with spectrum identifier.
#' @param peakId Peak identifier of the peak that is going to be plotted.
#'
#' @return This function returns a `ggplot2` object.
#' @export
#' @seealso
#'
Dulce_MS2FeaturePlot = function(MSpectra, distances, peakId){

  ymax = 1.1

  # Filter distances by peakId
  df = distances %>% dplyr::filter(peakId==!!peakId)
  df2 = data.frame(mz = MSpectra[[df$peakId[1]]]@mz,
                   intensity = MSpectra[[df$peakId[1]]]@intensity)

  # Plot
  ggplot2::ggplot() +
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
                        df$name_ion[1],", ",
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


