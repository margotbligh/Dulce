## Toolbox of the package


#' Fetch: find and group picks
#'
#' @description This functions is a robust pipeline that merges the 
#' \code{findChromPeaks()} function and \code{groupChromPeaks()} function. 
#' 
#' @param data  \code{XCMSnExp} object
#' @param cwp   \code{CentWaveParam()} object. If not given, default \code{CentWaveParam()} values are used. 
#' @param pdp   \code{PeakDensityParam()} object. If not given, default \code{PeakDensityParam()} values are used.
#' @param return_everything Logical value. If FALSE, (default) a \code{XCMSnExp} object is returned. 
#' If TRUE, a list, with (1) \code{XCMSnExp} object, (2) chromPeaks data.frame and (3) featureDefinitions data.frame, is returned.
#' 
#' @export
#'  
#' @details 
#' Dobie dobie doo diroo diroo
#' 
#' 
#' @examples
#' # Examples have to be made with a toy data object.
#' 
#' @seealso 
#' xcms::CentWaveParam()
#' 
Dulce_fetch = function(data, 
                       cwp=NULL, 
                       pdp=NULL, 
                       return_everything=F){
  
  if (!class(data) %in% c("OnDiskMSnExp","MSnExp")){
    stop("Dulce error: 'data' is not from 'OnDiskMSnExp' or 'MSnExp' class. Check for ?readMSData.")
  }
  
  # If no param objects defined, default is used.
  if (is.null(cwp)){cwp = xcms::CentWaveParam()}
  else if (!is.null(cwp) & class(cwp)!="CentWaveParam"){
    stop("Dulce error: 'cwp' is not an 'CentWaveParam' object. 
         Check for ?CentWaveParam or set it as NULL to use default arguments.")}
  
  if (is.null(pdp)){pdp = xcms::PeakDensityParam(sampleGroups=rep("Ungrouped", nrow(data)))}
  else if (class(pdp)!="PeakDensityParam"){
    stop("Dulce error: 'pdp' is not an 'PeakDensityParam' object. Check for ?PeakDensityParam or set it as NULL to use default arguments.")}
  
  # Defining a default "Ungrouped" category if no grouping is defined.
  if (length(pdp@sampleGroups)==0){
    pdp@sampleGroups = rep("Ungrouped", nrow(data))
    message("Dulce warning: No sample groups were defined. One group under the name of 'Ungrouped' will be created.")}
  
  # Fetch 
  suppressMessages({
  processed_data = xcms::findChromPeaks(data, param=cwp) %>% xcms::groupChromPeaks(param=pdp)
  })
  message("Peaks picked and grouped!")
  
  if (return_everything){
    peaks_data = as.data.frame(processed_data@msFeatureData[["chromPeaks"]])
    features_data = as.data.frame(processed_data@msFeatureData[["featureDefinitions"]])
    
    return(list(data=processed_data, 
                peaks=peaks_data, 
                features=features_data))
  }
  
  return(processed_data)
}



#' to_xcmsSet: transforms a \code{XCMSnExp} object into a \code{xcmsSet} object. 
#'
#' @description Da da da doobie doo da da
#' 
#' @param data  \code{XCMSnExp} object.
#' @param names Character vector with length equals to the amount of samples. It defines the identifier of each sample. 
#' If NULL (default), it generates names according to: \code{sprintf('sample_\%03d', 1:nrow(data))}.
#' @param classes Character vector with length equals to the amount of samples. It defines the main grouping of the samples.
#' If NULL (default), it defines an "Unclassified" group for all samples.
#' 
#' @return \code{xcmsSet} object.
#' @export
#'  
#' @details 
#' Dobie dobie doo diroo diroo
#' 
#' @examples
#' # Examples have to be made with a toy data object.
#' 
#' @seealso 
#' xcms::xcmsSet
#' 
Dulce_to_xcmsSet = function(data, names=NULL, classes=NULL){
  
  if (class(data)!="XCMSnExp"){stop("Dulce error: 'data' object is not from 'XCMSnExp' class.")}
  
  suppressMessages({data_xcmsSet = as(data, "xcmsSet")})
  
  if (is.null(names)){names = sprintf("sample_%03d", 1:nrow(data))}
  if (is.null(classes)){
    classes = "Unclassified"
    message("Dulce warning: No sample class given. One class under the name of 'Unclassified' will be created.")
  }
  
  sampnames(data_xcmsSet) = names
  sampclass(data_xcmsSet) = classes
  
  return(data_xcmsSet)
}



#' This function is the one I like the less.
#'
#' @description I will still use it in the main function, but major changes 
#' have to be made. 
#' 
#' @param bla bla bla bla
#' 
#' @return 
#' 
#' @export
#'  
#' @details Dobie dobie doo diroo diroo
#' 
#' @examples
#' # Examples have to be made with a toy data object.
#' 
#' @seealso
#' 
Dulce_find = function(data, isotopes=T, adducts=T, 
                      perfwhm=0.5, mzabs=0.01, cor_eic_th=0.75,
                      polarity=NULL){
  
  if (is.null(polarity)){
    stop("Dulce error: NULL polarity? Check it twice.")
  }
  
  data = CAMERA::xsAnnotate(data) %>% CAMERA::groupFWHM(perfwhm = perfwhm)
  
  if (isotopes){
    data = data %>% CAMERA::findIsotopes(mzabs=mzabs) %>% 
      CAMERA::groupCorr(cor_eic_th=cor_eic_th) 
  }
  
  if (adducts){
    data = data %>% CAMERA::findAdducts(polarity=polarity)
  }
  
  return(data)
}



#' Trim Isotopes
#'
#' @description This functions receives a \code{xsAnnotate} object and trims its isotopes,
#' leaving only the rows with no isotopes, and collapsing those with isotopes into one row 
#' (annotating those which were found in the "isotopes" column). 
#' 
#' @param data \code{xsAnnotate} object.
#' @param rtmin Minimum retention time value (in seconds) from which scans should be considered.
#' @param rtmax Maximum retention time value (in seconds) to which scans should be considered.
#' 
#' @return data.frame
#' 
#' @export
#'  
#' @details Dobie dobie doo diroo diroo
#' 
#' @examples
#' # Examples have to be made with a toy data object.
#' 
#' @seealso
#' 
Dulce_trimIsotopes = function(data, rtmin=0, rtmax=Inf){
  
  if (class(data)!="xsAnnotate"){stop("Dulce error: 'data' object is not from 'xsAnnotate' class.")}
  
  data = CAMERA::getPeaklist(data) %>% dplyr::filter(data.table::between(rt, !!rtmin, !!rtmax))
  
  data_isotopes = data %>% 
    dplyr::filter(isotopes!="") %>% 
    dplyr::mutate(isotope_group=sub(x=isotopes, pat="\\[M.*", rep="")) %>% 
    dplyr::group_by(isotope_group) %>% 
    dplyr::mutate(isotopes = paste(isotopes, collapse=",")) %>% 
    dplyr::distinct(isotope_group, .keep_all=T)
  
  data = data %>% 
    dplyr::filter(isotopes=="") %>% 
    dplyr::bind_rows(data_isotopes) %>% 
    dplyr::select(-isotope_group)
  
  return(data)
}



#' Annotate putative glycans
#' 
#' @include utils.r
#' 
#' @description This functions receives a \code{data.frame} object and merges it with the 
#' output of \code{glycanPredict::predictGlycans()} function. \code{data} must have a
#' 'mzmin' and a 'mzmax' columns. 
#' 
#' @param data \code{data.frame} object with 'mzmin' and 'mzmax' columns. 
#' @param pgp \code{predictGlycansParams()} object.  
#' @param ppm M/z range for putative glycans in ppm units.
#' @param mzabs M/z range for putative glycans in absolute units. Changed to \code{NULL} if both ppm and mzabs are specified.
#' 
#' @return data.frame
#' 
#' @export
#'  
#' @details Dobie dobie doo diroo diroo
#' 
#' @examples
#' # Examples have to be made with a toy data object.
#' 
#' @seealso
#' 
Dulce_annotate = function(data, pgp=NULL,
                          ppm=NULL, mzabs=NULL){
  
  if (!all(c("mzmin","mzmax") %in% colnames(data))){
    stop("Dulce error: check if 'mzmin' and 'mzmax' are columns in the 'data' object. They are needed to overlap the predictions.")
  }
  
  if (is.null(ppm) & is.null(mzabs)){stop("Dulce error: 'ppm' and 'mzabs' have NULL values. Please specify one.")}
  if ((!is.null(ppm) & !is.null(mzabs))){
    message("Dulce warning: 'ppm' and 'mzabs' were specified (it is one or the other, not both). Using only 'ppm'")
    mzabs = NULL
  }
  
  if (is.null(pgp)){pgp = predictGlycansParam()}
  else if (!is.null(pgp) & class(pgp)!="predictGlycansParam"){
    stop("Dulce error: 'pgp' is not an 'predictGlycansParam' object. Check for ?predictGlycansParam or set it as NULL to use default arguments.")}
  
  
  predicted = glycanPredict::predictGlycans(param=pgp) %>% 
    tidyr::pivot_longer(-c(name, dp, mass, formula), names_to = "ion", values_to = "mz") %>% 
    tidyr::drop_na()
  
  if (!is.null(mzabs)){
    predicted = predicted %>% 
      dplyr::mutate(mzmin=mz-mzabs, 
                    mzmax=mz+mzabs)
  } else if (!is.null(ppm)){
    predicted = predicted %>% 
      dplyr::mutate(mzmin=mz-ppm_to_mz(mz, ppm),
                    mzmax=mz+ppm_to_mz(mz, ppm))
  }
  
  data.table::setDT(predicted)
  data.table::setDT(data)
  data.table::setkey(predicted, mzmin, mzmax)
  
  predicted = data.table::foverlaps(data,predicted) %>% 
    tidyr::drop_na(name) %>% 
    dplyr::select(-i.mzmin, -i.mzmax, -i.mz) %>% 
    as.data.frame()
  
  return(predicted)
}



#' Annotate glycans in MS1 data
#' 
#' 
#' @description Run the complete pipeline. Lots of things are gonna be assumed, many
#' defaults values would be taken... but hey! You will get something at the end.
#' Does it have any value? maybe, only time will tell...
#' ...
#' And you, of course you have to tell if it has value or not (?). You are the one who should know!
#' Read some literature, don't be lazy.
#' 
#' @param data any class of object... Dulce message: Try me out baby.
#' @param output Character that determines what the function will return. 
#' If \code{"main"}, a list with an \code{'XCMSnExp'} object (if possible) and the last achieved product will be returned.
#' If \code{"all"}, a list with all the partial products along the pipeline is returned. 
#' If \code{"last"}, only the last achieved product will be returned. 
#' @return The last produced object before an error (or at least that is the idea)
#' 
#' @export
#'  
#' @details Dobie dobie doo diroo diroo
#' 
#' @examples
#' # Examples have to be made with a toy data object.
#' 
#' @seealso
#' 
Dulce_AnnotateMS1 = function(data, cwp=NULL, pdp=NULL,
                       names=NULL, classes="Unclassified",
                       isotopes=T, adducts=T,
                       perfwhm=0.5, mzabs.find=0.01, cor_eic_th=0.75,
                       polarity=NULL,
                       rtmin=0, rtmax=Inf,
                       pgp=NULL, ppm=NULL, mzabs=NULL,
                       output = "all"){
  
  data_all = list()
  data_main = list()
  
  if (!output %in% c("all", "main", "last")){
    stop("Dulce error: 'output' argument not recognized. The options are 'all', 'main' or 'last'. See ?Dulce_AnnotateMS1 for more information.")
  }
  
  BiocParallel::register(SerialParam())
  
  if (class(data) %in% c("OnDiskMSnExp","MSnExp")){
    message("Dulce note: executing Dulce_fetch function. Bip Bop... Bip Bop...")
    if (output == "all"){
      data = Dulce_fetch(data, cwp=cwp, pdp=pdp, return_everything=T)
      data_all$peaks = data$peaks
      data_all$features = data$features
      data_all$XCMSnExp = data$data
      data = data$data
    } else if (output == "main"){
      data = Dulce_fetch(data, cwp=cwp, pdp=pdp, return_everything=F)
      data_main$XCMSnExp = data
    } else if (output == "last"){
      data = Dulce_fetch(data, cwp=cwp, pdp=pdp, return_everything=F)
    }
  } else {
    message("Dulce warning: 'data' not from 'OnDiskMSnExp' or 'MSnExp' class.")
    message("Dulce warning: trying next step in pipeline, Dulce_to_xcmsSet().")
  }
  
  if (class(data)=="XCMSnExp"){
    message("Dulce note: executing Dulce_to_xcmsSet function. Diroo... Diroo Diroo... ")
    data = Dulce_to_xcmsSet(data, names, classes)
    if (output == "all"){data_all$xcmsSet = data}
    } else {
      message("Dulce warning: 'data' not from 'XCMSnExp' class.")
      message("Dulce warning: trying next step in pipeline, Dulce_find()")
    }
  
  if (class(data)=="xcmsSet"){
    message("Dulce note: executing Dulce_find function. Daroo bip... daroo bop... ")
    message("Dulce rants: I dont like this functions :c It needs to be improved.")
    data = Dulce_find(data, isotopes=isotopes, adducts=adducts, 
                     perfwhm=perfwhm, mzabs=mzabs.find, cor_eic_th=cor_eic_th,
                     polarity=polarity)
    if (output == "all"){data_all$xsAnnotate = data}
    } else {
      message("Dulce warning: 'data' not from 'xcmsSet' class.")
      message("Dulce warning: trying next step in pipeline, Dulce_trimIsotopes()")
    }
   
  if (class(data)=="xsAnnotate"){
    message("Dulce note: executing Dulce_trimIsotopes function. Beep... clink...")
    data = Dulce_trimIsotopes(data, rtmin=rtmin, rtmax=rtmax)
    if (output == "all"){data_all$trimmedIsotopes = data}
    } else {
      message("Dulce warning: 'data' not from 'xsAnnotate' class.")
      message("Dulce warning: trying next step in pipeline, Dulce_annotate()")
    }
    
  if (class(data)=="data.frame"){
    message("Dulce note: executing Dulce_annotate function. Wooosh! Birup... Birup... pa!")
    data = Dulce_annotate(data, pgp=pgp, ppm=ppm, mzabs=mzabs)
    if (output == "all"){data_all$Annotated = data}
    else if (output == "main"){data_main$Annotated = data} 
    
  } else {
      message("Dulce warning: 'data' not even from 'data.frame' class. What did you give me?")
  }
  
  if (output == "all"){
    message("Dulce note: returning a list with the product of all the steps I was able to do successfully.")
    return(data_all)
  } else if (output == "main"){
    message("Dulce note: returning a list with an 'XCMSnExp' object (if possible) and the last achieved product.")
    return(data_main)
  } else if (output == "last"){
    message("Dulce note: returning the last object I was able to produce.")
    return(data)
  }
}






