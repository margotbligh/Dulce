#' Pick and Group Peaks
#'
#' @description This function picks and groups peaks using the `findChromPeaks()`
#' and `groupChromPeaks()` functions from the `xcms` package.
#' @export
#' @param data  `XCMSnExp` object.
#' @param cwp   `CentWaveParam()` object. If not given, default
#' `CentWaveParam()` values are used.
#' @param pdp   `PeakDensityParam()` object. If not given, default
#' `PeakDensityParam()` values are used.
#' @param return_everything Logical value.
#'
#' @return If `return_everything` is set to FALSE (default), a `XCMSnExp` object is
#' returned. If TRUE, a list with (1) `XCMSnExp` object, (2) chromPeaks data.frame
#' and (3) featureDefinitions data.frame, is returned.
#'
#' @examples
#' # Examples have to be made with a toy data object.
#'
Dulce_pickGroupPeaks = function(data,
                                cwp=NULL,
                                pdp=NULL,
                                return_everything=F){

  # Check if data is an "OnDiskMSnExp" or "MSnExp" object.
  if (!class(data) %in% c("OnDiskMSnExp","MSnExp")){
    stop("Dulce error: 'data' is not from 'OnDiskMSnExp' or 'MSnExp' class. Check for ?readMSData.")
  }

  # Check for given parameters.
  # If no parameters are defined, default "CentWaveParam" parameters are used.
  # If parameters are not from "CentWaveParam" class, then function stops.
  if (is.null(cwp)){cwp = xcms::CentWaveParam()}
  else if (!is.null(cwp) & class(cwp)!="CentWaveParam"){
    stop("Dulce error: 'cwp' is not an 'CentWaveParam' object.
         Check for ?CentWaveParam or set it as NULL to use default arguments.")}

  # If no parameters are defined, default "PeakDensityParam" parameters are used.
  # If parameters are not from "PeakDensityParam" class, then function stops.
  if (is.null(pdp)){pdp = xcms::PeakDensityParam(sampleGroups=rep("Ungrouped", nrow(data)))}
  else if (class(pdp)!="PeakDensityParam"){
    stop("Dulce error: 'pdp' is not an 'PeakDensityParam' object. Check for ?PeakDensityParam or set it as NULL to use default arguments.")}

  # If no grouping is defined, default "Ungrouped" category is defined.
  if (length(pdp@sampleGroups)==0){
    pdp@sampleGroups = rep("Ungrouped", nrow(data))
    message("Dulce warning: No sample groups were defined. One group under the name of 'Ungrouped' will be created.")}

  # Pick and Group peaks
  data = xcms::findChromPeaks(data, param=cwp) %>% xcms::groupChromPeaks(param=pdp)

  # If return_everything set to TRUE, then also return a data.frame
  # for Peaks data and a data.frame for Groups data.
  if (return_everything){
    peaks_data = as.data.frame(data@msFeatureData[["chromPeaks"]])
    features_data = as.data.frame(data@msFeatureData[["featureDefinitions"]])

    return(list(data=data,
                peaks=peaks_data,
                features=features_data))
  }
  return(data)
}


#' Transform \code{XCMSnExp} object into a \code{xcmsSet} object.
#'
#' @description This functions transforms a `XCMSnExp` object into a `xcmsSet`
#' object. It will also classify and name the given samples if no name or
#' grouping class is specified, making the corresponding warnings.
#'
#' @param data `XCMSnExp` object
#' @param names Character vector with length equals to the amount of samples.
#' It defines the identifier of each sample. If `NULL` (default),
#' it generates names according to: `sprintf('sample_\%03d', 1:nrow(data))`
#' (sample_001, sample_002, etc).
#' @param classes Character vector with length equals to the amount of samples.
#' It defines the main grouping of the samples. If `NULL` (default), it defines
#' an "Unclassified" group for all samples.
#'
#' @return This function returns a `xcmsSet` object.
#'
#' @export
#'
#' @examples
#' # Examples have to be made with a toy data object.
#'
#' @seealso
#' xcms::xcmsSet
#' Dulce_AnnotateMS1
#'
Dulce_to_xcmsSet = function(data, names=NULL, classes=NULL){

  # Check if data is an "XCMSnExp" object
  if (class(data)!="XCMSnExp"){stop("Dulce error: 'data' object is not from 'XCMSnExp' class.")}

  # Transform data to "xcmsSet" object
  suppressMessages({data_xcmsSet = as(data, "xcmsSet")})

  # If not given, create sample names and sample class.
  if (is.null(names)){names = sprintf("sample_%03d", 1:nrow(data))}
  if (is.null(classes)){
    classes = "Unclassified"
    message("Dulce warning: No sample class given. One class under the name of 'Unclassified' will be created.")
  }

  # Set sample name and sample class
  xcms::sampnames(data_xcmsSet) = names
  xcms::sampclass(data_xcmsSet) = classes

  return(data_xcmsSet)
}



#' Find Isotopes and adducts.
#'
#' @description This function uses the CAMERA package to annotate the isotopes
#' and adducts of every feature.
#' @param data `xcmsSet` object.
#' @param perfwhm numeric value. It will be passed into `CAMERA::groupFWHM()` function
#' @param cor_eic_th numeric value. It will be passed into `CAMERA::groupCorr()` function
#' @param polarity character, either "positive" or "negative". It will be passed into
#' `CAMERA::groupCorr()` function
#'
#' @return
#' This function returns an `xsAnnotate` object with annotated isotopes and adducts.
#'
#' @export
#'
#' @examples
#' # Examples have to be made with a toy data object.
#'
#' @seealso
#' Dulce_AnnotateMS1
#'
Dulce_find = function(data,
                      perfwhm=0.5,
                      cor_eic_th=0.75,
                      polarity=NULL){

  # Check if data is an "xcmsSet" object
  if (class(data)!="xcmsSet"){stop("Dulce error: 'data' object is not from 'xcmsSet' class.")}

  if (is.null(polarity)){
    stop("Dulce error: Please specify polarity.")
  }

  # Apply CAMERA pipeline for finding isotopes and adducts.
  data = data %>%
    CAMERA::xsAnnotate() %>%
    CAMERA::groupFWHM(perfwhm=perfwhm) %>%
    CAMERA::findIsotopes() %>%
    CAMERA::groupCorr(cor_eic_th=cor_eic_th) %>%
    CAMERA::findAdducts(polarity=polarity)

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
#' @return This function returns a `data.frame` with collapsed information about the isotopic
#' signatures. A new column called "isotopes" is created, with a comma separated string
#' of the trimmed isotopes.
#'
#' @export
#'
#' @examples
#' # Examples have to be made with a toy data object.
#'
#' @seealso
#' Dulce_AnnotateMS1
#'
Dulce_trimIsotopes = function(data, rtmin=0, rtmax=Inf){

  # Check if "data" is an "xsAnnotate" object.
  if (class(data)!="xsAnnotate"){stop("Dulce error: 'data' object is not from 'xsAnnotate' class.")}

  # Filter peaks by custom retention time
  data = data %>%
    CAMERA::getPeaklist() %>%
    dplyr::filter(dplyr::between(rt, !!rtmin, !!rtmax))

  # Filter out rows without isotopes.
  # Collapse isotopes in one row and keep information about them.
  data_isotopes = data %>%
    dplyr::filter(isotopes!="") %>%
    dplyr::mutate(isotope_group=sub(x=isotopes, pat="\\[M.*", rep="")) %>%
    dplyr::group_by(isotope_group) %>%
    dplyr::mutate(isotopes = paste(isotopes, collapse=",")) %>%
    dplyr::distinct(isotope_group, .keep_all=T)

  # Join isotopes with no isotopes data.
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
#' @param pgp \code{predictGlycansParam()} object.
#' @param ppm M/z range for putative glycans in ppm units.
#' @param mzabs M/z range for putative glycans in absolute units. Changed to \code{NULL} if both ppm and mzabs are specified.
#'
#' @return This function returns a `data.frame` object. It is the result of
#' merging the putative glycans `data.frame` and the given features.
#'
#' @export
#'
#' @examples
#' # Examples have to be made with a toy data object.
#'
#' @seealso
#' Dulce_AnnotateMS1
#'
Dulce_featuresAnnotate = function(data, pgp=NULL,
                          ppm=NULL, mzabs=NULL){

  # Check if data has 'mzmin' and 'mzmax' columns
  if (!all(c("mzmin","mzmax") %in% colnames(data))){
    stop("Dulce error: check if 'mzmin' and 'mzmax' are columns in the 'data' object. They are needed to overlap the predictions.")
  }

  # Check if ppm or mzabs have valid values.
  if (is.null(ppm) & is.null(mzabs)){stop("Dulce error: 'ppm' and 'mzabs' have NULL values. Please specify one.")}
  if ((!is.null(ppm) & !is.null(mzabs))){
    message("Dulce warning: 'ppm' and 'mzabs' were specified (it is one or the other, not both). Using only 'ppm'")
    mzabs = NULL
  }

  # If no parameters are defined, default "predictGlycansParam" parameters are used.
  # If parameters are not from "predictGlycansParam" class, then function stops.
  if (is.null(pgp)){pgp = predictGlycansParam()}
  else if (!is.null(pgp) & class(pgp)!="predictGlycansParam"){
    stop("Dulce error: 'pgp' is not an 'predictGlycansParam' object. Check for ?predictGlycansParam or set it as NULL to use default arguments.")}

  # Use predictGlycans and turn into long format.
  predicted = glycanPredict::predictGlycans(param=pgp) %>%
    tidyr::pivot_longer(-c(name, dp, mass, formula), names_to = "ion", values_to = "mz") %>%
    tidyr::drop_na()

  # Create mz window for later overlap.
  if (!is.null(mzabs)){
    predicted = predicted %>%
      dplyr::mutate(mzmin=mz-mzabs,
                    mzmax=mz+mzabs)
  } else if (!is.null(ppm)){
    predicted = predicted %>%
      dplyr::mutate(mzmin=mz-ppm_to_mz(mz, ppm),
                    mzmax=mz+ppm_to_mz(mz, ppm))
  }

  # Overlap predicted glycans with features in data.
  data.table::setDT(predicted)
  data.table::setDT(data)
  data.table::setkey(predicted, mzmin, mzmax)

  predicted = data.table::foverlaps(data,predicted) %>%
    tidyr::drop_na(name) %>%
    as.data.frame() %>%
    dplyr::rename(mz.predicted=mz,
                  mzmin.predicted=mzmin,
                  mzmax.predicted=mzmax,
                  mz=i.mz,
                  mzmin=i.mzmin,
                  mzmax=i.mzmax)

  return(predicted)
}



#' Annotate glycans in MS1 data
#'
#'
#' @description This function constitutes the complete MS Data Annotation Pipeline. It is robust and smart
#' enough to *understand* what kind of data is being inputed and continue the pipeline from the appropriate
#' step. It is constructed in such a way that when an error occurs, the last correctly executed product is
#' saved and returned.
#'
#' @param data  Object of any of these classes: MSnExp, OnDiskMSnExp, XCMSnExp, xcmsSet,
#' xsAnnotate, data.frame.
#' @param cwp CentWaveParam() object. If not given, default CentWaveParam() values are used.
#' @param pdp PeakDensityParam() object. If not given, default PeakDensityParam() values are used.
#' @param names Character vector with length equals to the amount of samples. It defines
#' the identifier of each sample. If `NULL` (default), it generates names according
#' to: `sprintf('sample_\%03d', 1:nrow(data))` (sample_001, sample_002, etc).
#' @param classes Character vector with length equals to the amount of samples. It defines the main grouping of the samples. If `NULL` (default), it defines an "Unclassified" group for all samples.
#' @param perfwhm numeric value. It will be passed into `CAMERA::groupFWHM()` function.
#' @param cor_eic_th numeric value. It will be passed into `CAMERA::groupCorr()` function.
#' @param polarity character value, either "positive" or "negative". It will be passed into `CAMERA::findAdducts()` function.
#' @param rtmin Numeric value specifying the lower retention time boundary.
#' @param rtmax Numeric value specifying the higher retention time boundary.
#' @param pgp `predictGlycansParam()` object. If not given, default predictGlycansParam() values are used.
#' @param ppm  M/z range for putative glycans in ppm units.
#' @param mzabs  M/z range for putative glycans in absolute units. Changed to `NULL` if both, ppm and mzabs, are specified.
#' @param output Character object that determines what the function will return.
#' If \code{"main"}, a list with an \code{'XCMSnExp'} object (if possible) and the last achieved product will be returned.
#' If \code{"all"}, a list with all the partial products along the pipeline is returned.
#' If \code{"last"}, only the last achieved product will be returned.
#'
#' @return If no errors occurs during the execution of the function, then the output depends on th
#'  `output` argument. If `all`, the function will return every partial product achieved. If `main`,
#'  XCMSnExp object (with peaks picked and grouped) and an annotated putative glycans `data.frame`  are
#'  returned. If `last`, only the annotated putative glycans `data.frame` is returned. If an error occur
#'  during the execution of the function, it saves the last partial product achieved, and tries to return
#'  the specified output option.
#'
#' @export
#'
#' @details
#'
#' @examples
#' # Examples have to be made with a toy data object.
#'
#' @seealso
#'
Dulce_AnnotateMS1 = function(data, cwp=NULL, pdp=NULL,
                       names=NULL, classes="Unclassified",
                       perfwhm=0.5, cor_eic_th=0.75,
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

  # Try 'Dulce_pickGroupPeaks' function. If it does not work, then return the
  # imputed value and continue running the function.
  if (class(data) %in% c("OnDiskMSnExp","MSnExp")){
    message("Dulce note: executing Dulce_pickGroupPeaks function.")
    data = tryCatch(Dulce_pickGroupPeaks(data, cwp=cwp, pdp=pdp, return_everything=T),
                    error=function(e) {message("Dulce_pickGroupPeaks() did not work due to the following error:\n", e); data})


    if (output == "all"){
      data_all$peaks = data$peaks
      data_all$features = data$features
      data_all$XCMSnExp = data$data
      data = data$data
    } else if (output == "main"){
      data_main$XCMSnExp = data$data
      data = data$data
    } else if (output == "last"){
      data = data$data
    }

  } else {
    message("Dulce warning: 'data' not from 'OnDiskMSnExp' or 'MSnExp' class.")
    message("Dulce warning: trying next step in pipeline, Dulce_to_xcmsSet().")
  }

  # Try 'Dulce_to_xcmsSet' function. If it does not work, then return the imputed value, rise the error message and continue running the function.
  if (class(data)=="XCMSnExp"){
    message("Dulce note: executing Dulce_to_xcmsSet function.")
    data = tryCatch(Dulce_to_xcmsSet(data, names, classes),
                    error=function(e) {message("Dulce_to_xcms() did not work due to the following error:\n", e); data})


    if (output == "all"){data_all$xcmsSet = data}
    } else {
      message("Dulce warning: 'data' not from 'XCMSnExp' class.")
      message("Dulce warning: trying next step in pipeline, Dulce_find()")
    }

  # Try 'Dulce_find' function. If it does not work, then return the imputed value, rise the error message and continue running the function.
  if (class(data)=="xcmsSet"){
    message("Dulce note: executing Dulce_find function.")
    data = tryCatch(Dulce_find(data,
                              perfwhm=perfwhm, cor_eic_th=cor_eic_th,
                              polarity=polarity),
                    error=function(e) {message("Dulce_find() did not work due to the following error:\n", e); data})


    if (output == "all"){data_all$xsAnnotate = data}
    } else {
      message("Dulce warning: 'data' not from 'xcmsSet' class.")
      message("Dulce warning: trying next step in pipeline, Dulce_trimIsotopes()")
    }

  # Try 'Dulce_find' function. If it does not work, then return the imputed value, rise the error message and continue running the function.
  if (class(data)=="xsAnnotate"){
    message("Dulce note: executing Dulce_trimIsotopes function.")
    data = tryCatch(Dulce_trimIsotopes(data, rtmin=rtmin, rtmax=rtmax),
                    error=function(e) {message("Dulce_trimIsotopes() did not work due to the following error:\n", e); data})


    if (output == "all"){data_all$trimmedIsotopes = data}
    } else {
      message("Dulce warning: 'data' not from 'xsAnnotate' class.")
      message("Dulce warning: trying next step in pipeline, Dulce_annotate()")
    }

  # Try 'Dulce_featuresAnnotate' function. If it does not work, then return the imputed value, rise the error message and continue running the function.
  if (any(class(data)=="data.frame")){
    message("Dulce note: executing Dulce_featuresAnnotate function.")
    data = tryCatch(Dulce_featuresAnnotate(data, pgp=pgp, ppm=ppm, mzabs=mzabs),
                    error=function(e) {message("Dulce_featuresAnnotate() did not work due to the following error:\n", e); data})


  # Return the output depending on the "output" argument.
    if (output == "all"){data_all$Annotated = data
    } else if (output == "main"){data_main$Annotated = data}

  } else {
      message("Dulce warning: 'data' not even from 'data.frame' class. Check data object.")
  }

  if (output == "all"){
    message("Dulce note: returning a list with the product of all the steps that were done successfully.")
    return(data_all)
  } else if (output == "main"){
    message("Dulce note: returning a list with an 'XCMSnExp' object (if possible) and the last achieved product.")
    return(data_main)
  } else if (output == "last"){
    message("Dulce note: returning the last achieved object.")
    return(data)
  }
}






