#' Experimental data of Venus dimers dynamics in HEK-293 cells.
#'
#' This data set consists on a raster line scan performed over HEK-293 cells expressing dimers of the fluorescent protein Venus, also known as SEYFP-F46L.
#' The scan line is 64 pixels long, and the scanning direction is from the cytoplasm to the nucleus, across the nuclear envelope.
#' A pixel size of 50 nm was used, as well as a pixel dwell time of 12.5 us and a line scan time of 1.925 ms.
#' Fluorescence excitation was provided by a 488 nm laser at 0.1% power.
#' Fluorescence intensity data was collected using the photon-counting mode in an Olympus FV1000 Upright BX61WI confocal microscope.
#'
#' @docType data
#' @keywords datasets
#' @name v2DataSet
#' @usage data(v2DataSet)
#' @format A data frame with 64 rows and 20000 columns
#' @details In order to use the data in pcfData, the data.matrix() must be used to transform the data set into a matrix, and can be used in the examples.
NULL