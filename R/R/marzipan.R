#' @title Marzipan Dataset
#'
#' @description VIS-NIR Spectra for a compositional analysis of Marzipan samples. 32 samples were measured, from 9 different recipes.
#'
#' @format A list with 4 objects:
#' \describe{
#'   \item{ax}{Wavelength vector (400 to 2498 nm in 2nm intervals)}
#'   \item{X}{NIR Spectra of 32 marzipan samples, measured using a NIRSystems 6500 spectrometer}
#'   \item{Y}{Compositional data of the 32 samples; 1st Column contains the moisture content, second column gives the sugar content of each sample}
#'   \item{recipe}{recipe labels for each sample}
#' }
#' ' @examples
#' data(marzipan)
#' @references J Christensen, L Nørgaard, H Heimdal, JG Pedersen, SB Engelsen. Rapid Spectroscopic Analysis of Marzipan – Comparative Instrumentation. Journal of Near Infrared Spectroscopy, vol 12, No. 1 (2004).
#' @source <https://ucphchemometrics.com/datasets/>
#'
