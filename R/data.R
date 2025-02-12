#' Northern shoveler data
#'
#' Repeated count data of Northern shoveler with covariates, 
#' formatted for use with the \link[unmarked]{unmarked} package.
#'
#' @docType data
#' 
#' @format A list with three elements
#' \describe{
#'   \item{y}{A matrix with Northern shoveler counts}
#'   \item{site}{A data frame with site specific covariates}
#'   \item{obs}{A list containing observation specific covariates}
#' }
#' 
#' @references Knape et al. (2018) Methods in Ecology and Evolution, 9:2102-2114. \doi{10.1111/2041-210X.13062}
#' 
#' @examples
#' library(unmarked)
#' umf = unmarkedFramePCount(y = shoveler$y, obsCovs = shoveler$obs, siteCovs = shoveler$site)
"shoveler"