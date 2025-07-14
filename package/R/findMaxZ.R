#' Find Maximum Z Value
#'
#' @param s Starting index (integer)
#' @param t Ending index (integer)
#' @param mincpgs Minimum number of CpGs (integer)
#' @param XS Two dimensional array/matrix
#'
#' @return A vector of two integers representing the indices with maximum Z value
#' @export
#' @useDynLib AMRs.finder, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @examples
#' \dontrun{
#' # Example usage
#' XS <- matrix(rnorm(100), ncol = 2)
#' result <- findMaxZC(1, 50, 3, XS)
#' }
findMaxZC <- function(s, t, mincpgs, XS) {
    # Input validation
    if (!is.numeric(s) || !is.numeric(t) || !is.numeric(mincpgs)) {
        stop("s, t, and mincpgs must be numeric")
    }
    if (!is.matrix(XS) && !is.array(XS)) {
        stop("XS must be a matrix or array")
    }
    if (ncol(XS) < 2) {
        stop("XS must have at least 2 columns")
    }
    
    # Ensure inputs are integers
    s <- as.integer(s)
    t <- as.integer(t)
    mincpgs <- as.integer(mincpgs)
    
    # Call the C++ function
    .Call(`_AMRSfinder_findMaxZ`, s, t, mincpgs, XS)
}