#' Cosine similarity
#'
#' @description Calculate cosine similarity between two vectors of the same length, or row-wise
#' cosine similarity of two matrices of the same dimensions
#'
#' @param x A numeric vector or matrix
#' @param y A numeric vector or matrix
#'
#' @return A numeric vector or matrix of cosine similarities
#' @rdname cosSim
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(runif(200), nrow=20)
#' y <- matrix(runif(200), nrow=20)
#'
#' ## Works with matrices
#' cosSim(x, y)
#'
#' ## ...or vectors
#' cosSim(x[1,], y[1,])
#'
cosSim <- function(x, ...) {
   UseMethod("cosSim", x)
}

#' @rdname cosSim
#' @method cosSim default
#' @export
cosSim.default <- function(x, y){
   if(length(x)!=length(y)){ stop('`x` and `y` must be the same length') }

   ## Convert to numeric to prevent integer overflow
   x <- as.numeric(x)
   y <- as.numeric(y)

   sum(x*y)/sqrt(sum(x^2)*sum(y^2))
}

#' @rdname cosSim
#' @method cosSim matrix
#' @export
cosSim.matrix <- function(x, y) {
   if(!identical(dim(x), dim(y))){ stop('`x` and `y` must have the same dimensions') }

   ## Convert to numeric to prevent integer overflow
   x <- matrix(as.numeric(x), nrow=nrow(x))
   y <- matrix(as.numeric(y), nrow=nrow(y))

   rowSums(x*y) / sqrt( rowSums(x^2) * rowSums(y^2) )
}

#' @rdname cosSim
#' @method cosSim data.frame
#' @export
cosSim.data.frame <- function(x,y){
   x <- as.matrix(x)
   y <- as.matrix(y)
   cosSim.matrix(x, y)
}

