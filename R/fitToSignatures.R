#' Linear least-squares (non-negative) fitting
#'
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are the mutation
#' contexts and the columns are  the mutational signatures.
#' @param mut.contexts A vector of mutation contexts for fitting
#'
#' @description Calculate the contribution of each mutational signature in a sample, given a vector of mutation contexts.
#' The eps() and lsqnonneg() functions have been extracted and adapted from the pracma R package.
#'
#' @return A list containing, the absolute contribution of each signature (i.e. the number of mutations contributing to
#' each mutational signature), and the residual (amount of error from the fitting).
#' @export

## Machine epsilon (floating-point relative accuracy)
eps <- function(x = 1.0) {
   stopifnot(is.numeric(x))

   x <- max(abs(x))

   if (x < .Machine$double.xmin) {
      e <- .Machine$double.xmin
   } else {
      e <- 2^floor(log2(x)) * .Machine$double.eps
   }
   e
}

## Main function
fitToSignatures <- function(signature.profiles, mut.context.counts){

   if(!is.matrix(signature.profiles)){ stop("'signature.profiles' must be a matrix") }
   if(!is.vector(mut.context.counts)){ stop("'mut.context.counts' must be a vector") }

   m <- nrow(signature.profiles); n <- ncol(signature.profiles)
   if (m != length(mut.context.counts))
      stop("Arguments 'signature.profiles' and 'mut.context.counts' have nonconformable dimensions.")

   tol = 10 * eps() * norm(signature.profiles, type = "2") * (max(n, m) + 1)

   x  <- rep(0, n)             # initial point
   P  <- logical(n); Z <- !P   # non-active / active columns

   resid <- mut.context.counts - signature.profiles %*% x
   w <- t(signature.profiles) %*% resid
   wz <- numeric(n)

   # iteration parameters
   outeriter <- 0; it <- 0
   itmax <- 3 * n; exitflag <- 1

   while (any(Z) && any(w[Z] > tol)) {
      outeriter <- outeriter + 1
      z <- numeric(n)
      wz <- rep(-Inf, n)
      wz[Z] <- w[Z]
      im <- which.max(wz)
      P[im] <- TRUE; Z[im] <- FALSE
      z[P] <- qr.solve(signature.profiles[, P], mut.context.counts)

      while (any(z[P] <= 0)) {
         it <- it + 1
         if (it > itmax) stop("Iteration count exceeded.")

         Q <- (z <= 0) & P
         alpha <- min(x[Q] / (x[Q] - z[Q]))
         x <- x + alpha*(z - x)
         Z <- ((abs(x) < tol) & P) | Z
         P <- !Z
         z <- numeric(n)
         z[P] <- qr.solve(signature.profiles[, P], mut.context.counts)
      }
      x <- z
      resid <- mut.context.counts - signature.profiles %*% x
      w <- t(signature.profiles) %*% resid
   }
   return(list(x = x, resid.norm = sum(resid*resid)))
}
