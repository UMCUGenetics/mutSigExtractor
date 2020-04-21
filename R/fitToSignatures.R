#' Linear least-squares (non-negative) fitting
#'
#' @description Calculate the contribution of each mutational signature in a sample, given a vector of mutation contexts.
#' The eps() and lsqnonneg() functions have been extracted and adapted from the pracma R package.
#'
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are the mutation
#' contexts and the columns are  the mutational signatures.
#' @param mut.contexts A vector of mutation contexts for fitting
#' @param verbose Show messages?
#'
#' @return If vector is provided to mut.contexts, a vector returned containing the the absolute
#' contribution of each signature (i.e. the number of mutations contributing to each mutational
#' signature). If a matrix is provided, a matrix of absolute contributions is returned
#' @export
fitToSignatures <- function(signature.profiles, mut.context.counts, verbose=F){

   ## Machine epsilon (Get floating-point relative accuracy)
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

   ## Least squares fitting algorithm
   lsqnonneg <- function(signature.profiles, mut.context.counts){

      m <- nrow(signature.profiles); n <- ncol(signature.profiles)
      if (m != length(mut.context.counts))
         stop("Number of values in mut.context.counts does not match number of signatures in signature.profiles")

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
      return(x)
   }

   signature.profiles <- as.matrix(signature.profiles) ## Force signature.profiles into matrix if is dataframe

   ## Return
   if(is.vector(mut.context.counts)){
      if(verbose){ message('Fitting input vector to signature profiles...') }
      out <- lsqnonneg(signature.profiles, mut.context.counts)
      names(out) <- colnames(signature.profiles)
   }

   if(is.matrix(mut.context.counts) | is.data.frame(mut.context.counts)){

      mut.context.counts <- as.matrix(mut.context.counts)
      if( !identical(colnames(mut.context.counts), rownames(signature.profiles)) ){
         message('Contexts of mut context matrix and signature profiles do not match. Attempting to select revelant contexts...')
         mut.context.counts <- mut.context.counts[,rownames(signature.profiles)]
      }

      if(verbose){ message('Fitting input matrix to signature profiles...') }
      out <- t(apply(mut.context.counts, 1, function(i){
         lsqnonneg(signature.profiles, i)
      }))
      colnames(out) <- colnames(signature.profiles)
      rownames(out) <- rownames(mut.context.counts)
   }
   return(out)
}
