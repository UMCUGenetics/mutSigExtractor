####################################################################################################
#' Linear least-squares (non-negative) fitting
#'
#' @rdname fitToSignatures
#'
#' @description Calculate the contribution of each mutational signature in a sample, given a vector
#' of mutation contexts. This function is an adaptation of the lsqnonneg() function in the pracma
#' R package.
#'
#' When method=='lsq' fitToSignatures() simplify performs the least-squares fitting, i.e. find a
#' linear non-negative combination of mutation signatures that reconstructs the mutation matrix
#'
#' fitToSignatures() is however prone to overfitting. When method=='strict', this is solved by
#' removing the signature with the lowest contribution and least-squares fitting is repeated. This
#' is done in an iterative fashion. Each time the cosine distance between the original and
#' reconstructed profile is calculated. Iterations are stopped when cosine distance > max/delta. The
#' second-last set of signatures is then returned.
#'
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#' @param mut.contexts A vector of mutation contexts for fitting
#' @param verbose Show messages?
#'
#' @return If vector is provided to mut.contexts, a vector returned containing the the absolute
#' contribution of each signature (i.e. the number of mutations contributing to each mutational
#' signature). If a matrix is provided, a matrix of absolute contributions is returned
#' @export
#'
fitToSignatures <- function (mut.context.counts, ...) {
   UseMethod("fitToSignatures", mut.context.counts)
}

## Machine epsilon (floating-point relative accuracy). Adapted from eps() function from pracma R package
EPS <- (function(x=1.0) {
   x <- max(abs(x))

   if (x < .Machine$double.xmin) {
      .Machine$double.xmin
   } else {
      2^floor(log2(x)) * .Machine$double.eps
   }
})()


cosSim <- function(x, y) { x %*% y / (sqrt(x %*% x) * sqrt(y %*% y)) }

#' @rdname fitToSignatures
#' @method fitToSignatures default
fitToSignatures.default <- function(
   mut.context.counts, signature.profiles,
   method='lsq', max.delta=0.01,
   verbose=F
){
   if(F){
      signature.profiles=SBS_SIGNATURE_PROFILES_V3
      method='strict'
      max.delta=0.01
   }

   ## Ripped least squares fitting function from pracma
   lsqnonneg <- function(mut.context.counts, signature.profiles){
      m <- nrow(signature.profiles)
      n <- ncol(signature.profiles)

      if (m != length(mut.context.counts)){
         stop(
            "No. of contexts in mut.context.counts (",ncol(mut.context.counts),")",
            "does not match no. of contexts in signature.profiles (",m,")"
         )
      }

      tol <- 10 * EPS * norm(signature.profiles, type = "2") * (max(n, m) + 1)

      x <- rep(0, n)             # initial point
      P <- logical(n); Z <- !P   # non-active / active columns

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

      names(x) <- colnames(signature.profiles)

      return(x)
   }

   ## Different fitting methods
   if(method=='lsq'){
      out <- lsqnonneg(mut.context.counts, signature.profiles)
   } else if(method=='strict'){
      v_init <- lsqnonneg(mut.context.counts, signature.profiles)

      v_template <- structure(
         rep(0,ncol(signature.profiles)),
         names=colnames(signature.profiles)
      )

      v1 <- v_init
      cossim1 <- 1
      for(i in 1:(ncol(signature.profiles)-1)){
         ## Remove the lowest contributing signature and signatures with 0 contribution
         v1_modified <- v1
         v1_modified[v1_modified==0] <- Inf ## Ensures which.min() skips 0 values

         excl_indexes <- unique(c(
            which.min(v1_modified),
            which(v1==0)
         ))
         signature_profiles_2 <- signature.profiles[,-excl_indexes, drop=F]

         if(ncol(signature_profiles_2)==0){ break } ## Exception when v1 only has contribution in one signature

         ## Do fit again
         v2 <- v_template ## Using the template ensures v2 has the same contexts as v1
         v2_pre <- lsqnonneg(mut.context.counts, signature_profiles_2)
         v2[names(v2_pre)] <- v2_pre
         #rbind(v1,v2)

         ## Calculate difference in cosine similarity
         cossim2 <- cosSim(v1,v2)[1]
         delta <- cossim1 - cossim2
         #print(paste0(cossim2,'_',delta))
         if(delta >= max.delta){ break }

         v1 <- v2
         cossim1 <- cossim2
      }

      out <- v1
   } else {
      stop("Method must be 'lsq' or 'strict'")
   }

   names(out) <- colnames(signature.profiles)
   return(out)
}

#' @rdname fitToSignatures
#' @method fitToSignatures matrix
fitToSignatures.matrix <- function(mut.context.counts, signature.profiles, verbose=F, ...){

   mut.context.counts <- as.matrix(mut.context.counts)

   if( !identical(colnames(mut.context.counts), rownames(signature.profiles)) ){
      warning('Contexts of mut.context.counts and signature.profiles do not match. Attempting to select revelant contexts...')
      mut.context.counts <- mut.context.counts[,rownames(signature.profiles)]
   }

   if(verbose){
      message('Fitting input matrix to signature profiles...')
      counter <- 0
      pb <- txtProgressBar(max=nrow(mut.context.counts), style=2)
   }
   out <- apply(mut.context.counts, 1, function(i){
      if(verbose){ counter <<- counter + 1; setTxtProgressBar(pb, counter) }
      #if(verbose){ counter <<- counter + 1; print(rownames(mut.context.counts)[counter]) }
      fitToSignatures.default(i, signature.profiles, ...)
   })
   if(verbose){ message('\n') }
   out <- t(out)
   colnames(out) <- colnames(signature.profiles)
   rownames(out) <- rownames(mut.context.counts)

   return(out)
}

#' @rdname fitToSignatures
#' @method fitToSignatures data.frame
fitToSignatures.data.frame <- fitToSignatures.matrix







