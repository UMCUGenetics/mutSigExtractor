####################################################################################################
#' Linear least-squares (non-negative) fitting
#'
#' @rdname fitToSignatures
#'
#' @description Calculate the contribution of each mutational signature in a sample, given a vector
#' of mutation contexts. This is done by least-squares fitting, i.e. find a linear non-negative
#' combination of mutation signatures that optimally reconstructs the mutation matrix
#'
#' fitToSignatures() is a slow R implementation of least-squares fitting based on the lsqnonneg()
#' function in the pracma from the R package.
#'
#' fitToSignaturesFast() is a wrapper for NNLM::nnlm(), a least-squares fitting function implemented
#' in C++.
#'
#' Only performing least-squares fitting is however potentially prone to overfitting.
#' fitToSignaturesFastStrict() tries to solve this problem by performing an initial fit,
#' then iteratively removing the signature with the lowest contribution and with fitting then
#' repeated. Each time the cosine distance between the original and reconstructed mutation context
#' profile is calculated. Iterations are stopped when cosine distance > max.delta. The second-last
#' set of signatures is then returned.
#'
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are the mutational signatures.
#' @param mut.contexts A vector of mutation contexts for fitting
#' @param max.delta See description.
#' @param detailed.output Only for fitToSignaturesFastStrict(). Also return results from the
#' iterative fitting? Includes: cosine similarity between the original and reconstructed mutation
#' context; signatures removed at each iteration. Useful for plotting fitting performance.
#' @param verbose Show progress messages?
#'
#' @return If vector is provided to mut.contexts, a vector returned containing the the absolute
#' contribution of each signature (i.e. the number of mutations contributing to each mutational
#' signature). If a matrix is provided, a matrix of absolute contributions is returned
#' @export
#'
fitToSignatures <- function (mut.context.counts, ...) {
   UseMethod("fitToSignatures", mut.context.counts)
}

## Fast implementations via NNLM ###################################################################
#' @rdname fitToSignatures
fitToSignaturesFast <- function(mut.context.counts, signature.profiles, verbose=F){
   # if(F){
   #    mut.context.counts=contexts$indel
   #    mut.context.counts=contexts$indel[1,]
   #    signature.profiles=INDEL_SIGNATURE_PROFILES
   # }

   ## Checks --------------------------------
   if(is.vector(mut.context.counts)){ mut.context.counts <- t(mut.context.counts) }
   mut.context.counts <- as.matrix(mut.context.counts)
   if(!is.numeric(mut.context.counts)){ stop('mut.context.counts must be a numeric matrix or dataframe') }

   signature.profiles <- as.matrix(signature.profiles)
   if(!is.numeric(mut.context.counts)){ stop('signature.profiles must be a numeric matrix or dataframe') }

   if( nrow(mut.context.counts)!=nrow(signature.profiles) ){
      stop(
         "No. of contexts in mut.context.counts (",nrow(mut.context.counts),")",
         "does not match no. of contexts in signature.profiles (",nrow(signature.profiles),")"
      )
   }

   if( !(identical(rownames(mut.context.counts), rownames(signature.profiles))) ){
      warning("Context names of mut.context.counts and signature.profiles do not match. Fitting may not be correct.")
   }

   ## Main --------------------------------
   sigs <- NNLM::nnlm(signature.profiles, mut.context.counts)$coefficients
   sigs <- t(sigs)

   return(sigs)
}

#' @rdname fitToSignatures
fitToSignaturesFastStrict <- function(
   mut.context.counts, signature.profiles, max.delta=0.004,
   detailed.output=F, verbose=F
){

   if(F){
      mut.context.counts=contexts$snv[1:2,]
      signature.profiles=SBS_SIGNATURE_PROFILES_V3
      max.delta=0.004
      verbose=T
   }

   ## Checks --------------------------------
   require(NNLM)

   if(is.vector(mut.context.counts)){ mut.context.counts <- t(mut.context.counts) }
   mut.context.counts <- as.matrix(mut.context.counts)
   if(!is.numeric(mut.context.counts)){ stop('mut.context.counts must be a numeric matrix or dataframe') }

   signature.profiles <- as.matrix(signature.profiles)
   if(!is.numeric(mut.context.counts)){ stop('signature.profiles must be a numeric matrix or dataframe') }

   if( nrow(mut.context.counts)!=nrow(signature.profiles) ){
      stop(
         "No. of contexts in mut.context.counts (",nrow(mut.context.counts),")",
         "does not match no. of contexts in signature.profiles (",nrow(signature.profiles),")"
      )
   }

   if( !(identical(rownames(mut.context.counts), rownames(signature.profiles))) ){
      warning("Context names of mut.context.counts and signature.profiles do not match. Fitting may not be correct.")
   }

   ## --------------------------------
   if(verbose){ message('Performing first fit...') }

   mut.context.counts <- t(mut.context.counts)
   fit_init <- nnlm(signature.profiles, mut.context.counts)$coefficients
   sig_pres <- rowSums(fit_init) != 0
   my_signatures_total <- signature.profiles[, sig_pres, drop=FALSE]

   ## --------------------------------
   if(verbose){ message('Performing signature selection per sample...') }

   cosSim <- function(x, y) { sum(x*y)/sqrt(sum(x^2)*sum(y^2)) }

   n_sigs <- ncol(my_signatures_total)
   n_samples <- ncol(mut.context.counts)

   l_contribs <- list()
   l_sims <- list()
   l_removed_sigs <- list()

   if(verbose==2){ pb <- txtProgressBar(max=n_samples, style=3) }

   for (i in 1:n_samples) {
      #i=1

      if(verbose==1){ message(' [',i,'/',n_samples,'] ', colnames(mut.context.counts)[i]) }
      if(verbose==2){ setTxtProgressBar(pb, i) }

      my_signatures <- my_signatures_total
      mut_mat_sample <- mut.context.counts[, i, drop=FALSE]

      ## Fit again
      fit_res <- list()
      fit_res$contribution <- nnlm(my_signatures, mut_mat_sample)$coefficients
      fit_res$reconstructed <- my_signatures %*% fit_res$contribution

      sim <- cosSim(fit_res$reconstructed, mut_mat_sample)

      ## Keep track of the cosine similarity and which signatures are removed.
      sims <- rep(NA, n_sigs); sims[[1]] <- sim
      removed_sigs <- rep(NA, n_sigs)

      ## Sequentially remove the signature with the lowest contribution
      for (j in 2:n_sigs) {
         #j=2

         # Remove signature with the weakest relative contribution
         contri_order <- order(fit_res$contribution[,1] / sum(fit_res$contribution[,1]))
         weakest_sig_index <- contri_order[1]
         weakest_sig <- colnames(my_signatures)[weakest_sig_index]
         removed_sigs[[j]] <- weakest_sig
         signatures_sel <- my_signatures[, -weakest_sig_index, drop=FALSE]

         # Fit with new signature selection
         fit_res <- list()
         fit_res$contribution <- nnlm(signatures_sel, mut_mat_sample)$coefficients
         fit_res$reconstructed <- signatures_sel %*% fit_res$contribution

         sim_new <- cosSim(fit_res$reconstructed, mut_mat_sample)

         if(is.na(sim_new)){
            sim_new <- 0
            # if(verbose){
            #    warning("New similarity between the original and the reconstructed
            #             spectra after the removal of a signature was NaN.
            #             It has been converted into a 0.
            #             This happened with the following fit_res:")
            #    print(fit_res)
            # }
         }
         sims[[j]] <- sim_new

         # Check if the loss in cosine similarity between the original vs reconstructed after removing the signature is below the cutoff.
         if(sim-sim_new <= max.delta){
            my_signatures <- signatures_sel
            sim <- sim_new
         } else {
            break
         }
      }

      ## Fit with the final set of signatures
      ## Fill in 0 for absent signatures
      contrib_pre <- NNLM::nnlm(my_signatures, mut_mat_sample)$coefficients[,1]
      contrib <- structure(
         rep(0, ncol(signature.profiles)),
         names=colnames(signature.profiles)
      )
      contrib[names(contrib_pre)] <- contrib_pre

      l_contribs[[i]] <- contrib
      l_sims[[i]] <- sims
      l_removed_sigs[[i]] <- removed_sigs
   }

   ## --------------------------------
   if(verbose==2){ message('\n') }
   if(verbose){ message('Returning output...') }

   m_contribs <- do.call(rbind, l_contribs)
   rownames(m_contribs) <- colnames(mut.context.counts)

   if(!detailed.output){ return(m_contribs) }

   m_sim <- do.call(rbind, l_sims)
   rownames(m_sim) <- colnames(mut.context.counts)

   m_removed_sigs <- do.call(rbind, l_removed_sigs)
   rownames(m_removed_sigs) <- colnames(mut.context.counts)

   list(
      sig_contrib=m_contribs,
      removed_sigs_common=names(sig_pres)[ !sig_pres ],
      removed_sigs_iter=m_removed_sigs,
      cos_sims=m_sim
   )
}


## Slow R implementations ##########################################################################
## Machine epsilon (floating-point relative accuracy). Adapted from eps() function from pracma R package
EPS <- (function(x=1.0) {
   x <- max(abs(x))

   if (x < .Machine$double.xmin) {
      .Machine$double.xmin
   } else {
      2^floor(log2(x)) * .Machine$double.eps
   }
})()

#' Non-negative least squares fitting
#'
#' @description Ripped least squares fitting function from pracma
#'
#' @param mut.context.counts A numeric vector of mutation context counts
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are the mutational signatures.
#'
#' @return A vector returned containing the the absolute contribution of each signature
#'
lsqnonneg <- function(mut.context.counts, signature.profiles){
   m <- nrow(signature.profiles)
   n <- ncol(signature.profiles)

   if (m != length(mut.context.counts)){
      stop(
         "No. of contexts in mut.context.counts (",length(mut.context.counts),")",
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

#' @rdname fitToSignatures
#' @method fitToSignatures default
fitToSignatures.default <- function(
   mut.context.counts, signature.profiles,
   method='lsq', max.delta=0.01,
   verbose=F
){
   # if(F){
   #    signature.profiles=SBS_SIGNATURE_PROFILES_V3
   #    method='strict'
   #    max.delta=0.01
   # }

   out <- lsqnonneg(mut.context.counts, signature.profiles)
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







