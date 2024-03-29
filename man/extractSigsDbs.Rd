% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractSigsDbs.R
\name{extractSigsDbs}
\alias{extractSigsDbs}
\title{Extract doublet substitution signatures}
\usage{
extractSigsDbs(
  vcf.file = NULL,
  df = NULL,
  output = "contexts",
  sample.name = NULL,
  signature.profiles = DBS_SIGNATURE_PROFILES,
  ref.genome = NULL,
  verbose = F,
  ...
)
}
\arguments{
\item{vcf.file}{Path to the vcf file}

\item{df}{A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to
vcf.file}

\item{output}{Output the absolute signature contributions (default, 'signatures'), the
DBS contexts ('contexts'), or an annotated bed-like dataframe ('df')}

\item{sample.name}{If a character is provided, the header for the output matrix will be named to
this. If none is provided, the basename of the vcf file will be used.}

\item{signature.profiles}{A matrix containing the mutational signature profiles, where rows are
the mutation contexts and the columns are  the mutational signatures.}

\item{ref.genome}{Deprecated. Argument only kept for compatibility.}

\item{verbose}{Print progress messages?}
}
\value{
A 1-column matrix containing the context counts or signature contributions
}
\description{
Will output a 1-column matrix containing: (if output=='signatures') the absolute
signature contributions (i.e. the number of mutations contributing to each mutational signature),
or (if output=='contexts') the mutation contexts, or (if output=='df') a dataframe with each
mutation annotated by context
}
