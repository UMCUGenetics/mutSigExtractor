% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detHomology.R
\name{nBasesMH}
\alias{nBasesMH}
\title{Calculate the number of bases that are homologous to the 3' flanking sequence}
\usage{
nBasesMH(indel.seq, flank.seq, indel.len)
}
\arguments{
\item{indel.seq}{The indel sequence as a character}

\item{flank.seq}{The flanking sequence as a character}
}
\value{
An integer stating the number of bases in microhomology
}
\description{
Helper function for extractSigsIndel(). Scans a maximum of 1 indel length from
5' to 3' in the flanking sequence to detect identical bases. Stops prematurely if a
non-identical base is found.

DSBs can be repaired using 3' microhomology, which can be present relative to either the +
or - strand. Therefore, both left and right flanking sequences of the indel need to be
considered. When calculating left flanking homology (i.e. homology in the 3' direction for
the - strand), the reverse complement of the indel sequence and flanking sequence need to
be taken. However, the reverse can be taken for the calculation to save computation.
}
