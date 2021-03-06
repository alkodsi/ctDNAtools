% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_read_counts.R
\name{get_read_counts}
\alias{get_read_counts}
\title{Gets read counts for a specific locus in the genome}
\usage{
get_read_counts(
  chr,
  pos,
  bam,
  tag = "",
  min_base_quality = 20,
  max_depth = 1e+05,
  min_mapq = 30
)
}
\arguments{
\item{chr}{chromosome name}

\item{pos}{genomic coordinate}

\item{bam}{path to bam file}

\item{tag}{the RG tag if the bam has more than one sample}

\item{min_base_quality}{minimum base quality for a read to be counted}

\item{max_depth}{maximum depth above which sampling will happen}

\item{min_mapq}{the minimum mapping quality for a read to be counted}
}
\value{
a list, number of reads for each of the four basepairs
}
\description{
Uses samtools pileup to get the read counts for each base in the genomic position specified
}
\keyword{internal}
