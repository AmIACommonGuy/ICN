% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NICE_fast_cut.R
\name{NICE_fast}
\alias{NICE_fast}
\title{NICE_fast}
\usage{
NICE_fast(corr_m, thres, cutter = 1)
}
\arguments{
\item{corr_m}{correlation matrix}

\item{thres}{an integer threshold value between (0, 1)}

\item{cutter}{A sequence jump during iteration of K, default 1.}
}
\value{
Cindx: the cluster index of every non-isolated node

CID: the cluster index of every cluster in a power
descending order. i.e. CID(1) will be the cluster index of the
most concentrated cluster

Clist: the reordered node index, nodes in the same
cluster are permuted together.

K_selected: number of communities identified

t.overall: Time counter for NICE

Prp_net: evaluation index of K
}
\description{
NICE_fast
}
