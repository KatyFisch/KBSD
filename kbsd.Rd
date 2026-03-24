\name{kbsd}
\alias{kbsd}
\title{Kernel-based Sparsity Diagnostic}
\description{
A diagnostic to assess sparsity issues in causal inference estimation
}

\usage{
kbsd(data, int_data_list, disthalf_vec, col="darkgreen", type="Rfast",
                          minval_vec=NA, kernel="rbf.normal", plot.out=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame with the observed data used for estimation. 
}
  \item{int_data_list}{
A list of data.frames with the counterfactual data for each estimand considered.
}
  \item{disthalf_vec}{
A named vector mapping every variable in the data to a so called half distance, the distance along this one variable where an observation would be assumed to contain half as much information as a point exactly matching the point of interest.
}
  \item{col}{
  A color for the plot.
}
  \item{type}{
  Either "Rfast", "minval", "harmonicmean", or "default". Rfast is vectorized and computes the normal diagnostic based on the Gaussian kernel (ie the parameter kernel is ignored). minval and harmonicmean compute the corresponding variants of the diagnostic. Rfast needs the package Rfast to be loaded.
}
  \item{minval_vec}{
  A vector of minimum values used only if type="minval".
}
  \item{plot.out}{
If TRUE, the function returns a plot, if FALSE, the function returns a data frame including the diagnostic values. Default is TRUE. For the plot, the ggplot2 package needs to be loaded.
}
 }
  \value{
Depending on plot.out either:
  \item{data.frame}{The diagnostic values}
  \item{ggplot}{A modified boxplot graph for the different estimands. The whiskers end at the 5th and 95th perentile in order to provide more insights into the lower diagnostic values.}
}
  \references{
Katharina Ring, Michael Schomaker: "A Diagnostic to Find and Help Combat Positivity Issues -- with a Focus on Continuous Treatments", arxiv, 2025

}

  \author{
Michael Schomaker, Philipp Baumann, Han Bao, Katharina Ring, Christoph Wiederkehr.
}


  \seealso{

}
  \examples{
library(simcausal)
n=100

## CREATE DAG
D.base <- DAG.empty() + 
  node("L", distr = "rbern", prob = 0.3) +
  node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14)
Dset <- simcausal::set.DAG(D.base)

# get data for plotting
Odat <- simcausal::sim(DAG = Dset, n = n, rndseed = 1, verbose=F) 
Odat$ID <- NULL

Odat_int1 <- Odat
Odat_int1$A = 1

Odat_int2 <- Odat
Odat_int2$A = Odat_int2$A + 0.1

Odat_int3 <- Odat
Odat_int3$A = Odat_int3$A + 0.5


data_result <- kbsd(data=Odat, 
                             int_data_list=list(Odat_int1, Odat_int2, Odat_int3), 
                             disthalf_vec=c(L=0.5, A=0.1)
                             )
}
  \keyword{sparsity}
