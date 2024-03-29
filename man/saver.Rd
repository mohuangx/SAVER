% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saver.R
\name{saver}
\alias{saver}
\title{Runs SAVER}
\usage{
saver(
  x,
  do.fast = TRUE,
  ncores = 1,
  size.factor = NULL,
  npred = NULL,
  pred.cells = NULL,
  pred.genes = NULL,
  pred.genes.only = FALSE,
  null.model = FALSE,
  mu = NULL,
  estimates.only = FALSE
)
}
\arguments{
\item{x}{An expression count matrix. The rows correspond to genes and
the columns correspond to cells. Can be sparse.}

\item{do.fast}{Approximates the prediction step. Default is TRUE.}

\item{ncores}{Number of cores to use. Default is 1.}

\item{size.factor}{Vector of cell size normalization factors.
If \code{x} is already normalized or normalization is not desired, use
\code{size.factor = 1}. Default uses mean library size normalization.}

\item{npred}{Number of genes for regression prediction. Selects the top
\code{npred} genes in terms of mean expression for regression prediction.
Default is all genes.}

\item{pred.cells}{Indices of cells to perform regression prediction.
Default is all cells.}

\item{pred.genes}{Indices of specific genes to perform regression
prediction. Overrides \code{npred}. Default is all genes.}

\item{pred.genes.only}{Return expression levels of only \code{pred.genes}.
Default is FALSE (returns expression levels of all genes).}

\item{null.model}{Whether to use mean gene expression as prediction.}

\item{mu}{Matrix of prior means.}

\item{estimates.only}{Only return SAVER estimates. Default is FALSE.}
}
\value{
If `estimates.only = TRUE`, then a matrix of SAVER estimates.

If `estimates.only = FALSE`, a list with the following components
\item{\code{estimate}}{Recovered (normalized) expression.}
\item{\code{se}}{Standard error of estimates.}
\item{\code{info}}{Information about dataset.}

The \code{info} element is a list with the following components:
\item{\code{size.factor}}{Size factor used for normalization.}
\item{\code{maxcor}}{Maximum absolute correlation for each gene. 2 if not
calculated}
\item{\code{lambda.max}}{Smallest value of lambda which gives the null
model.}
\item{\code{lambda.min}}{Value of lambda from which the prediction model is
used}
\item{\code{sd.cv}}{Difference in the number of standard deviations in
deviance between the model with lowest cross-validation error and the null
model}
\item{\code{pred.time}}{Time taken to generate predictions.}
\item{\code{var.time}}{Time taken to estimate variance.}
\item{\code{maxcor}}{Maximum absolute correlation cutoff used to determine
if a gene should be predicted.}
\item{\code{lambda.coefs}}{Coefficients for estimating lambda with lowest
cross-validation error.}
\item{\code{total.time}}{Total time for SAVER estimation.}
}
\description{
Recovers expression using the SAVER method.
}
\details{
The SAVER method starts by estimating the prior mean and variance for the
true expression level for each gene and cell. The prior mean is obtained
through predictions from a LASSO Poisson regression for each gene
implemented using the \code{glmnet} package. Then, the variance is estimated
through maximum likelihood assuming constant variance, Fano factor, or
coefficient of variation variance structure for each gene. The posterior
distribution is calculated and the posterior mean is reported as the SAVER
estimate.
}
\examples{
data("linnarsson")

\dontrun{
system.time(linnarsson_saver <- saver(linnarsson, ncores = 12))
}

# predictions for top 5 highly expressed genes
\dontrun{
saver2 <- saver(linnarsson, npred = 5)
}

# predictions for certain genes
\dontrun{
genes <- c("Thy1", "Mbp", "Stim2", "Psmc6", "Rps19")
genes.ind <- which(rownames(linnarsson) \%in\% genes)
saver3 <- saver(linnarsson, pred.genes = genes.ind)
}

# return only certain genes
\dontrun{
saver4 <- saver(linnarsson, pred.genes = genes.ind, pred.genes.only = TRUE)
}

}
