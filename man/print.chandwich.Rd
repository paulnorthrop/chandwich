% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_summary_print.R
\name{print.chandwich}
\alias{print.chandwich}
\title{Print method for objects of class "chandwich"}
\usage{
\method{print}{chandwich}(x, ...)
}
\arguments{
\item{x}{an object of class "chandwich", a result of a call to
\code{\link{adjust_loglik}}.}

\item{...}{Additional optional arguments. At present no optional
arguments are used.}
}
\value{
The argument \code{x}, invisibly, as for all
  \code{\link[base]{print}} methods.
}
\description{
\code{print} method for class "chandwich".
}
\details{
Just prints the original call to \code{\link{adjust_loglik}}
  and a character vector giving the names of the attributes
  (produced using \code{ls(attributes(x))}) to the function returned
  from \code{\link{adjust_loglik}}.
  To view an individual attribute called \code{att_name} use
  \code{attr(x, "att_name")} or \code{attributes(x)$att_name}.
}
\seealso{
\code{\link{summary.chandwich}}: \code{summary} method for
  class "chandwich".

\code{\link{adjust_loglik}} to adjust a user-supplied
  loglikelihood.
}
