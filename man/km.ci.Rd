\name{km.ci}
\alias{km.ci}

\title{Confidence intervals for the Kaplan-Meier estimator.}

\description{
 Computes pointwise and simultaneous confidence intervals for the Kaplan-Meier estimator.
}
\usage{
km.ci(survi, conf.level=0.95, tl=NA, tu=NA, method="rothman")
}

\arguments{
  \item{survi}{A survival object for which the new confidence limits should be computed. This can be
  built using the "Surv" and the "survfit" function in the R package "survival". "km.ci" modifies the
  confidence limits in this object.}
  \item{conf.level}{The level for a two-sided confidence interval on the
          survival curve. Default is 0.95. }
  \item{tl}{The lower time boundary for the simultaneous confidence limits. If
  it is missing the smallest event time is used.}
   \item{tu}{The upper time boundary for the simultaneous confidence limits. If
  it is missing the largest event time is used.}
    \item{method}{One of '"peto"', '"linear"', '"log"', "loglog"', 
  '"rothman"', "grunkemeier"', '"hall-wellner"', '"loghall"', "epband"', "logep"
  }
  }
\details{
A simulation study showed, that three confidence intervals produce
 satisfying confidence limits. One is the "loglog" confidence interval, an interval
which is based on the log of the hazard. The other competitive confidence concept
was introduced by Rothman (1978) and is using the assumption that the survival estimator
follows a binomial distribution. Another good confidence concept was invented by
 Thomas and Grunkemeier (1975) and is derived by minimizing the likelihood function
 under certain constraints. Special thanks goes to Robert Gentleman for providing code for
  the confidence interval by Thomas and Grunkemeier.
 
 
 The confidence interval using Peto's variance can not be recommended since it yields
 confidence limits outside the admissible range [0;1] as well as the "linear" and the
 "log" (which is based on the logarithm of S(t)). 

The function can produce simultaneous confidence bands, too. The Hall-Wellner band (1980)
and the Equal Precision band by Nair (1984) together with their log-transformed
counterpart. From all simultaneous confidence intervals only the log-transformed Equal Precision
"logep" band can be recommended. The limits are computed according to the statistical 
tables in Klein and Moeschberger (2002).

}
\value{
  a 'survfit' object; 
  
   see the help on 'survfit.object' for details.
}
\author{Strobl, R.}
\references{Strobl, R., Dirschedl, P. and Mansmann, U.. 
Comparison of simultaneous and pointwise confidence intervals for survival functions.
 (2005, submitted to Biom. J.).}

\seealso{\code{\link[survival]{survfit}}, \code{\link[survival]{print.survfit}}, 
  \code{\link[survival]{plot.survfit}}, \code{\link[survival]{lines.survfit}}, 
  \code{\link[survival]{summary.survfit}}, \code{\link[survival]{survfit.object}}, 
 \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}}, \code{\link[survival]{strata}}. }
\examples{
require(survival)
data(rectum.dat)

# fit a Kaplan-Meier and plot it
fit <- survfit(Surv(time, status) ~ 1, data=rectum.dat)
plot(fit)
fit2 <- km.ci(fit)
plot(fit2)
}
\keyword{survival}
