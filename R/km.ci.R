#' Confidence Intervals for the Kaplan-Meier Estimator.
#' 
#' Computes pointwise and simultaneous confidence intervals for the
#' Kaplan-Meier estimator.
#' 
#' A simulation study showed, that three confidence intervals produce
#' satisfying confidence limits. One is the "loglog" confidence interval, an
#' interval which is based on the log of the hazard. The other competitive
#' confidence concept was introduced by Rothman (1978) and is using the
#' assumption that the survival estimator follows a binomial distribution.
#' Another good confidence concept was invented by Thomas and Grunkemeier
#' (1975) and is derived by minimizing the likelihood function under certain
#' constraints. Special thanks goes to Robert Gentleman for providing code for
#' the confidence interval by Thomas and Grunkemeier.
#' 
#' The confidence interval using Peto's variance can not be recommended since
#' it yields confidence limits outside the admissible range [0;1] as well as
#' the "linear" and the "log" (which is based on the logarithm of S(t)).
#' 
#' The function can produce simultaneous confidence bands, too. The
#' Hall-Wellner band (1980) and the Equal Precision band by Nair (1984)
#' together with their log-transformed counterpart. From all simultaneous
#' confidence intervals only the log-transformed Equal Precision "logep" band
#' can be recommended. The limits are computed according to the statistical
#' tables in Klein and Moeschberger (2002).
#' 
#' @param survi A survival object for which the new confidence limits should be
#' computed. This can be built using the "Surv" and the "survfit" function in
#' the R package "survival". "km.ci" modifies the confidence limits in this
#' object.
#' @param conf.level The level for a two-sided confidence interval on the
#' survival curve. Default is 0.95.
#' @param tl The lower time boundary for the simultaneous confidence limits. If
#' it is missing the smallest event time is used.
#' @param tu The upper time boundary for the simultaneous confidence limits. If
#' it is missing the largest event time is used.
#' @param method One of '"peto"', '"linear"', '"log"', "loglog"', '"rothman"',
#' "grunkemeier"', '"hall-wellner"', '"loghall"', "epband"', "logep"
#' @return a 'survfit' object;
#' 
#' see the help on 'survfit.object' for details.
#' @author Strobl, R.
#' @seealso \code{\link[survival]{survfit}},
#' \code{\link[survival]{print.survfit}}, \code{\link[survival]{plot.survfit}},
#' \code{\link[survival]{lines.survfit}},
#' \code{\link[survival]{summary.survfit}},
#' \code{\link[survival]{survfit.object}}, \code{\link[survival]{coxph}},
#' \code{\link[survival]{Surv}}, \code{\link[survival]{strata}}.
#' @references Strobl, R., Dirschedl, P. and Mansmann, U..  Comparison of
#' simultaneous and pointwise confidence intervals for survival functions.
#' (2005, submitted to Biom. J.).
#' @keywords survival
#' @examples
#' 
#' require(survival)
#' data(rectum.dat)
#' 
#' # fit a Kaplan-Meier and plot it
#' fit <- survfit(Surv(time, status) ~ 1, data=rectum.dat)
#' plot(fit)
#' fit2 <- km.ci(fit)
#' plot(fit2)
#' 
#' @export km.ci
"km.ci" <-
function(survi,conf.level=0.95, tl=NA, tu=NA, method="rothman")
{
    # This function can compute the most desirable confidence bands.
    # The method "log" is implemented as "log" in R survfit.
    # The method "loglog" is implemented as "log-log" in R survfit.
    # The method "linear" is called "plain" in R survfit.

    if(conf.level < 0 || conf.level > 1)
      stop("confidence level must be between 0 and 1")
    if (data.class(survi)!="survfit")
            stop("Survi must be a survival object")
    method <- match.arg(method,c( "peto", "linear", "log" ,"loglog", "rothman","grunkemeier",
                "epband", "logep", "hall-wellner","loghall"))

    if(method=="grunkemeier")
    {
        result <- grunk.all.fun(survi,1-conf.level)
        result$conf.type <- "Grunkemeier"
    }

    if(method=="linear")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$linear$lower
        result$upper <- cf$linear$upper
        result$conf.type <- "Linear"
    }

    if(method=="rothman")
    {
        result <- rothman.fun(survi,conf.level)$surv.object
        result$conf.type <- "Rothman"
    }
    if(method=="peto")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$peto$lower
        result$upper <- cf$peto$upper
        result$conf.type <- "Peto"
    }
    if(method=="log")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$greenwood$lower
        result$upper <- cf$greenwood$upper
        result$conf.type <- "Log"
    }
    if(method=="loglog")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$log$lower
        result$upper <- cf$log$upper
        result$conf.type <- "Log-Log"
    }
    if(method=="hall-wellner")
    {
        result <- hall.wellner.fun(survi, tl=tl, tu=tu, conf.lev=conf.level)
        result$conf.type <- "Hall-Wellner"
    }
    if(method=="loghall")
    {
        result <- hall.wellner.fun(survi,tl=tl, tu=tu, method="log", conf.lev=conf.level)
        result$conf.type <- "Log(Hall-Wellner)"
    }

    if(method=="epband")
    {
        result <- epband.fun(survi, tl=tl, tu=tu, conf.lev=conf.level)
        result$conf.type <- "Equal Precision"
    }
    if(method=="logep")
    {
        result <- epband.fun(survi, tl=tl, tu=tu, method="log",conf.lev=conf.level)
        result$conf.type <- "Log(Equal Precision)"
    }
    return(result)
}

