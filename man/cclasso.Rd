\name{cclasso}
\alias{cclasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Measure correlation by cclasso function
%%  ~~function to do ... ~~
}
\description{ CCLasso function, originally created by Fang Huaying and copied here
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cclasso(x, counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-04, 1), k_max = 20, n_boot = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ n x p Matrix (row/column is sample/variable)
  n samples & p compositional variables
%%     ~~Describe \code{x} here~~
}
  \item{counts}{ Is compositional data matrix a count matrix? Default: False
%%     ~~Describe \code{counts} here~~
}
  \item{pseudo}{ pseudo count if counts=TRUE. Default: 0.5
%%     ~~Describe \code{pseudo} here~~
}
  \item{k_cv}{ Number of k-fold cross validations. Default: 3
%%     ~~Describe \code{k_cv} here~~
}
  \item{lam_int}{ tuning parameter interval
%%     ~~Describe \code{lam_int} here~~
}
  \item{k_max}{ Maximum iterations for golden section method. Default:20
%%     ~~Describe \code{k_max} here~~
}
  \item{n_boot}{ Bootstrap times. Default:20
%%     ~~Describe \code{n_boot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns list of variance matrix, correlation matrix, p-values, lambda and all info
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-04,
    1), k_max = 20, n_boot = 20)
{
    n <- nrow(x)
    p <- ncol(x)
    if (counts) {
        x <- x + pseudo
        x <- x/rowSums(x)
    }
    x <- log(x)
    vx2 <- stats::var(x)
    rmean_vx2 <- rowMeans(vx2)
    wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) +
        mean(rmean_vx2))
    wd2 <- sqrt(wd)
    rho <- 1
    u_f <- eigen(diag(p) - 1/p)$vectors
    wd_u <- (t(u_f) \%*\% (wd * u_f))[-p, -p]
    wd_u_eig <- eigen(wd_u)
    d0_wd <- 1/((rep(wd_u_eig$values, each = p - 1) + wd_u_eig$values)/(2 *
        rho) + 1)
    u0_wd <- wd_u_eig$vectors
    sigma <- vx2
    lam_int2 <- log10(range(lam_int))
    a1 <- lam_int2[1]
    b1 <- lam_int2[2]
    lams <- NULL
    fvals <- NULL
    a2 <- a1 + 0.382 * (b1 - a1)
    b2 <- a1 + 0.618 * (b1 - a1)
    fb2 <- cv_loss_cclasso(lambda2 = 10^b2/rho, x = x, k_cv = k_cv,
        sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
        wd2 = wd2)
    lams <- c(lams, b2)
    fvals <- c(fvals, fb2$cv_loss)
    fa2 <- cv_loss_cclasso(lambda2 = 10^a2/rho, x = x, k_cv = k_cv,
        sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd,
        d0_wd = d0_wd, wd2 = wd2)
    lams <- c(lams, a2)
    fvals <- c(fvals, fa2$cv_loss)
    err_lam2 <- 0.1 * max(1, lam_int2)
    err_fval <- 1e-04
    err <- b1 - a1
    k <- 0
    while (err > err_lam2 && k < k_max) {
        fval_max <- max(fa2$cv_loss, fb2$cv_loss)
        if (fa2$cv_loss > fb2$cv_loss) {
            a1 <- a2
            a2 <- b2
            fa2 <- fb2
            b2 <- a1 + 0.618 * (b1 - a1)
            fb2 <- cv_loss_cclasso(lambda2 = 10^b2/rho, x = x,
                k_cv = k_cv, sigma = fa2$sigma, wd = wd, u_f = u_f,
                u0_wd = u0_wd, d0_wd = d0_wd, wd2 = wd2)
            lams <- c(lams, b2)
            fvals <- c(fvals, fb2$cv_loss)
        }
        else {
            b1 <- b2
            b2 <- a2
            fb2 <- fa2
            a2 <- a1 + 0.382 * (b1 - a1)
            fa2 <- cv_loss_cclasso(lambda2 = 10^a2/rho, x = x,
                k_cv = k_cv, sigma = fb2$sigma, wd = wd, u_f = u_f,
                u0_wd = u0_wd, d0_wd = d0_wd, wd2 = wd2)
            lams <- c(lams, a2)
            fvals <- c(fvals, fa2$cv_loss)
        }
        fval_min <- min(fa2$cv_loss, fb2$cv_loss)
        k <- k + 1
        err <- b1 - a1
        if (abs(fval_max - fval_min)/(1 + fval_min) <= err_fval) {
            break
        }
    }
    info_cv <- list(lams = lams, fvals = fvals, k = k + 2, lam_int = 10^c(a1,
        b1))
    if (a1 == lam_int2[1] || b1 == lam_int2[2]) {
        cat("WARNING:\n", "\tOptimal lambda is near boundary! ([",
            10^a1, ",", 10^b1, "])\n", sep = "")
    }
    lambda <- 10^((a2 + b2)/2)
    lambda2 <- lambda/rho
    info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2,
        n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd)
    return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w,
        p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
