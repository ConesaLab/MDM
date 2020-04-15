\name{filterTaxonMatrix}
\alias{filterTaxonMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Filter taxa in an abundance matrix
%%  ~~function to do ... ~~
}
\description{ function from seqtime package, Discard taxa with less than the given minimum number of occurrences
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filterTaxonMatrix(x,minocc = 0, dependency = FALSE, keepSum = FALSE,
    return.filtered.indices = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ taxon abundance matrix, rows are taxa, columns are samples
%%     ~~Describe \code{x} here~~
}
\item{minocc}{minimum occurrence (minimum number of samples w/ non-zero taxon abundance)}
\item{dependency}{if true, remove all taxa with a slope above -0.5 or a non-linear slope in the periodogram in log-scale (samples are supposed to represent equidistant time points)}
\item{keepSum}{if keepSum is true,  the discarded rows are summed and the sum is added as a row with name: summed-nonfeat-rows}
\item{return.filtered.indices}{if true, return an object with the filtered abundance matrix in mat and the indices of removed taxa in the original matrix in filtered.indices}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{filtered abundance matrix
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
function (x, minocc=0,dependency = FALSE, keepSum = FALSE,
    return.filtered.indices = FALSE)
{
    toFilter = c()
    xcopy = x
    xcopy[xcopy > 0] = 1
    rowsums = apply(xcopy, 1, sum)
    toFilter = which(rowsums < minocc)
    if (dependency == TRUE) {
        nt = identifyNoisetypes(x, epsilon = 0.5)
        toKeep = c(nt$pink, nt$brown, nt$black)
        toFilter = c(toFilter, setdiff(c(1:nrow(x)), toKeep))
    }
    indices.tokeep = setdiff(c(1:nrow(x)), toFilter)
    if (keepSum == TRUE) {
        filtered = x[toFilter, ]
        x = x[indices.tokeep, ]
        rownames = rownames(x)
        sums.filtered = apply(filtered, 2, sum)
        x = rbind(x, sums.filtered)
        rownames = append(rownames, "summed-nonfeat-rows")
        rownames(x) = rownames
    }
    else {
        x = x[indices.tokeep, ]
    }
    if (return.filtered.indices == TRUE) {
        res = list(x, toFilter)
        names(res) = c("mat", "filtered.indices")
        return(res)
    }
    else {
        return(x)
    }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
