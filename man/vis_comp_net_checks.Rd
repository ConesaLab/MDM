\name{vis_comp_net_checks}
\alias{vis_comp_net_checks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Analyze and visualize network measure changes, checking for sample or network tool bias
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vis_comp_net_checks(df, perc_sample_options = NULL, metric_options = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{perc_sample_options}{
%%     ~~Describe \code{perc_sample_options} here~~
}
  \item{metric_options}{
%%     ~~Describe \code{metric_options} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
function (df, perc_sample_options = NULL, metric_options = NULL)
{
    spec_df = df
    if (!is.null(perc_sample_options)) {
        spec_df_median <- spec_df \%>\% dplyr::group_by(type, rank_level,
            measure, percent_samples_included) \%>\% dplyr::summarise(Median = median(data)) \%>\%
            as.data.frame()
        g1 <- ggplot(spec_df_median, aes(rank_level, Median,
            color = type)) + geom_point(aes(shape = percent_samples_included)) +
            facet_wrap(~ measure, scales = "free") +
            theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    else if (!is.null(metric_options)) {
        spec_df_median <- spec_df \%>\% dplyr::group_by(type, rank_level,
            measure, metric) \%>\% dplyr::summarise(Median = median(data)) \%>\%
            as.data.frame()
        g1 <- ggplot(spec_df_median, aes(rank_level, Median,
            color = type)) + geom_point(aes(shape = metric)) +
            facet_wrap(~measure, scales = "free") + theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1,
                size = 7))
    }
    return(g1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
