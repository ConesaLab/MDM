\name{vis_comp_net_meas_boxplots2}
\alias{vis_comp_net_meas_boxplots2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Analyze and visualize network measure changes for Original, Without Unknown and Bootstrap networks
%%  ~~function to do ... ~~
}
\description{ Creates customized boxplot of network measure changes across 3 network types for all ranks from Phylum to Genus for a given environment. Single or combinations of measures and ranks may also be selected. The statistical significance can be shown using stars (p.signif) or actual p-value scores (p.format).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vis_comp_net_meas_boxplots2(df, met_name, meas_name_options = NULL, rank_level_options = NULL, p_label = c("p.signif", "p.format"), p_size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{ Dataframe of all scores per measure, per network type, and per rank for all nodes, created from comp_by_deleting_random_knowns_t_v3 function
%%     ~~Describe \code{df} here~~
}
  \item{met_name}{ Name of environment
%%     ~~Describe \code{met_name} here~~
}
  \item{meas_name_options}{ Network measures to compare. If left blank, will compare all network measures (degree, betweenness, and closeness). Otherwise, can select a single measure or combination of measures.
%%     ~~Describe \code{meas_name_options } here~~
}
  \item{rank_level_options}{ Classification ranks to compare. If left blank, will compare network measure changes for all ranks from Phylum to Genus. Otherwise, can select a single rank or combination of ranks
%%     ~~Describe \code{rank_level_options} here~~
}
  \item{p_label}{ Label of statistical significance. Can choose between 2 options - p.signif (using stars to denote significance) or p.format (using numerical p-value to denote significance)
%%     ~~Describe \code{p_label} here~~
}
  \item{p_size}{ Size of statistical significance indicators. Default is 2.
%%     ~~Describe \code{p_size} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns boxplot ggplot object comparing network measure changes across network types for the indicated rank levels.
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
function (df, met_name, meas_name_options = NULL, rank_level_options = NULL,
    p_label = c("p.signif", "p.format"), p_size = 2)
{
    df = df[df$rank_level != "Kingdom" & df$rank_level != "Species",
        ]
    if (is.null(meas_name_options) && is.null(rank_level_options)) {
        g1 <- ggplot(df, aes(type, data, color = type)) + facet_grid(measure ~
            rank_level, scales = "free") + geom_boxplot(outlier.shape = NA) +
            ylab("Centrality score") + ggpubr::stat_compare_means(label = p_label,
            size = p_size, comparisons = list(c("Original", "Without_unknown"),
                c("Without_unknown", "Bootstrap"), c("Bootstrap",
                  "Original"))) + theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
    }
    else if (is.null(meas_name_options) && !(is.null(rank_level_options))) {
        spec_df = df[df$rank_level \%in\% rank_level_options, ]
        g1 <- ggplot(spec_df, aes(type, data, color = type)) +
            facet_grid(measure ~ rank_level, scales = "free") +
            geom_boxplot(outlier.shape = NA) + theme(axis.title.x = element_blank(),
            axis.text.x = element_blank()) + ggpubr::stat_compare_means(label = p_label,
            size = p_size, comparisons = list(c("Original", "Without_unknown"),
                c("Without_unknown", "Bootstrap"), c("Bootstrap",
                  "Original")))
    }
    else if (!(is.null(meas_name_options)) && is.null(rank_level_options)) {
        spec_df = df[df$measure \%in\% meas_name_options, ]
        g1 <- ggplot(spec_df, aes(type, data, color = type)) +
            facet_grid(measure ~ rank_level, scales = "free") +
            geom_boxplot(outlier.shape = NA) + theme(axis.title.x = element_blank(),
            axis.text.x = element_blank()) + ggpubr::stat_compare_means(label = p_label,
            size = p_size, comparisons = list(c("Original", "Without_unknown"),
                c("Without_unknown", "Bootstrap"), c("Bootstrap",
                  "Original")))
    }
    else {
        spec_df = df[df$measure \%in\% meas_name_options & df$rank_level \%in\%
            rank_level_options, ]
        g1 <- ggplot(spec_df, aes(type, data, color = type)) +
            geom_boxplot(outlier.shape = NA)
        if (length(meas_name_options) > 1 && length(rank_level_options) >
            1) {
            g1 <- g1 + facet_grid(measure ~ rank_level, scales = "free") +
                ggpubr::stat_compare_means(label = p_label, size = p_size,
                  comparisons = list(c("Original", "Without_unknown"),
                    c("Without_unknown", "Bootstrap"), c("Bootstrap",
                      "Original"))) + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank())
        }
        else if (length(meas_name_options) > 1 && length(rank_level_options) ==
            1) {
            g1 <- g1 + facet_wrap(~measure, scales = "free") +
                ggpubr::stat_compare_means(label = p_label, size = p_size,
                  comparisons = list(c("Original", "Without_unknown"),
                    c("Without_unknown", "Bootstrap"), c("Bootstrap",
                      "Original"))) + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank())
        }
        else if (length(meas_name_options) == 1 && length(rank_level_options) >
            1) {
            g1 <- g1 + facet_wrap(~rank_level, scales = "free") +
                ggpubr::stat_compare_means(label = p_label, size = p_size,
                  comparisons = list(c("Original", "Without_unknown"),
                    c("Without_unknown", "Bootstrap"), c("Bootstrap",
                      "Original"))) + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank())
        }
        else if (length(meas_name_options) == 1 && length(rank_level_options) ==
            1) {
            g1 <- g1 + ylab(meas_name_options) + ggtitle(rank_level_options) +
                ggpubr::stat_compare_means(label = p_label, comparisons = list(c("Original",
                  "Without_unknown"), c("Without_unknown", "Bootstrap"),
                  c("Bootstrap", "Original"))) + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank())
        }
    }
    g1 <- g1 + ggtitle(paste(met_name, "Centrality Score Comparison",
        sep = " ")) + ylab("Centrality score") + theme(plot.title = element_text(hjust = 0.5))
    return(g1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
