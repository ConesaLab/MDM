\name{get_pearson_info}
\alias{get_pearson_info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get graph and networks using pearson to measure correlation between taxa co-occurrence
%%  ~~function to do ... ~~
}
\description{ Measure taxa co-occurrence relationships using pearson method and construct networks of taxa with significant direct co-occurrence relationships by this method.An igraph object to visualize networks or calculate network measures is returned as output. Class network and hub network plots are also plotted and returned as output.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_pearson_info(orig_phylo, met_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_phylo}{ Phyloseq including all taxa meeting threshold, produced from get_back_res_meeting_min_occ function
%%     ~~Describe \code{orig_phylo} here~~
}
  \item{met_name}{ Name of environment
%%     ~~Describe \code{met_name} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns List of 2 objects:
1. pearson igraph. 2. List of 2 plots: Hub and Class network plots
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{pearson_igraph}{Igraph of pearson network}
 \item{Plot_list}{List of 2 plots: 1. Class network plot, 2. Hub network plot}
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
function (orig_phylo, met_name)
{
    print(met_name)
    hy_test_pearson <- Hmisc::rcorr(as.matrix(t(otu_table(orig_phylo))),
        type = "pearson")
    hy_test_pearson_cormat <- as.matrix(hy_test_pearson$r)
    hy_test_pearson_pmat <- as.matrix(hy_test_pearson$P)
    hy_test_pearson_pmat[is.na(hy_test_pearson_pmat)] <- 1
    hy_test_pearson_cormat1 <- ifelse(hy_test_pearson_pmat <=
        0.05, hy_test_pearson_cormat, 0)
    hy_test_pearson_cormat2 <- abs(hy_test_pearson_cormat1) >=
        0.6
    hy_test_pearson_cormat2[1:10, 1:10]
    hy_test_pearson_adjmat <- Matrix::Matrix(hy_test_pearson_cormat2,
        sparse = TRUE)
    hy_test_pearson_igraph <- igraph::adj2igraph(hy_test_pearson_adjmat,
        vertex.attr = list(name = taxa_names(orig_phylo)))
    print(hy_test_pearson_igraph)
    plot_a <- phyloseq::plot_network(hy_test_pearson_igraph,
        orig_phylo, type = "taxa", color = "Class", label = NULL,
        title = paste(met_name, "Pearson network", sep = " ")) +
        theme(legend.position = "none")
    plot_b <- phyloseq::plot_network(hy_test_pearson_igraph,
        orig_phylo, type = "taxa", color = "Genus", label = NULL,
        title = paste(met_name, "Pearson hub network", sep = " "),
        point_size = hub_score(hy_test_pearson_igraph)$vector *
            10) + theme(legend.position = "none")
    all_plots <- list(plot_a, plot_b)
    all_info <- list(hy_test_pearson_igraph, all_plots)
    return(all_info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
