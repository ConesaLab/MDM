\name{get_cclasso_info}
\alias{get_cclasso_info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get graph and networks using cclasso to measure correlation between taxa co-occurrence
%%  ~~function to do ... ~~
}
\description{ Measure taxa co-occurrence relationships using cclasso method and construct networks of taxa with significant direct co-occurrence relationships by this method.An igraph object to visualize networks or calculate network measures is returned as output. Class network and hub network plots are also plotted and returned as output.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_cclasso_info(orig_phylo, met_name)
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
1. cclasso igraph. 2. List of 2 plots: Hub and Class network plots
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{cclasso_igraph}{Igraph of cclasso network}
 \item{Plot_list}{List of 2 plots: 1. Class network plot, 2. Hub network plot}
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
    hy_cclasso_m <- as.matrix(t(otu_table(orig_phylo)))
    hy_test_cclasso_f <- cclasso(hy_cclasso_m, counts = TRUE)
    hy_cclasso_corr_m <- as.matrix(hy_test_cclasso_f$cor_w)
    rownames(hy_cclasso_corr_m) = colnames(hy_cclasso_corr_m) = rownames(otu_table(orig_phylo))
    hy_cclasso_pval_m <- as.matrix(hy_test_cclasso_f$p_vals)
    dim(hy_cclasso_pval_m)
    rownames(hy_cclasso_pval_m) = colnames(hy_cclasso_pval_m) = rownames(otu_table(orig_phylo))
    hy_cclasso_cormat1 <- ifelse(hy_cclasso_pval_m <= 0.05, hy_cclasso_corr_m,
        0)
    diag(hy_cclasso_cormat1) <- 0
    hy_cclasso_cormat1_corr_thresh_0.6 <- abs(hy_cclasso_cormat1) >=
        0.6
    hy_cclasso_cormat1_adjmat <- Matrix::Matrix(hy_cclasso_cormat1_corr_thresh_0.6,
        sparse = TRUE)
    library(SpiecEasi)
    hy_cclasso_igraph <- igraph::adj2igraph(hy_cclasso_cormat1_adjmat,
        vertex.attr = list(name = taxa_names(orig_phylo)))
    met_cclasso_graph = hy_cclasso_igraph
    plot_a <- phyloseq::plot_network(hy_cclasso_igraph, orig_phylo,
        type = "taxa", color = "Class", label = NULL, title = paste(met_name,
            "CCLASSO network", sep = " ")) + theme(legend.position = "none")
    plot_b <- phyloseq::plot_network(hy_cclasso_igraph, orig_phylo,
        type = "taxa", color = "Genus", label = NULL, title = paste(met_name,
            "CCLASSO hub network", sep = " "), point_size = hub_score(hy_cclasso_igraph)$vector *
            10) + theme(legend.position = "none")
    all_plots <- list(plot_a, plot_b)
    all_info <- list(met_cclasso_graph, all_plots)
    return(all_info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
