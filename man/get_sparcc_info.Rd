\name{get_sparcc_info}
\alias{get_sparcc_info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get graph and networks using sparcc to measure correlation between taxa co-occurrence
%%  ~~function to do ... ~~
}
\description{ Measure taxa co-occurrence relationships using sparcc method and construct networks of taxa with significant direct co-occurrence relationships by this method.An igraph object to visualize networks or calculate network measures is returned as output. Class network and hub network plots are also plotted and returned as output.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_sparcc_info(orig_phylo, met_name)
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
\value{ Returns List of 2 objects: 1. sparcc igraph for network. 2. Net and hub plots
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{sparcc_igraph}{Igraph of sparcc network}
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
    otu_matrix <- as(otu_table(orig_phylo), "matrix")
    print(dim(otu_matrix))
    otu_matrix_t <- t(otu_matrix)
    print(dim(otu_matrix_t))
    orig_sparcc <- SpiecEasi::sparcc(otu_matrix_t)
    rownames(orig_sparcc$Cor) = colnames(orig_sparcc$Cor) = taxa_names(otu_table(orig_phylo))
    val = 0.3
    new_sparcc.graph <- abs(orig_sparcc$Cor) >= val
    diag(new_sparcc.graph) <- 0
    new_sparcc.graph <- Matrix::Matrix(new_sparcc.graph, sparse = TRUE)
    new_sparcc_igraph <- igraph::adj2igraph(new_sparcc.graph)
    print(new_sparcc_igraph)
    V(new_sparcc_igraph)$name <- rownames(orig_sparcc$Cor)
    print(V(new_sparcc_igraph)$name)
    plot_a <- phyloseq::plot_network(new_sparcc_igraph, orig_phylo,
        type = "taxa", color = "Class", label = NULL, title = paste(met_name,
            "SparCC network", sep = " ")) + theme(legend.position = "none")
    plot_b <- phyloseq::plot_network(new_sparcc_igraph, orig_phylo,
        type = "taxa", color = "Genus", label = NULL, title = paste(met_name,
            "SparCC hub network", sep = " "), point_size = hub_score(new_sparcc_igraph)$vector *
            10) + theme(legend.position = "none")
    all_plots <- list(plot_a, plot_b)
    all_info <- list(new_sparcc_igraph, all_plots)
    return(all_info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
