\name{get_hub_plots_all_ranks}
\alias{get_hub_plots_all_ranks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create 5 network plots with nodes sized by hub score, each colored by taxonomic classification
%%  ~~function to do ... ~~
}
\description{ Networks are created for each rank from Phylum to Genus, coloring nodes by their taxonomic classification at that rank and sizing nodes by their hub score. MDM nodes are always featured in gray.Larger nodes indicate larger hub score value.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_hub_plots_all_ranks(orig_graph, orig_phylo, met_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_graph}{ Network including all taxa, produced from get_back_res_meeting_min_occ function
%%     ~~Describe \code{orig_graph} here~~
}
  \item{orig_phylo}{Phyloseq including all taxa meeting threshold, produced from get_back_res_meeting_min_occ function
%%     ~~Describe \code{orig_phylo} here~~
}
  \item{met_name}{ Name of environment
%%     ~~Describe \code{met_name} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns list of 5 ggplot objects visualizing networks at each rank from Phylum to Genus, resizing nodes by hub score.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{hub_plot1 }{Phylum Hub Network}
\item{hub_plot2 }{Class Hub Network}
\item{hub_plot3 }{Order Hub Network}
\item{hub_plot4 }{Family Hub Network}
\item{hub_plot5}{Genus Hub Network}
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
function (orig_graph, orig_phylo, met_name)
{
    all_plots_l = list()
    for (rank_num in seq_along(colnames(tax_table(orig_phylo))[2:6])) {
        print(rank_num)
        tax_rank_names = c("Phylum", "Class", "Order", "Family",
            "Genus")
        rank_name = colnames(tax_table(orig_phylo))[2:6][rank_num]
        print(rank_name)
        rank_name2 = tax_rank_names[rank_num]
        print(rank_name2)
        g_net_plot <- phyloseq::plot_network(orig_graph, orig_phylo,
            type = "taxa", color = rank_name, label = NULL, point_size = hub_score(orig_graph)$vector *
                10, title = paste(met_name, rank_name2, sep = " ")) +
            theme(legend.position = "none")
        all_plots_l[[rank_name]] = g_net_plot
    }
    return(all_plots_l)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
