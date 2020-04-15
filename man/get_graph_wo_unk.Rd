\name{get_graph_wo_unk}
\alias{get_graph_wo_unk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create igraph object without MDM nodes at each rank
%%  ~~function to do ... ~~
}
\description{ For each taxonomic rank, all taxa designated as MDM are removed from the original network (igraph) using names of the MDM found by the met_wo_unk function. A new network (igraph) for each rank is produced, where only taxa with known taxonomic classification are present.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_graph_wo_unk(orig_graph, met_wo_unk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_graph}{ Network (graph) of all taxa, including MDM
%%     ~~Describe \code{orig_graph} here~~
}
  \item{met_wo_unk}{ List of 7 otu tables, produced from met_wo_unk function, each corresponding to the taxa with known classification at that rank. Each OTU table shows the taxa ids and counts for these known taxa.
%%     ~~Describe \code{met_wo_unk} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns List of 7 igraph objects.
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{wo_unk_graph1}{Network (igraph) of taxa with known Kingdom classification}
 \item{wo_unk_graph2 }{Network (igraph) of taxa with known Phylum classification}
 \item{wo_unk_graph3}{Network (igraph) of taxa with known Class classification}
 \item{wo_unk_graph4}{Network (igraph) of taxa with known Order classification'}
 \item{wo_unk_graph5}{Network (igraph) of taxa with known Family classification}
 \item{wo_unk_graph6}{Network (igraph) of taxa with known Genus classification}
 \item{wo_unk_graph7 }{Network (igraph) of taxa with known Species classification}

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
function (orig_graph, met_wo_unk)
{
    wo_unk_graph_l = list()
    for (name in names(met_wo_unk)) {
        print(name)
        tax_level = name
        wo_unk_graph <- igraph::delete_vertices(orig_graph, which(!names(V(orig_graph)) \%in\%
            rownames(met_wo_unk[[tax_level]])))
        print(wo_unk_graph)
        wo_unk_graph_l[[name]] = wo_unk_graph
    }
    return(wo_unk_graph_l)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
