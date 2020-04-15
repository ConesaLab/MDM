\name{get_degree_df_wo_unk}
\alias{get_degree_df_wo_unk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get degree df without unk
%%  ~~function to do ... ~~
}
\description{ Calculate degree, betweenness, and closeness scores for each node in network (igraph), excluding MDM taxa at each rank, and store results in a dataframe. 7 dataframes are returned as output, 1 for each rank.  Nodes are represented by row. Each column represents a separate network centrality measure.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_degree_df_wo_unk(graph_wo_unk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph_wo_unk}{ List of 7 graphs, each only including known taxa for that taxonomic rank (from Kingdom to Species), produced by get_graph_wo_unk function.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns List of 7 dataframes of degree, betweenness, and closeness scores for all nodes in network.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{df1}{Dataframe of nodes with known Kingdom network measure scores}
\item{df2}{Dataframe of nodes with known Phylum network measure scores}
\item{df3}{Dataframe of nodes with known Class network measure scores}
\item{df4}{Dataframe of nodes with known Order network measure scores}
\item{df5}{Dataframe of nodes with known Family network measure scores}
\item{df6}{Dataframe of nodes with known Genus network measure scores}
\item{df7}{Dataframe of nodes with known Species network measure scores}
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
function(graph_wo_unk){
  met_degree_df_wo_unk_l = list()
  for(i in 1:length(graph_wo_unk)){
    graph_tax_level = names(graph_wo_unk)[[i]]
    print(graph_tax_level)
    graph = graph_wo_unk[[i]]
    graph_degree_df <- degree_calc_f(graph)
    met_degree_df_wo_unk_l[[graph_tax_level]] = graph_degree_df
  }
  return(met_degree_df_wo_unk_l)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
