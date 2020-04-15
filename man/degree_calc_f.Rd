\name{degree_calc_f}
\alias{degree_calc_f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate and store degree, betweenness, and closeness scores in dataframe
%%  ~~function to do ... ~~
}
\description{ Calculate degree, betweenness, and closeness scores for each node in network (igraph) and store results in dataframe. Nodes are represented by row. Each column represents a separate network centrality measure.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
degree_calc_f(orig_graph)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_graph}{ Graph object of taxa co-occurrence relationships
%%     ~~Describe \code{orig_graph} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns dataframe of degree, betweenness, and closeness scores for all nodes in network
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
function (orig_graph)
{
    degree_df <- as.data.frame(igraph::degree(orig_graph))
    degree_df$names <- rownames(degree_df)
    colnames(degree_df)[1] = "degree"
    degree_df$bw <- as.vector(igraph::betweenness(orig_graph))
    degree_df$closeness <- as.vector(igraph::closeness(orig_graph))
    print(head(degree_df))
    return(degree_df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
