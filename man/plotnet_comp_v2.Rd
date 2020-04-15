\name{plotnet_comp_v2}
\alias{plotnet_comp_v2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Modified Plot_network function to use specified graph layout
%%  ~~function to do ... ~~
}
\description{ To compare network attributes and shape between etworks created by different network methods, plot networks in same manner, using one graph as default layout for all networks. 2 graphs are used as input: graph of one network, graph of network used to compare between. Output is network in same layout as comparison network.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotnet_comp_v2(g, g_orig_for_layout, physeq = NULL, type = "samples", color = NULL, shape = NULL, point_size = 4, alpha = 1, label = "value", hjust = 1.35, line_weight = 0.5, line_color = color, line_alpha = 0.4, layout.method = layout.fruchterman.reingold, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{ Graph
%%     ~~Describe \code{g} here~~
}
  \item{g_orig_for_layout}{ Graph used to compare with.
%%     ~~Describe \code{g_orig_for_layout} here~~
}
  \item{physeq}{ Phyloseq
%%     ~~Describe \code{physeq} here~~
}
  \item{type}{ What to plot (choose between samples or taxa). Default: samples.
%%     ~~Describe \code{type} here~~
}
  \item{color}{ Color of nodes (choose between rank levels)
%%     ~~Describe \code{color} here~~
}
  \item{shape}{ Shape of nodes
%%     ~~Describe \code{shape} here~~
}
  \item{point_size}{ Size of nodes
%%     ~~Describe \code{point_size} here~~
}
  \item{alpha}{ Transparency of node color
%%     ~~Describe \code{alpha} here~~
}
  \item{label}{ Label of nodes
%%     ~~Describe \code{label} here~~
}
  \item{hjust}{ Adjust plot
%%     ~~Describe \code{hjust} here~~
}
  \item{line_weight}{ Weight of edges
%%     ~~Describe \code{line_weight} here~~
}
  \item{line_color}{ Color of edges
%%     ~~Describe \code{line_color} here~~
}
  \item{line_alpha}{ Transparency of edge color
%%     ~~Describe \code{line_alpha} here~~
}
  \item{layout.method}{ Layout method for network
%%     ~~Describe \code{layout.method} here~~
}
  \item{title}{ Title of network plot
%%     ~~Describe \code{title} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns network plot (ggplot object) in specified graph layout
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
function (g, g_orig_for_layout, physeq = NULL, type = "samples",
    color = NULL, shape = NULL, point_size = 4, alpha = 1, label = "value",
    hjust = 1.35, line_weight = 0.5, line_color = color, line_alpha = 0.4,
    layout.method = layout.fruchterman.reingold, title = NULL)
{
    if (vcount(g) < 2) {
        stop("The graph you provided, `g`, has too few vertices. \n         Check your graph, or the output of `make_network` and try again.")
    }
    if (type \%in\% c("taxa", "species", "OTUs", "otus", "otu")) {
        type <- "taxa"
    }
    edgeDF <- data.frame(get.edgelist(g))
    edgeDF$id <- 1:length(edgeDF[, 1])
    vertDF <- layout.method(g_orig_for_layout)
    colnames(vertDF) <- c("x", "y")
    vertDF <- data.frame(value = get.vertex.attribute(g, "name"),
        vertDF)
    if (!is.null(physeq)) {
        extraData <- NULL
        if (type == "samples" & !is.null(sample_data(physeq,
            FALSE))) {
            extraData = data.frame(sample_data(physeq))[as.character(vertDF$value),
                , drop = FALSE]
        }
        else if (type == "taxa" & !is.null(tax_table(physeq,
            FALSE))) {
            extraData = data.frame(tax_table(physeq))[as.character(vertDF$value),
                , drop = FALSE]
        }
        if (!is.null(extraData)) {
            vertDF <- data.frame(vertDF, extraData)
        }
    }
    graphDF <- merge(reshape2::melt(edgeDF, id = "id"), vertDF,
        by = "value")
    p <- ggplot(vertDF, aes(x, y))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank())
    p <- p + geom_point(aes_string(color = color, shape = shape),
        size = point_size, na.rm = TRUE)
    if (!is.null(label)) {
        p <- p + geom_text(aes_string(label = label), size = 2,
            hjust = hjust, na.rm = TRUE)
    }
    p <- p + geom_line(aes_string(group = "id", color = line_color),
        graphDF, size = line_weight, alpha = line_alpha, na.rm = TRUE)
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
