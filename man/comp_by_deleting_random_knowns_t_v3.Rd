\name{comp_by_deleting_random_knowns_t_v3}
\alias{comp_by_deleting_random_knowns_t_v3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create bootstrap networks and compare network measures for Original, Without Unknown, and Bootstrap networks for all taxonomic rank level
%%  ~~function to do ... ~~
}
\description{ Bootstrap networks are generated, removing an equal number of random knowns as MDM at each rank and all network measures are recalculated for nodes in these networks.To compare network measures for Original, Without Unknown, and Bootstrap networks for all taxonomic rank level, a dataframe is created, showing the node values for degree, betweenness, and closeness for nodes present in each network type at each taxonomic rank.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comp_by_deleting_random_knowns_t_v3(orig_graph, new_wo_unk_graph,
        orig_phylo, orig_df, degree_df_wo_unk, iter = 100, mc.coreval=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_graph}{ Network including all taxa (Original network), produced from get_back_res_meeting_min_occ function
%%     ~~Describe \code{orig_graph} here~~
}
  \item{new_wo_unk_graph}{ Network excluding MDM taxa (Without Unknown network), produced by get_graph_wo_unk function
%%     ~~Describe \code{new_wo_unk_graph} here~~
}
  \item{orig_phylo}{Phyloseq including all taxa meeting threshold, produced from get_back_res_meeting_min_occ function
%%     ~~Describe \code{orig_phylo} here~~
}
  \item{orig_df}{ Dataframe of degree, betweenness, and closeness scores for all nodes in Original network, produced by degree_calc_df function
%%     ~~Describe \code{orig_df} here~~
}
  \item{degree_df_wo_unk}{ Dataframe of degree, betweenness, and closeness scores for all nodes in Without Unknown Network, produced by get_degree_df_wo_unk function
%%     ~~Describe \code{degree_df_wo_unk} here~~
}
  \item{iter}{ Number of iterations of random subsampling and removal of knowns. Default is 100 iterations.
%%     ~~Describe \code{iter} here~~
}
\item{mc.coreval}{ Number of cores to use for this function to parallelize this long process. Default is 2 to allow this function to run on a computer.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns dataframe of 4 columns, respectively indicating taxonomic rank (Kingdom to Species), network type (Original, Bootstrap, or Without_unk), measure (degree, betweenness, closeness), and score for each node.
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
function (orig_graph, new_wo_unk_graph, orig_phylo, orig_df,
    degree_df_wo_unk, iter = 100, mc.coreval=2)
{
    orig_vertices <- as.vector(names(V(orig_graph)))
    measurements <- c("degree", "bw", "closeness")
    all_net_tax_level_plot_l <- do.call(rbind, parallel::mclapply(seq_along(new_wo_unk_graph),
        mc.cores = mc.coreval, function(i) {
            rank_level = names(new_wo_unk_graph)[[i]]
            print(paste0("Rank level: ", rank_level))
            num_nodes_to_remove <- length(V(orig_graph)) - length(V(new_wo_unk_graph[[rank_level]]))
            unk_names <- setdiff(names(V(orig_graph)), names(V(new_wo_unk_graph[[rank_level]])))
            orig_vertices_wo_unk <- intersect(names(V(orig_graph)),
                names(V(new_wo_unk_graph[[rank_level]])))
            bstrap_df_l <- lapply(seq_len(iter), function(i) {
                vecnames <- sample(x = orig_vertices_wo_unk,
                  size = num_nodes_to_remove, replace = TRUE)
                wo_unk_graph <- igraph::delete_vertices(orig_graph,
                  which(names(V(orig_graph)) \%in\% vecnames))
                boot_degree_df <- degree_calc_f(wo_unk_graph)
                if (!identical(length(rownames(boot_degree_df)),
                  length(rownames(degree_df_wo_unk[[rank_level]])))) {
                  boot_degree_df = boot_degree_df[seq_along(rownames(degree_df_wo_unk[[rank_level]])),
                    ]
                }
                boot_degree_df[, measurements]
            })
            do.call(rbind, lapply(measurements, function(act_meas) {
                orig_meas <- orig_df[[act_meas]]
                wo_unk_meas <- degree_df_wo_unk[[rank_level]][[act_meas]]
                bstrap_meas <- rowMeans(do.call(cbind, lapply(bstrap_df_l,
                  function(x) x[, act_meas])))
                data.frame(rank_level = rank_level, type = factor(c(rep("Original",
                  length(orig_meas)), rep("Without_unknown",
                  length(wo_unk_meas)), rep("Bootstrap", length(bstrap_meas))),
                  levels = c("Original", "Without_unknown", "Bootstrap")),
                  measure = act_meas, data = c(orig_meas, wo_unk_meas,
                    bstrap_meas))
            }))
        }))
    return(all_net_tax_level_plot_l)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
