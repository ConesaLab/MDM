\name{get_back_res_meeting_min_occ}
\alias{get_back_res_meeting_min_occ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create phyloseq and igraph objects of taxa meeting certain sample percent threshold
%%  ~~function to do ... ~~
}
\description{ To create suitable networks for each environment, we model the most prevalent taxa, excluding rarer taxa that do not meet a strict sample prevalence threshold. Using the original OTU table information, we set a threshold of sample prevalence percentage and discard any taxa not meeting this threshold. All taxa that remain are saved to a new phyloseq object containing the updated, reduced OTU table information. The direct co-occurrence relationships are then measured for remaining taxa using the SpiecEasi MB neighborhood algorithm and stored in an igraph object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_back_res_meeting_min_occ(phylo, filter_val_percent = 0.4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{ Phyloseq object of OTU table and taxonomic information
%%     ~~Describe \code{phylo} here~~
}
  \item{filter_val_percent}{ Minimum sample percent threshold (from 0 to 1) taxa must meet to be included in the environmental network. Default threshold  is 0.4 (40 percent of all samples). Using this default threshold would mean that all taxa present in less than 40 percent of all samples will be discarded from network creation and subsequent steps.
%%     ~~Describe \code{filter_val_percent} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns list of 2 objects.
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{silva_phyloseqobj.final }{Phyloseq object including only taxa meeting sample prevalence threshold}
  \item{silva_spiec.graph }{Igraph object of direct co-occurrence relationships between these taxa calculated by SpiecEasi MB algorithm.}
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
function (phylo, filter_val_percent = 0.4)
{
    silva_otus = phyloseq::otu_table(phylo)
    newsilva_taxa = phyloseq::tax_table(phylo)
    silva_filterobj = seqtime::filterTaxonMatrix(silva_otus,
        minocc = filter_val_percent * length(rownames(sample_data(phylo))),
        keepSum = TRUE, return.filtered.indices = TRUE)
    silva_otus.f = silva_filterobj$mat
    silva_taxa.f = newsilva_taxa[setdiff(1:nrow(newsilva_taxa),
        silva_filterobj$filtered.indices), ]
    dummyTaxonomy = c("D_0__dummy", "D_1__", "D_2__", "D_3__",
        "D_4__", "D_5__", "D_6__")
    silva_taxa.f = rbind(silva_taxa.f, dummyTaxonomy)
    rownames(silva_taxa.f)[nrow(silva_taxa.f)] = "0"
    rownames(silva_otus.f)[nrow(silva_otus.f)] = "0"
    silva_updatedotus = otu_table(silva_otus.f, taxa_are_rows = TRUE)
    silva_updatedtaxa = tax_table(silva_taxa.f)
    silva_phyloseqobj.final = phyloseq::phyloseq(silva_updatedotus,
        silva_updatedtaxa)
    silva_spiec.out = SpiecEasi::spiec.easi(silva_phyloseqobj.final,
        method = "mb", icov.select.params = list(rep.num = 20))
    silva_spiec.graph = igraph::adj2igraph(silva_spiec.out$refit,
        vertex.attr = list(name = taxa_names(silva_phyloseqobj.final)))
    phylo_and_graph <- list(silva_phyloseqobj.final, silva_spiec.graph)
    return(phylo_and_graph)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
