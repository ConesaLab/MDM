\name{orig_phylo_w_MDM_names}
\alias{orig_phylo_w_MDM_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create final phyloseq object with ambiguous, unassigned, and uncultured taxa reassigned to MDM
%%  ~~function to do ... ~~
}
\description{ Before beginning network analysis, we reassign any OTUs designated as ambiguous, unassigned or uncultured at each rank to MDM. Rank names are changed as well from Rank1 --> Kingdom, Rank2 --> Phylum, and so forth for Kingdom to Species.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
orig_phylo_w_MDM_names(phylo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{ Original OTU table in phyloseq format
%%     ~~Describe \code{phylo} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns
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
function (phylo, mapping)
{
    unassignedOTUs_MDM <- c("NA", "D_1__uncultured", "D_1__uncultured bacterium",
        "D_1__Unknown Phylum", "D_2__uncultured", "D_2__uncultured bacterium",
        "D_2__Unknown Class", "Ambiguous_taxa", "Unassigned",
        "D_3__uncultured", "D_3__uncultured bacterium", "D_3__Unknown Order",
        "D_4__uncultured", "D_4__uncultured bacterium", "D_4__Unknown Family",
        "D_5__uncultured", "D_5__uncultured bacterium", "D_5__Unknown Genus",
        "D_6__uncultured", "D_6__uncultured bacterium", "D_6__Unknown Species")
    silva_otus = phyloseq::otu_table(phylo)
    silva_taxa = phyloseq::tax_table(phylo)
    silva_taxa_df <- as.data.frame(silva_taxa)
    colnames(silva_taxa_df) = c("Kingdom", "Phylum", "Class",
        "Order", "Family", "Genus", "Species")
    silva_taxa_df$Phylum[silva_taxa_df$Phylum \%in\% unassignedOTUs_MDM] <- "<NA>"
    silva_taxa_df$Class[silva_taxa_df$Class \%in\% unassignedOTUs_MDM] <- "<NA>"
    silva_taxa_df$Order[silva_taxa_df$Order \%in\% unassignedOTUs_MDM] <- "<NA>"
    silva_taxa_df$Family[silva_taxa_df$Family \%in\% unassignedOTUs_MDM] <- "<NA>"
    silva_taxa_df$Genus[silva_taxa_df$Genus \%in\% unassignedOTUs_MDM] <- "<NA>"
    silva_taxa_df$Species[silva_taxa_df$Species \%in\% unassignedOTUs_MDM] <- "<NA>"
    silva_taxa_df[] <- lapply(silva_taxa_df, as.character)
    silva_taxa_df[] <- lapply(silva_taxa_df, function(x) {
        stringr::str_replace_all(x, "uncultured", NA_character_)
    })
    newsilvataxtable2_m <- as.matrix(silva_taxa_df)
    newsilva_taxa <- tax_table(newsilvataxtable2_m)
    new_silva_phylo = phyloseq(silva_otus, newsilva_taxa, mapping)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
