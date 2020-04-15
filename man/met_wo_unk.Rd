\name{met_wo_unk}
\alias{met_wo_unk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ For each rank from Phylum to Genus, remove all MDM
%%  ~~function to do ... ~~
}
\description{ For each rank from Phylum to Genus, remove all MDM from the original OTU table and retain only taxa with known classifications. A reduced OTU table is produced for each rank, retaining only known taxa.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
met_wo_unk(orig_phylo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_phylo}{ Phyloseq including all taxa meeting threshold, produced from get_back_res_meeting_min_occ function
%%     ~~Describe \code{orig_phylo} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns list of 7 OTU tables
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{table1}{OTU table of all known Kingdom taxa}
\item{table2}{OTU table of all known Phylum taxa}
\item{table3}{OTU table of all known Class taxa}
\item{table4}{OTU table of all known Order taxa}
\item{table5 }{OTU table of all known Family taxa}
\item{table6}{OTU table of all known Genus taxa}
\item{table7 }{OTU table of all known Species taxa}

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
function (orig_phylo)
{
    nodes_list_m = list()
    met_tax_table = phyloseq::tax_table(orig_phylo)
    print(dim(met_tax_table))
    met_tax_table <- as.data.frame(met_tax_table)
    met_tax_table[] <- lapply(met_tax_table, as.character)
    met_tax_table1 <- met_tax_table \%>\% dplyr::mutate_all(funs(replace(.,
        is.na(.), "Unassigned")))
    rownames(met_tax_table1) = rownames(met_tax_table)
    print(dim(met_tax_table1))
    for (i in 1:length(colnames(met_tax_table1))) {
        print(i)
        unassignedOTUs_MDM <- c("NA", "D_1__uncultured", "D_1__uncultured bacterium",
            "D_1__Unknown Phylum", "D_2__uncultured", "D_2__uncultured bacterium",
            "D_2__Unknown Class", "Ambiguous_taxa", "Unassigned",
            "D_3__uncultured", "D_3__uncultured bacterium", "D_3__Unknown Order",
            "D_4__uncultured", "D_4__uncultured bacterium", "D_4__Unknown Family",
            "D_5__uncultured", "D_5__uncultured bacterium", "D_5__Unknown Genus",
            "D_6__uncultured", "D_6__uncultured bacterium", "D_6__Unknown Species")
        node_a <- met_tax_table1[met_tax_table1[, i] \%in\% unassignedOTUs_MDM,
            ]
        print(dim(node_a))
        nodes_list_m[[i]] = node_a
    }
    final_otus = phyloseq::otu_table(orig_phylo)
    print(dim(final_otus))
    wo_unk_otus_at_tax_level_l = list()
    for (t in 1:length(nodes_list_m)) {
        print(t)
        unk_otus_at_tax_level <- final_otus[!rownames(final_otus) \%in\%
            rownames(nodes_list_m[[t]]), ]
        print(dim(unk_otus_at_tax_level))
        wo_unk_otus_at_tax_level_l[[t]] = unk_otus_at_tax_level
    }
    names(wo_unk_otus_at_tax_level_l) = colnames(tax_table(orig_phylo))
    return(wo_unk_otus_at_tax_level_l)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
