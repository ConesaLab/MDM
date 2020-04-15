\name{make_mdt}
\alias{make_mdt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make melted dataframe of relative abundance of each taxa
%%  ~~function to do ... ~~
}
\description{ Using as input the OTU table biom produced by QIIME analysis, make melted dataframe of relative abundance of each taxa.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_mdt(biom_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biom_file}{ Full path name of OTU table biom file
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
function(biom_file){
  orig_phylo <- import_biom(biom_file)
  otutab = as(otu_table(orig_phylo), "matrix")
    otudt = data.table(otutab, keep.rownames = TRUE)
    setnames(otudt, "rn", "TaxaID")
    otudt[, TaxaIDchar := as.character(TaxaID)]
    otudt[, TaxaID := NULL]
    setnames(otudt, "TaxaIDchar", "TaxaID")
    # Melt count table
    mdt = melt.data.table(otudt,
                          id.vars = "TaxaID",
                          variable.name = "SampleID",
                          value.name = "count")

    mdt <- mdt[count > 0]
    # Omit NAs
    mdt <- mdt[!is.na(count)]
    mdt[, RelativeAbundance := count / sum(count), by = SampleID] #find relative abundance by sample IDs
    taxdt = data.table(as(tax_table(orig_phylo, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt] #dataframe, 11 columns = last 2 columns = count, RelativeAbundance, columns 1-7 = taxonomy
    return(mdt)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
