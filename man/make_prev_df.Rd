\name{make_prev_df}
\alias{make_prev_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make dataframe of taxa prevalence, abundance, and classification
%%  ~~function to do ... ~~
}
\description{ Using the resultant OTU table biom file produced by QIIME, find the relative abundance and sample prevalence of all OTUs and store results in dataframe. Resultant dataframe will show taxonomic classification, abundance, and prevalence for all OTUs.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_prev_df(biom_file)
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
\value{ Returns dataframe combining sample prevalence, abundance, and taxonomic classification for all OTUs.
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
  phylo <- import_biom(biom_file)
  prevdf = apply(X = otu_table(phylo),MARGIN = ifelse(taxa_are_rows(ps), yes=1, no=2), FUN=function(x){sum(x>0)}) #find prevalence of each taxa per sample
  #bind together taxonomic information and prevalence and total abundance of each taxa into one dataframe
  prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(phylo), tax_table(phylo))
  prevdf <- prevdf[prevdf$Prevalence > 0,] #remove taxa not present in samples
  return(prevdf)
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
