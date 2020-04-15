\name{get_back_counts_for_lineplotf}
\alias{get_back_counts_for_lineplotf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find prevalence of known and MDM OTUs across a range of sample sizes
%%  ~~function to do ... ~~
}
\description{ Calculate prevalence of taxa across range of sample sizes, from 1 to 100.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_back_counts_for_lineplotf(prev_df, met_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prev_df}{ Dataframe of sample prevalence and taxonomic classification for all OTUs in OTU table
%%     ~~Describe \code{x} here~~
}
\item{met_name}{ Name of environment (title of plot)
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns lineplot of prevalence of known and MDM taxa.
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
function(prev_df, met_name) {
  newl <- lapply(list(1, 25, 50, 75, 100), function(val) {
  prev_df$Rank6 <- as.character(prev_df$Rank6)
    prev_df[grep("uncultured", prev_df$Rank6), ]$Rank6 <- "MDM" # convert all uncultured to MDM
    prev_df[grep("Ambiguous", prev_df$Rank6), ]$Rank6 <- "MDM" # convert all ambiguous taxa to MDM
    prev_df[is.na(prev_df$Rank6), ]$Rank6 <- "MDM" # convert all unassigned taxa to MDM
    prev_df[prev_df$Rank6 != "MDM", ]$Rank6 <- "Known" # convert rest of taxa (non-MDM) to Known
    total_val <- nrow(prev_df[prev_df$Prevalence > val, ])
    prev_df_above_val <- prev_df[prev_df$Prevalence > val, ] # subset to most prevalent taxa meeting val threshold
    num_MDM <- nrow(prev_df_above_val[prev_df_above_val$Rank6 == "MDM", ])
    print(num_MDM) # calculate number of MDM present in at least [val] samples
    num_K <- nrow(prev_df_above_val) - num_MDM
    print(num_K) # calculate number of known taxa present in at least [val] samples
    MDM_or_known_per_tax_level_df <- data.frame(MDM_type = c("MDM", "Known"), Num_OTU = c(num_MDM, num_K)) # create dataframe of # MDM/Known taxa present for each val threshold
    MDM_or_known_per_tax_level_df$Prev_val <- val
    return(MDM_or_known_per_tax_level_df)
  })
  newl_df <- do.call(rbind, newl) # bind together all prevalence threshold results in one dataframe
  newl_df$Met_Name <- met_name # create new column of environment - useful when comparing against multiple environments
  # return(newl_df)
  # create lineplot of prevalence of OTUs present in 1,25,50,75,100 samples
  g <- ggplot(newl_df, aes(Prev_val, Num_OTU, color = MDM_type)) +
    geom_line(aes(linetype = MDM_type)) + scale_y_log10() +
    theme_classic() + ylab("Number of OTUs (log scaled)") + xlab("Sample Prevalence") +
    ggtitle(paste(met_name, "OTU Prevalence", sep = " "))
  print(g)
  things_to_keep <- list(df = newl_df, prev_plot = g)
  return(things_to_keep)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
