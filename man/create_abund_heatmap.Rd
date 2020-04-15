\name{create_abund_heatmap}
\alias{create_abund_heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create heatmap of OTU abundance
%%  ~~function to do ... ~~
}
\description{ Creates heatmap of taxa prevalence across samples for a given environment. OTU Counts are clr-normalized and shown using a blue-yellow color-mapping scale. Color legends are provided for the Class rank, category (MDM or Known), country (geographic location), and project.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_abund_heatmap(orig_phylo_w_map, met_type_title, filepathname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_phylo_w_map}{ Phyloseq object including sample and taxonomic information
%%     ~~Describe \code{orig_phylo_w_map} here~~
}
  \item{met_type_title}{ Name of title for heatmap
%%     ~~Describe \code{met_type_title} here~~
}
  \item{filepathname}{ Name of file path where heatmap image will be produced
%%     ~~Describe \code{met_type_title} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns pdf file of the OTU abundance heatmap
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
function (orig_phylo_w_map, met_type_title, filepathname)
{
    manualcolors <- c("black", "forestgreen", "red2", "orange",
        "cornflowerblue", "magenta", "darkolivegreen4", "indianred1",
        "tan4", "darkblue", "mediumorchid1", "firebrick4", "yellowgreen",
        "lightsalmon", "tan3", "tan1", "darkgray", "wheat4",
        "#DDAD4B", "chartreuse", "seagreen1", "moccasin", "mediumvioletred",
        "seagreen", "cadetblue1", "darkolivegreen1", "tan2",
        "tomato3", "#7CE3D8", "gainsboro")
    met_otu_table_mat <- as(otu_table(orig_phylo_w_map)[-nrow(otu_table(orig_phylo_w_map)),
        ], "matrix")
    met_otu_table_mat_clr <- SpiecEasi::clr(met_otu_table_mat)
    dendo_otunames2 <- hclust(dist(met_otu_table_mat_clr))
    dendo_samplenames2 <- hclust(dist(t(met_otu_table_mat_clr)))
    otu_label2 <- as.data.frame(dendo_otunames2$labels)
    met_tax <- as.data.frame(tax_table(orig_phylo_w_map))
    met_tax$ID <- rownames(met_tax)
    otu_label2_df2 <- merge(otu_label2, met_tax, by.x = "dendo_otunames2$labels",
        by.y = "ID", all.x = TRUE)
    otu_label2_df2[, ] <- lapply(otu_label2_df2[, ], as.character)
    otu_label2_df2[is.na(otu_label2_df2)] <- "Unknown"
    new_otu_colors_clr <- with(otu_label2_df2, data.frame(otu_Class_label = levels(factor(Class)),
        color = rainbow(length(unique(Class)))))
    otu_label2_df2_clr <- merge(otu_label2_df2, new_otu_colors_clr,
        by.x = "Class", by.y = "otu_Class_label", all.x = TRUE)
    colnames(otu_label2_df2_clr)[9] = "otu_colors"
    otu_label2_df2_clr$otu_colors <- as.character(otu_label2_df2_clr$otu_colors)
    otu_label2_df2_clr[otu_label2_df2_clr$Class == "Unknown",
        ]$otu_colors <- "black"
    otu_label2_df2_clr$org_type <- "white"
    otu_label2_df2_clr[otu_label2_df2_clr$Class == "Unknown",
        ]$org_type <- "black"
    new_otu_sample_col_clr <- as.matrix(rbind(otu_label2_df2_clr$otu_colors,
        otu_label2_df2_clr$org_type))
    rownames(new_otu_sample_col_clr) <- c("Class", "Category")
    sample_label_clr <- as.data.frame(dendo_samplenames2$labels)
    sample_data_from_map <- as.data.frame(as.matrix(sample_data(orig_phylo_w_map)))
    sample_label_df2_clr <- merge(sample_label_clr, sample_data_from_map,
        by.x = "dendo_samplenames2$labels", by.y = "SampleID",
        all.x = TRUE)
    num_of_proj <- length(unique(sample_label_df2_clr$ProjectName))
    if (num_of_proj > 11) {
        new_proj_colors_clr <- new_proj_colors_clr <- with(sample_label_df2_clr,
            data.frame(Project = levels(factor(ProjectName)),
                color = sample(manualcolors, num_of_proj, replace = FALSE)))
    }
    else {
        new_proj_colors_clr <- with(sample_label_df2_clr, data.frame(Project = levels(factor(ProjectName)),
            color = brewer.pal(num_of_proj, "Spectral")))
    }
    new_sample_label_df2_clr <- merge(sample_label_df2_clr, new_proj_colors_clr,
        by.x = "ProjectName", by.y = "Project", all.x = TRUE)
    colnames(new_sample_label_df2_clr)[16] = "colors_of_samples"
    new_sample_label_df2_clr$Location <- as.character(new_sample_label_df2_clr$Location)
    hscountrynames <- gsub("\\:.*", "", new_sample_label_df2_clr$Location)
    hscountrynames2 <- gsub("\\_.*", "", hscountrynames)
    new_sample_label_df2_clr$Country = hscountrynames2
    num_to_change <- length(unique(new_sample_label_df2_clr$Country))
    new_country_colors_clr <- with(new_sample_label_df2_clr,
        data.frame(Country = levels(factor(Country)), color = brewer.pal(num_to_change,
            "Paired")))
    new_sample_label_df2_clr2 <- merge(new_sample_label_df2_clr,
        new_country_colors_clr, by.x = "Country", by.y = "Country",
        all.x = TRUE)
    sample_label_colors2_clr <- as.matrix(cbind(as.character(new_sample_label_df2_clr2$colors_of_samples),
        as.character(new_sample_label_df2_clr2$color)))
    colnames(sample_label_colors2_clr) = c("Project", "Country")
    par(mar = c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
    pdf(paste0(filepathname, met_type_title, "_", "clr_transformed", ".pdf"))
    heatmap.3(x = met_otu_table_mat_clr, col = viridis::viridis, scale = "none",
      ColSideColors = sample_label_colors2_clr, ColSideColorsSize = 2,
      RowSideColors = new_otu_sample_col_clr, RowSideColorsSize = 2,
      keysize = 1, KeyValueName = "OTU Count (after clr)", labCol = FALSE,
      labRow = FALSE, xlab = "Samples", ylab = "OTUs", main = met_type_title)
    legend("topright", legend = c(unique(new_sample_label_df2_clr2$Country)),
        fill = c(unique(as.character(new_sample_label_df2_clr2$color))),
        bty = "n", y.intersp = 0.7, cex = 0.7, ncol = 2, yjust = 0,
        xjust = 1.5, xpd = TRUE, inset = c(0, -0.02))
    legend("bottomleft", legend = c("Known", "Unknown"), fill = c("White",
        "Black"), border = TRUE, bty = "n", y.intersp = 0.7,
        cex = 0.7, xjust = 0)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
