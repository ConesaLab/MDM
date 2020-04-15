\name{get_bar_and_dens_class_interactions}
\alias{get_bar_and_dens_class_interactions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Do node neighbor comparison for all taxa within network
%%  ~~function to do ... ~~
}
\description{ Find Class neighbors of nodes for all nodes in network and calculate frequency of these interactions (frequency of shared edges). Create bar chart of the frequency of each neighbor interaction. Calculate self-self interactions and create density plots to compare all self-self/all possible interactions vs MDM self-self interactions
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_bar_and_dens_class_interactions(orig_graph, orig_phylo, met_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_graph}{ Graph object of all taxa present in network
%%     ~~Describe \code{orig_graph} here~~
}
  \item{orig_phylo}{ Phyloseq object of all taxa in network
%%     ~~Describe \code{orig_phylo} here~~
}
  \item{met_type}{ Name of environment
%%     ~~Describe \code{met_type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns List of 4 (list) objects
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{barplot_of_interactions_l}{For each unique class node, the frequency of interactions with all neighbors are shown in bar chart format}
\item{pie_b_plot_interactions_l}{For each unique class node, the frequency of interactions with all neighbors are shown in pie chart format}
\item{dens_plot_all_poss2}{Density plots showing all possible class interactions vs MDM self-self interactions}
\item{dens_plot_just_sameclass_poss2}{Density plots showing all self-self class interactions vs MDM self-self interactions}

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
function (orig_graph, orig_phylo, met_type)
{
    edgeDF = as.data.frame(igraph::as_edgelist(orig_graph))
    colnames(edgeDF) = c("X1", "X2")
    edgeDF_w_alltaxnames <- (merge(tax_table(orig_phylo), edgeDF,
        by.x = "row.names", by.y = "X1"))
    edgeDF_w_alltaxnames_updated <- merge(edgeDF_w_alltaxnames,
        tax_table(orig_phylo), by.x = "X2", by.y = "row.names",
        all.x = TRUE)
    print(nrow(edgeDF_w_alltaxnames_updated))
    identical(nrow(edgeDF_w_alltaxnames_updated), nrow(edgeDF))
    colnames(edgeDF_w_alltaxnames_updated)[1:2] = c("Y", "X")
    edgeDF_w_alltaxnames_updated$Class.x <- as.character(edgeDF_w_alltaxnames_updated$Class.x)
    edgeDF_w_alltaxnames_updated$Class.y <- as.character(edgeDF_w_alltaxnames_updated$Class.y)
    edgeDF_w_alltaxnames_updated[is.na(edgeDF_w_alltaxnames_updated$Class.x),
        ]$Class.x = "Unknown"
    edgeDF_w_alltaxnames_updated[is.na(edgeDF_w_alltaxnames_updated$Class.y),
        ]$Class.y = "Unknown"
    edgeDF_w_alltaxnames_updated[edgeDF_w_alltaxnames_updated$Y ==
        "0", ]$Class.y = "Unknown"
    edgeDF_w_alltaxnames_updated$combinedXY = with(edgeDF_w_alltaxnames_updated,
        paste0(Class.x, "&", Class.y))
    barplot_of_interactions_l <- list()
    pie_b_plot_interactions_l <- list()
    for (x_class in unique(edgeDF_w_alltaxnames_updated$Class.x)) {
        edge_df_for_class = as.data.frame(table(edgeDF_w_alltaxnames_updated[grep(x_class,
            edgeDF_w_alltaxnames_updated$combinedXY), ]$combinedXY))
        x_class_at_end = paste0("&", x_class)
        all_to_rep <- edge_df_for_class$Var1[grep(x_class_at_end,
            edge_df_for_class$Var1)]
        test_st <- strsplit(as.character(all_to_rep), "&")
        test_st2 <- lapply(test_st, function(x) {
            new_var <- paste0(x_class, "&", x[1])
            return(new_var)
        })
        edge_df_for_class$Var1 <- as.character(edge_df_for_class$Var1)
        edge_df_for_class$Var1[all_to_rep] = (unlist(test_st2))
        edge_df_for_class <- edge_df_for_class \%>\% mutate(size = Freq/sum(Freq))
        same_class_var = paste0(x_class, "&", x_class)
        edge_df_for_class$Var1 <- as.factor(edge_df_for_class$Var1)
        b_plot <- ggplot(edge_df_for_class, aes(Var1, Freq, fill = Var1)) +
            geom_bar(stat = "identity")
        b_plot <- b_plot + geom_bar(stat = "identity", data = edge_df_for_class[edge_df_for_class$Var1 \%in\%
            same_class_var, ], aes(Var1), size = 1, fill = "grey50",
            color = "black") + coord_flip() + theme(legend.position = "none",
            axis.title.y = element_blank())
        b_plot <- b_plot + ggtitle(paste("Neighbor Interactions for Class ",
            x_class, sep = " "))
        custom_colors = rainbow(length(unique(edge_df_for_class$Var1)))
        num_to_change_to_grey = which(levels(edge_df_for_class$Var1) ==
            same_class_var)
        custom_colors[num_to_change_to_grey] = "grey50"
        pie_b_plot <- ggplot(edge_df_for_class, aes(x = "", y = size,
            fill = Var1)) + geom_bar(width = 1, stat = "identity",
            position = position_fill()) + scale_fill_manual(values = custom_colors) +
            coord_polar("y", start = 0) + geom_text(aes(label = paste0(round(size *
            100), "\%")), position = position_fill(vjust = 0.6),
            size = 3) + theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(), panel.border = element_blank(),
            panel.grid = element_blank(), panel.background = element_blank(),
            axis.ticks = element_blank(), axis.text.x = element_blank()) +
            theme(legend.text = element_text(size = 7)) + guides(color = guide_legend(override.aes = list(size = 1)))
        barplot_of_interactions_l[[x_class]] = b_plot
        pie_b_plot_interactions_l[[x_class]] = pie_b_plot
    }
    combinedXY_df = as.data.frame(table(edgeDF_w_alltaxnames_updated$combinedXY))
    val_of_unk1 = combinedXY_df[which(combinedXY_df$Var1 == "Unknown&Unknown"),
        ]$Freq
    next_smallest_val1 = val_of_unk1 - 1
    p_val_sig_dif = length(which(combinedXY_df$Freq > next_smallest_val1))/length(combinedXY_df$Freq)
    p_val_sig_dif = round(p_val_sig_dif, digits = 5)
    d <- density(combinedXY_df$Freq)
    y_axis_max_dens <- max(d$y)
    dens_plot_all_poss <- ggplot(combinedXY_df, aes(Freq)) +
        geom_density() + geom_vline(xintercept = val_of_unk1,
        linetype = "dashed", color = "red", size = 1) + ggtitle(paste(met_type,
        "All Possible Interactions", sep = " "))
    dens_plot_all_poss2 <- dens_plot_all_poss + geom_label(aes(x = val_of_unk1,
        y = y_axis_max_dens, label = paste("P-value =", p_val_sig_dif,
            sep = " ")), color = "black", angle = 0, hjust = "inward")
    dens_plot_all_poss2 <- dens_plot_all_poss2 + geom_label(aes(x = val_of_unk1,
        y = y_axis_max_dens/2, label = paste("Unknown-Unknown"),
        hjust = "inward"))
    print(dens_plot_all_poss2)
    combinedXY_just_sameclass_df = as.data.frame(table(edgeDF_w_alltaxnames_updated[edgeDF_w_alltaxnames_updated$Class.x ==
        edgeDF_w_alltaxnames_updated$Class.y, ]$combinedXY))
    val_of_unk2 = combinedXY_just_sameclass_df[which(combinedXY_just_sameclass_df$Var1 ==
        "Unknown&Unknown"), ]$Freq
    next_smallest_val2 = val_of_unk2 - 1
    p_val_sig_dif2 = length(which(combinedXY_just_sameclass_df$Freq >
        next_smallest_val2))/length(combinedXY_just_sameclass_df$Freq)
    p_val_sig_dif2 = round(p_val_sig_dif2, digits = 5)
    d2 <- density(combinedXY_just_sameclass_df$Freq)
    y_axis_max_dens2 <- max(d2$y)
    dens_plot_just_sameclass_poss <- ggplot(combinedXY_just_sameclass_df,
        aes(Freq)) + geom_density() + geom_vline(xintercept = val_of_unk2,
        linetype = "dashed", color = "red", size = 1) + ggtitle(paste(met_type,
        "Self-self Class Interactions", sep = " "))
    dens_plot_just_sameclass_poss2 <- dens_plot_just_sameclass_poss +
        geom_label(aes(x = val_of_unk2, y = y_axis_max_dens2,
            label = paste("P-value =", p_val_sig_dif2, sep = " ")),
            color = "black", angle = 0, hjust = "inward")
    dens_plot_just_sameclass_poss2 <- dens_plot_just_sameclass_poss2 +
        geom_label(aes(x = val_of_unk2, y = y_axis_max_dens2/2,
            label = paste("Unknown-Unknown"), hjust = "inward"))
    print(dens_plot_just_sameclass_poss2)
    bar_pie_and_dens_plots_l <- list(barplot_of_interactions_l,
        pie_b_plot_interactions_l, dens_plot_all_poss2, dens_plot_just_sameclass_poss2)
    return(bar_pie_and_dens_plots_l)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
