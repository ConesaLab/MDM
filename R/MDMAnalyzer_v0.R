#####################################################
#### mdmnets package                          ###
######################################################

# #Load required libraries
# library(data.table)
# library(reshape2)
# library(phyloseq)
# library(dplyr)
# library(SpiecEasi)
# library(seqtime)
# library(ggplot2)
# library(igraph)
# library(gridExtra)
# library(ggpubr)
# library(gplots)
# library(devtools)
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
# # library(RColorBrewer)
# library(Hmisc)
# library(Matrix)
#
# load("~/NetsforShiny.RData")
# load("~/Updated_phyloseqobj_final.RData")

### to create package:
# library(devtools)
# library(roxygen2)
# devtools::create("mdmnets")

# load biom/mapping data into phyloseq format in R
# library(phyloseq)
# met_biom <- import_biom("otu_table_mc2_w_tax.biom")
# met_sample_data <- import_qiime_sample_data("Sample_Map.txt")
#mapping_MDM <- import_qiime_sample_data("~/FDMap_new1.txt")
# met_phylo <- merge_phyloseq(met_biom, met_sample_data)
# x= met_phylo

# unassignedOTUs_MDM = all uncultured/unknown, ambiguous,NA terms
unassignedOTUs_MDM <- c(
  "NA", "D_1__uncultured", "D_1__uncultured bacterium", "D_1__Unknown Phylum",
  "D_2__uncultured", "D_2__uncultured bacterium", "D_2__Unknown Class",
  "Ambiguous_taxa", "Unassigned",
  "D_3__uncultured", "D_3__uncultured bacterium", "D_3__Unknown Order",
  "D_4__uncultured", "D_4__uncultured bacterium",
  "D_4__Unknown Family", "D_5__uncultured", "D_5__uncultured bacterium",
  "D_5__Unknown Genus", "D_6__uncultured",
  "D_6__uncultured bacterium",
  "D_6__Unknown Species"
)
#

#' @import igraph
#' @importFrom ggpubr stat_compare_means
#' @importFrom dplyr group_by summarise "%>%" mutate_all funs mutate
#' @importFrom SpiecEasi spiec.easi adj2igraph clr sparcc
#' @importFrom stringr str_replace_all
#' @importFrom stats var
#' @importFrom Matrix Matrix
#' @importFrom Hmisc rcorr
#' @importFrom parallel mclapply
#' @import RColorBrewer
#' @importFrom reshape2 melt
#' @importFrom viridis viridis
#' @import phyloseq
#' @import ggplot2
#' @export
get_back_res_meeting_min_occ <- function(phylo, filter_val_percent = 0.4) {
  # phylo = phyloseq object, filter_val_percent = percentage of samples taxa must be present
  silva_otus <- otu_table(phylo)
  newsilva_taxa <- tax_table(phylo)
  silva_filterobj <- filterTaxonMatrix(silva_otus, minocc = filter_val_percent * length(rownames(sample_data(phylo))), keepSum = TRUE, return.filtered.indices = TRUE)
  silva_otus.f <- silva_filterobj$mat
  # print(dim(silva_otus.f))
  silva_taxa.f <- newsilva_taxa[setdiff(1:nrow(newsilva_taxa), silva_filterobj$filtered.indices), ]
  # print(dim(silva_taxa.f))
  dummyTaxonomy <- c("D_0__dummy", "D_1__", "D_2__", "D_3__", "D_4__", "D_5__", "D_6__")
  silva_taxa.f <- rbind(silva_taxa.f, dummyTaxonomy)
  rownames(silva_taxa.f)[nrow(silva_taxa.f)] <- "0"
  rownames(silva_otus.f)[nrow(silva_otus.f)] <- "0"
  silva_updatedotus <- otu_table(silva_otus.f, taxa_are_rows = TRUE)
  silva_updatedtaxa <- tax_table(silva_taxa.f)

  silva_phyloseqobj.final <- phyloseq(silva_updatedotus, silva_updatedtaxa)
  silva_spiec.out <- spiec.easi(silva_phyloseqobj.final, method = "mb", icov.select.params = list(rep.num = 20))
  silva_spiec.graph <- adj2igraph(silva_spiec.out$refit, vertex.attr = list(name = taxa_names(silva_phyloseqobj.final)))
  phylo_and_graph <- list(silva_phyloseqobj.final, silva_spiec.graph)
  return(phylo_and_graph)
  # silva_spiec_plot <- plot_network(silva_spiec.graph, silva_phyloseqobj.final, type='taxa', color="Rank3", label=NULL)
  # return(silva_spiec_plot)
}
#
tax_rank_names <- c("Phylum", "Class", "Order", "Family", "Genus")
#
# ### change MDM-related taxa to NA and change rank names to taxonomic classification names
# #after checking proportions of uncultured, unassigned, ambiguous, change these to NA accordingly
# #change to NA so that in igraph/net plots, all MDM nodes will be gray
# #to get back igraph to create network for each Met (HS, ARC = minocc = 0.4, HY, DS = minocc = 0.3)
#' @export
orig_phylo_w_MDM_names <- function(phylo) {
  unassignedOTUs_MDM <- c(
    "NA", "D_1__uncultured", "D_1__uncultured bacterium", "D_1__Unknown Phylum",
    "D_2__uncultured", "D_2__uncultured bacterium", "D_2__Unknown Class",
    "Ambiguous_taxa", "Unassigned",
    "D_3__uncultured", "D_3__uncultured bacterium", "D_3__Unknown Order",
    "D_4__uncultured", "D_4__uncultured bacterium",
    "D_4__Unknown Family", "D_5__uncultured", "D_5__uncultured bacterium",
    "D_5__Unknown Genus", "D_6__uncultured",
    "D_6__uncultured bacterium",
    "D_6__Unknown Species"
  )
  #
  silva_otus <- otu_table(phylo)
  silva_taxa <- tax_table(phylo)
  silva_taxa_df <- as.data.frame(silva_taxa)
  colnames(silva_taxa_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  silva_taxa_df$Phylum[silva_taxa_df$Phylum %in% unassignedOTUs_MDM] <- "<NA>"
  silva_taxa_df$Class[silva_taxa_df$Class %in% unassignedOTUs_MDM] <- "<NA>"
  silva_taxa_df$Order[silva_taxa_df$Order %in% unassignedOTUs_MDM] <- "<NA>"
  silva_taxa_df$Family[silva_taxa_df$Family %in% unassignedOTUs_MDM] <- "<NA>"
  silva_taxa_df$Genus[silva_taxa_df$Genus %in% unassignedOTUs_MDM] <- "<NA>"
  silva_taxa_df$Species[silva_taxa_df$Species %in% unassignedOTUs_MDM] <- "<NA>"
  silva_taxa_df[] <- lapply(silva_taxa_df, as.character)
  silva_taxa_df[] <- lapply(silva_taxa_df, function(x) {
    str_replace_all(x, "uncultured", NA_character_)
  })
  newsilvataxtable2_m <- as.matrix(silva_taxa_df)
  # rownames(newsilvataxtable2_m) <- rownames(silva_updatedtaxa)
  newsilva_taxa <- tax_table(newsilvataxtable2_m)
  new_silva_phylo <- phyloseq(silva_otus, newsilva_taxa, mapping)
}


#' @export
filterTaxonMatrix <- function (x, minocc = 0, dependency = FALSE, keepSum = FALSE,
          return.filtered.indices = FALSE)
{
  toFilter = c()
  xcopy = x
  xcopy[xcopy > 0] = 1
  rowsums = apply(xcopy, 1, sum)
  toFilter = which(rowsums < minocc)
  if (dependency == TRUE) {
    nt = identifyNoisetypes(x, epsilon = 0.5)
    toKeep = c(nt$pink, nt$brown, nt$black)
    toFilter = c(toFilter, setdiff(c(1:nrow(x)), toKeep))
  }
  indices.tokeep = setdiff(c(1:nrow(x)), toFilter)
  if (keepSum == TRUE) {
    filtered = x[toFilter, ]
    x = x[indices.tokeep, ]
    rownames = rownames(x)
    sums.filtered = apply(filtered, 2, sum)
    x = rbind(x, sums.filtered)
    rownames = append(rownames, "summed-nonfeat-rows")
    rownames(x) = rownames
  }
  else {
    x = x[indices.tokeep, ]
  }
  if (return.filtered.indices == TRUE) {
    res = list(x, toFilter)
    names(res) = c("mat", "filtered.indices")
    return(res)
  }
  else {
    return(x)
  }
}


# create net plot at each rank for each environment
# returns 5 ggplots - each of the 5 is a network where nodes are colored by the respective rank (phylum to genus)
#' @export
get_net_plots_all_ranks <- function(orig_graph, orig_phylo, met_name) {
  all_plots_l <- list()
  for (rank_num in seq_along(colnames(tax_table(orig_phylo))[2:6])) {
    print(rank_num)
    tax_rank_names <- c("Phylum", "Class", "Order", "Family", "Genus")
    rank_name <- colnames(tax_table(orig_phylo))[2:6][rank_num]
    print(rank_name)
    rank_name2 <- tax_rank_names[rank_num]
    print(rank_name2)
    g_net_plot <- plot_network(orig_graph, orig_phylo, type = "taxa", color = rank_name, label = NULL, title = paste(met_name, rank_name2, sep = " ")) + theme(legend.position = "none")
    all_plots_l[[rank_name]] <- g_net_plot
  }
  return(all_plots_l)
}


## now do same for hub
# create hub networks, where nodes are sized by hub score and colored by rank
# returns 5 ggplots, each a hub network where nodes are colored by respective rank (phylum to genus)
#' @export
get_hub_plots_all_ranks <- function(orig_graph, orig_phylo, met_name) {
  all_plots_l <- list()
  # for(rank_name in tax_rank_names){
  for (rank_num in seq_along(colnames(tax_table(orig_phylo))[2:6])) {
    print(rank_num)
    tax_rank_names <- c("Phylum", "Class", "Order", "Family", "Genus")
    rank_name <- colnames(tax_table(orig_phylo))[2:6][rank_num]
    print(rank_name)
    rank_name2 <- tax_rank_names[rank_num]
    print(rank_name2)
    g_net_plot <- plot_network(orig_graph, orig_phylo,
      type = "taxa", color = rank_name, label = NULL,
      point_size = hub_score(orig_graph)$vector * 10,
      title = paste(met_name, rank_name2, sep = " ")
    ) + theme(legend.position = "none")
    all_plots_l[[rank_name]] <- g_net_plot
  }
  return(all_plots_l)
}



# use degree_calc_f to get degree, betweenness, and closeness score for each node within each met network
#' @export
degree_calc_f <- function(orig_graph) {
  degree_df <- as.data.frame(degree(orig_graph)) # calculate degree for each node in orig_graph
  degree_df$names <- rownames(degree_df) # keep taxa names of nodes as separate (2nd) column
  colnames(degree_df)[1] <- "degree" # first column renamed as degree for future analyses
  degree_df$bw <- as.vector(betweenness(orig_graph)) # calculate betweenness for each node/store in 3rd column "bw"
  degree_df$closeness <- as.vector(closeness(orig_graph)) # calculate closeness centrality for each node/store in 4th column "closeness"
  print(head(degree_df))
  return(degree_df)
}

#' @export
met_wo_unk <- function(orig_phylo) {
  nodes_list_m <- list()
  met_tax_table <- tax_table(orig_phylo)
  print(dim(met_tax_table))
  met_tax_table <- as.data.frame(met_tax_table)
  met_tax_table[] <- lapply(met_tax_table, as.character)
  met_tax_table1 <- met_tax_table %>% mutate_all(funs(replace(., is.na(.), "Unassigned")))
  rownames(met_tax_table1) <- rownames(met_tax_table)
  print(dim(met_tax_table1))
  for (i in 1:length(colnames(met_tax_table1))) {
    print(i)
    unassignedOTUs_MDM <- c(
      "NA", "D_1__uncultured", "D_1__uncultured bacterium", "D_1__Unknown Phylum",
      "D_2__uncultured", "D_2__uncultured bacterium", "D_2__Unknown Class",
      "Ambiguous_taxa", "Unassigned",
      "D_3__uncultured", "D_3__uncultured bacterium", "D_3__Unknown Order",
      "D_4__uncultured", "D_4__uncultured bacterium",
      "D_4__Unknown Family", "D_5__uncultured", "D_5__uncultured bacterium",
      "D_5__Unknown Genus", "D_6__uncultured",
      "D_6__uncultured bacterium",
      "D_6__Unknown Species"
    )
    node_a <- met_tax_table1[met_tax_table1[, i] %in% unassignedOTUs_MDM, ]
    print(dim(node_a))
    nodes_list_m[[i]] <- node_a
    # print(nodes_list_m[[i]])
  }
  final_otus <- otu_table(orig_phylo)
  print(dim(final_otus))
  wo_unk_otus_at_tax_level_l <- list()
  for (t in 1:length(nodes_list_m)) {
    print(t)
    unk_otus_at_tax_level <- final_otus[!rownames(final_otus) %in% rownames(nodes_list_m[[t]]), ]
    print(dim(unk_otus_at_tax_level))
    wo_unk_otus_at_tax_level_l[[t]] <- unk_otus_at_tax_level
  }
  names(wo_unk_otus_at_tax_level_l) = colnames(tax_table(orig_phylo))
  return(wo_unk_otus_at_tax_level_l)
}

# use get_graph_wo_unk to get back for each environment at each taxonomic level, a graph without unknowns at that level
#' @export
get_graph_wo_unk <- function(orig_graph, met_wo_unk) {
  wo_unk_graph_l <- list()
  for (name in names(met_wo_unk)) {
    print(name)
    tax_level <- name
    wo_unk_graph <- delete_vertices(orig_graph, which(!names(V(orig_graph)) %in% rownames(met_wo_unk[[tax_level]])))
    print(wo_unk_graph)
    wo_unk_graph_l[[name]] <- wo_unk_graph
  }
  return(wo_unk_graph_l)
}

# #use get_degree_df_wo_unk to get degree,betweenness, closeness scores for each of the graphs without unknowns at each taxonomic level
#' @export
get_degree_df_wo_unk <- function(graph_wo_unk){
  met_degree_df_wo_unk_l = list()
  for(i in 1:length(graph_wo_unk)){
    graph_tax_level = names(graph_wo_unk)[[i]]
    print(graph_tax_level)
    graph = graph_wo_unk[[i]]
    graph_degree_df <- degree_calc_f(graph)
    met_degree_df_wo_unk_l[[graph_tax_level]] = graph_degree_df
  }
  return(met_degree_df_wo_unk_l)
}
#
# hs_degree_df_wo_unk <- get_degree_df_wo_unk(hs_graph_wo_unk)
# hy_degree_df_wo_unk <- get_degree_df_wo_unk(hy_graph_wo_unk)
# arc_degree_df_wo_unk <- get_degree_df_wo_unk(arc_graph_wo_unk)
# ds_degree_df_wo_unk <- get_degree_df_wo_unk(ds_graph_wo_unk)

### generate bootstrap networks and compare degree,bw, closeness measures for all 3 net types across all levels
# input= orig_graph, new_wo_unk_graph, orig_phylo, orig_df, degree_df_wo_unk
#' @export
comp_by_deleting_random_knowns_t_v3 <- function(orig_graph, new_wo_unk_graph,
                                                orig_phylo, orig_df,
                                                degree_df_wo_unk, iter = 100, mc.coreval=2) {
  orig_vertices <- as.vector(names(V(orig_graph)))
  measurements <- c("degree", "bw", "closeness")
  all_net_tax_level_plot_l <- do.call(rbind, mclapply(seq_along(new_wo_unk_graph), mc.cores = mc.coreval, function(i) { # for(i in 1:length(new_wo_unk_graph)){
    # i <- 1;
    # for each rank level:
    rank_level <- names(new_wo_unk_graph)[[i]]
    print(paste0("Rank level: ", rank_level))
    num_nodes_to_remove <- length(V(orig_graph)) - length(V(new_wo_unk_graph[[rank_level]]))
    # do bootstrap permutation of node removal
    unk_names <- setdiff(names(V(orig_graph)), names(V(new_wo_unk_graph[[rank_level]]))) # from the original remove the knowns
    orig_vertices_wo_unk <- intersect(names(V(orig_graph)), names(V(new_wo_unk_graph[[rank_level]]))) # the knowns

    # BOOTSTRAP PERMUTATION
    # for each iteration, random knowns are removed (same number as removed when removing unknown OTUs)
    bstrap_df_l <- lapply(seq_len(iter), function(i) { # for (i in  1:iter) {
      # i <- 1;
      vecnames <- sample(x = orig_vertices_wo_unk, size = num_nodes_to_remove, replace = TRUE)
      wo_unk_graph <- delete_vertices(orig_graph, which(names(V(orig_graph)) %in% vecnames))
      boot_degree_df <- degree_calc_f(wo_unk_graph)
      # test out using phylum of hot springs
      if (!identical(length(rownames(boot_degree_df)), length(rownames(degree_df_wo_unk[[rank_level]])))) {
        boot_degree_df <- boot_degree_df[seq_along(rownames(degree_df_wo_unk[[rank_level]])), ]
      }
      boot_degree_df[, measurements]
    })

    do.call(rbind, lapply(measurements, function(act_meas) {
      # act_meas <- measurements[[1]];
      orig_meas <- orig_df[[act_meas]]
      wo_unk_meas <- degree_df_wo_unk[[rank_level]][[act_meas]]
      bstrap_meas <- rowMeans(do.call(cbind, lapply(bstrap_df_l, function(x) x[, act_meas])))
      data.frame(
        rank_level = rank_level,
        type = factor(c(
          rep("Original", length(orig_meas)),
          rep("Without_unknown", length(wo_unk_meas)),
          rep("Bootstrap", length(bstrap_meas))
        ), levels = c("Original", "Without_unknown", "Bootstrap")),
        measure = act_meas,
        data = c(orig_meas, wo_unk_meas, bstrap_meas)
      )
    }))
  }))
  return(all_net_tax_level_plot_l) # returned object should be a list of 7 (for each of classification levels) and for each of the levels, a list of 3 (a boxplot comparison for each measure (degree, bw, closeness))
}

# Get neighbor frequency interactions and create density plots
#' @export
get_bar_and_dens_class_interactions <- function(orig_graph, orig_phylo, met_type) {
  # edgeDF = as_edgelist(hs_graph)
  edgeDF <- as.data.frame(as_edgelist(orig_graph))
  colnames(edgeDF) <- c("X1", "X2")
  # edgeDF_w_alltaxnames <- (merge(tax_table(new_hs_phyloseqobj_final), edgeDF, by.x = "row.names", by.y = "X1"))
  edgeDF_w_alltaxnames <- (merge(tax_table(orig_phylo), edgeDF, by.x = "row.names", by.y = "X1"))
  # edgeDF_w_alltaxnames_updated <- merge(edgeDF_w_alltaxnames, tax_table(new_hs_phyloseqobj_final), by.x = "X2", by.y = "row.names", all.x = TRUE)
  edgeDF_w_alltaxnames_updated <- merge(edgeDF_w_alltaxnames, tax_table(orig_phylo), by.x = "X2", by.y = "row.names", all.x = TRUE)
  print(nrow(edgeDF_w_alltaxnames_updated))
  identical(nrow(edgeDF_w_alltaxnames_updated), nrow(edgeDF))
  # dim(edgeDF_w_alltaxnames_updated) = 552 18 = correct
  colnames(edgeDF_w_alltaxnames_updated)[1:2] <- c("Y", "X")
  edgeDF_w_alltaxnames_updated$Class.x <- as.character(edgeDF_w_alltaxnames_updated$Class.x)
  edgeDF_w_alltaxnames_updated$Class.y <- as.character(edgeDF_w_alltaxnames_updated$Class.y)
  edgeDF_w_alltaxnames_updated[is.na(edgeDF_w_alltaxnames_updated$Class.x), ]$Class.x <- "Unknown"
  edgeDF_w_alltaxnames_updated[is.na(edgeDF_w_alltaxnames_updated$Class.y), ]$Class.y <- "Unknown"
  edgeDF_w_alltaxnames_updated[edgeDF_w_alltaxnames_updated$Y == "0", ]$Class.y <- "Unknown"
  edgeDF_w_alltaxnames_updated$combinedXY <- with(edgeDF_w_alltaxnames_updated, paste0(Class.x, "&", Class.y))

  barplot_of_interactions_l <- list()
  pie_b_plot_interactions_l <- list()
  for (x_class in unique(edgeDF_w_alltaxnames_updated$Class.x)) {
    edge_df_for_class <- as.data.frame(table(edgeDF_w_alltaxnames_updated[grep(x_class, edgeDF_w_alltaxnames_updated$combinedXY), ]$combinedXY))
    # edge_df_for_class = as.data.frame(table(edgeDF_w_alltaxnames_updated[grep(x_class, edgeDF_w_alltaxnames_updated$Class.x),]$combinedXY))
    x_class_at_end <- paste0("&", x_class)
    all_to_rep <- edge_df_for_class$Var1[grep(x_class_at_end, edge_df_for_class$Var1)]
    test_st <- strsplit(as.character(all_to_rep), "&")
    test_st2 <- lapply(test_st, function(x) {
      new_var <- paste0(x_class, "&", x[1])
      return(new_var)
    })
    edge_df_for_class$Var1 <- as.character(edge_df_for_class$Var1)
    edge_df_for_class$Var1[all_to_rep] <- (unlist(test_st2))

    edge_df_for_class <- edge_df_for_class %>% mutate(size = Freq / sum(Freq))
    same_class_var <- paste0(x_class, "&", x_class)
    edge_df_for_class$Var1 <- as.factor(edge_df_for_class$Var1)
    # b_plot <- ggplot(as.data.frame(table(edgeDF_w_alltaxnames_updated[grep(x_class, edgeDF_w_alltaxnames_updated$Class.x),]$combinedXY)), aes(Var1, Freq, fill=Var1)) + geom_bar(stat="identity") + coord_flip() + ggtitle(paste("Neighbor Interactions for " ,x_class, sep = " ")) + theme(legend.position = "none")
    b_plot <- ggplot(edge_df_for_class, aes(Var1, Freq, fill = Var1)) + geom_bar(stat = "identity")
    # b_plot <- b_plot + geom_bar(stat="identity",data=edge_df_for_class[edge_df_for_class$Var1 %in% same_class_var, ], aes(Var1), alpha=0, size=1, color="black") + coord_flip() + theme(legend.position = "none")
    b_plot <- b_plot + geom_bar(stat = "identity", data = edge_df_for_class[edge_df_for_class$Var1 %in% same_class_var, ], aes(Var1), size = 1, fill = "grey50", color = "black") + coord_flip() + theme(legend.position = "none", axis.title.y = element_blank())
    b_plot <- b_plot + ggtitle(paste("Neighbor Interactions for Class ", x_class, sep = " "))

    ## pie chart
    # get back # colors
    # library(viridis)
    # custom_colors = viridis(length(unique(edge_df_for_class$Var1)))
    custom_colors <- rainbow(length(unique(edge_df_for_class$Var1)))
    num_to_change_to_grey <- which(levels(edge_df_for_class$Var1) == same_class_var)
    # num_to_change_to_grey = which(factor(edge_df_for_class$Var1) == same_class_var)
    custom_colors[num_to_change_to_grey] <- "grey50"
    pie_b_plot <- ggplot(edge_df_for_class, aes(x = "", y = size, fill = Var1)) +
      geom_bar(width = 1, stat = "identity", position = position_fill()) + scale_fill_manual(values = custom_colors) +
      coord_polar("y", start = 0) +
      geom_text(aes(label = paste0(round(size * 100), "%")),
        position = position_fill(vjust = 0.6), size = 3
      ) +
      theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank()
      ) + theme(legend.text = element_text(size = 7)) + guides(color = guide_legend(override.aes = list(size = 1)))
    barplot_of_interactions_l[[x_class]] <- b_plot
    pie_b_plot_interactions_l[[x_class]] <- pie_b_plot
  }

  combinedXY_df <- as.data.frame(table(edgeDF_w_alltaxnames_updated$combinedXY))
  val_of_unk1 <- combinedXY_df[which(combinedXY_df$Var1 == "Unknown&Unknown"), ]$Freq
  next_smallest_val1 <- val_of_unk1 - 1
  p_val_sig_dif <- length(which(combinedXY_df$Freq > next_smallest_val1)) / length(combinedXY_df$Freq)
  p_val_sig_dif <- round(p_val_sig_dif, digits = 5)
  d <- density(combinedXY_df$Freq)
  y_axis_max_dens <- max(d$y)

  dens_plot_all_poss <- ggplot(combinedXY_df, aes(Freq)) + geom_density() + geom_vline(xintercept = val_of_unk1, linetype = "dashed", color = "red", size = 1) + ggtitle(paste(met_type, "All Possible Interactions", sep = " "))
  dens_plot_all_poss2 <- dens_plot_all_poss + geom_label(aes(x = val_of_unk1, y = y_axis_max_dens, label = paste("P-value =", p_val_sig_dif, sep = " ")), color = "black", angle = 0, hjust = "inward")
  dens_plot_all_poss2 <- dens_plot_all_poss2 + geom_label(aes(x = val_of_unk1, y = y_axis_max_dens / 2, label = paste("Unknown-Unknown"), hjust = "inward"))

  print(dens_plot_all_poss2)
  combinedXY_just_sameclass_df <- as.data.frame(table(edgeDF_w_alltaxnames_updated[edgeDF_w_alltaxnames_updated$Class.x == edgeDF_w_alltaxnames_updated$Class.y, ]$combinedXY)) # plot just self interactions with same class (ex.Betaproteobacteria-Betaproteobacteria, Unknown-Unknown interactions )
  val_of_unk2 <- combinedXY_just_sameclass_df[which(combinedXY_just_sameclass_df$Var1 == "Unknown&Unknown"), ]$Freq
  next_smallest_val2 <- val_of_unk2 - 1
  p_val_sig_dif2 <- length(which(combinedXY_just_sameclass_df$Freq > next_smallest_val2)) / length(combinedXY_just_sameclass_df$Freq)
  p_val_sig_dif2 <- round(p_val_sig_dif2, digits = 5)
  d2 <- density(combinedXY_just_sameclass_df$Freq)
  y_axis_max_dens2 <- max(d2$y)

  dens_plot_just_sameclass_poss <- ggplot(combinedXY_just_sameclass_df, aes(Freq)) + geom_density() + geom_vline(xintercept = val_of_unk2, linetype = "dashed", color = "red", size = 1) + ggtitle(paste(met_type, "Self-self Class Interactions", sep = " "))
  dens_plot_just_sameclass_poss2 <- dens_plot_just_sameclass_poss + geom_label(aes(x = val_of_unk2, y = y_axis_max_dens2, label = paste("P-value =", p_val_sig_dif2, sep = " ")), color = "black", angle = 0, hjust = "inward")
  dens_plot_just_sameclass_poss2 <- dens_plot_just_sameclass_poss2 + geom_label(aes(x = val_of_unk2, y = y_axis_max_dens2 / 2, label = paste("Unknown-Unknown"), hjust = "inward"))
  print(dens_plot_just_sameclass_poss2)
  # bar_and_dens_plots_l <- list(barplot_of_interactions_l, dens_plot_all_poss2, dens_plot_just_sameclass_poss2)
  # return(bar_and_dens_plots_l)
  bar_pie_and_dens_plots_l <- list(barplot_of_interactions_l, pie_b_plot_interactions_l, dens_plot_all_poss2, dens_plot_just_sameclass_poss2)
  return(bar_pie_and_dens_plots_l)
}

# hs_bar_and_dens <- get_bar_and_dens_class_interactions(hs_graph, new_hs_phyloseqobj_final, "Hot Springs")
# hy_bar_and_dens <- get_bar_and_dens_class_interactions(hy_graph, new_hy_phyloseqobj_final, "Hypersaline")
# ds_bar_and_dens <- get_bar_and_dens_class_interactions(ds_graph, new_ds_phyloseqobj_final, "Deep Sea")
# arc_bar_and_dens <- get_bar_and_dens_class_interactions(arc_graph, new_arc_phyloseqobj_final, "Polar")


# #library(dplyr)
# library(ggsignif)
# library(ggpubr)

#' @export
vis_comp_net_meas_boxplots2 <- function(df, met_name, meas_name_options = NULL, rank_level_options = NULL,
                                        p_label = c("p.signif", "p.format"), p_size = 2) {
  # parameters = df, meas_name_options, rank_level_options, p_label, p_size
  # optional parameters = sample_level_options
  df <- df[df$rank_level != "Kingdom" & df$rank_level != "Species", ]
  if (is.null(meas_name_options) && is.null(rank_level_options)) {
    g1 <- ggplot(df, aes(type, data, color = type)) + facet_grid(measure ~ rank_level, scales = "free") +
      geom_boxplot(outlier.shape = NA) + ylab("Centrality score") +
      stat_compare_means(
        label = p_label, size = p_size,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      ) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      )
  }
  else if (is.null(meas_name_options) && !(is.null(rank_level_options))) {
    spec_df <- df[df$rank_level %in% rank_level_options, ]
    g1 <- ggplot(spec_df, aes(type, data, color = type)) + facet_grid(measure ~ rank_level, scales = "free") +
      geom_boxplot(outlier.shape = NA) + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      stat_compare_means(
        label = p_label, size = p_size,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      )
  }
  else if (!(is.null(meas_name_options)) && is.null(rank_level_options)) {
    spec_df <- df[df$measure %in% meas_name_options, ]
    g1 <- ggplot(spec_df, aes(type, data, color = type)) + facet_grid(measure ~ rank_level, scales = "free") +
      geom_boxplot(outlier.shape = NA) + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      stat_compare_means(
        label = p_label, size = p_size,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      )
  }
  else {
    spec_df <- df[df$measure %in% meas_name_options & df$rank_level %in% rank_level_options, ]
    # if(p_size != 2.5){
    g1 <- ggplot(spec_df, aes(type, data, color = type)) + geom_boxplot(outlier.shape = NA) #+ stat_compare_means(label= p_label, size=p_size,
    # comparisons = list(c("Original", "Without_unknown"),
    #                   c("Without_unknown", "Bootstrap"),
    #                  c("Bootstrap", "Original"))) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    if (length(meas_name_options) > 1 && length(rank_level_options) > 1) {
      g1 <- g1 + facet_grid(measure ~ rank_level, scales = "free") + stat_compare_means(
        label = p_label, size = p_size,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      ) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    else if (length(meas_name_options) > 1 && length(rank_level_options) == 1) {
      g1 <- g1 + facet_wrap(~measure, scales = "free") + stat_compare_means(
        label = p_label, size = p_size,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      ) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    else if (length(meas_name_options) == 1 && length(rank_level_options) > 1) {
      g1 <- g1 + facet_wrap(~rank_level, scales = "free") + stat_compare_means(
        label = p_label, size = p_size,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      ) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
    else if (length(meas_name_options) == 1 && length(rank_level_options) == 1) {
      g1 <- g1 + ylab(meas_name_options) + ggtitle(rank_level_options) + stat_compare_means(
        label = p_label,
        comparisons = list(
          c("Original", "Without_unknown"),
          c("Without_unknown", "Bootstrap"),
          c("Bootstrap", "Original")
        )
      ) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
  }
  g1 <- g1 + ggtitle(paste(met_name, "Centrality Score Comparison", sep = " ")) + ylab("Centrality score") +
    theme(plot.title = element_text(hjust = 0.5))
  return(g1)
}

# test_vis_comp_plot <- vis_comp_net_meas_boxplots(df = all_met_comp_net_df_all_centrality_meas,
#                                                  meas_name_options = c("closeness", "bw"),
#                                                  rank_level_options = c("Family"),
#                                                  p_label = "p.format", p_size=2.5)




#### vis - compare metric or sample
# library(dplyr)
# library(ggpubr)
# library(ggsignif)
#' @export
vis_comp_net_checks <- function(df, perc_sample_options = NULL, metric_options = NULL) {
  spec_df <- df
  if (!is.null(perc_sample_options)) {
    spec_df %>%
      group_by(type, rank_level, measure, percent_samples_included) %>%
      summarise(Median = median(data)) %>%
      as.data.frame() -> spec_df_median

    g1 <- ggplot(spec_df_median, aes(rank_level, Median, color = type)) +
      geom_point(aes(shape = percent_samples_included)) +
      facet_wrap(~ measure, scales = "free") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
      )
  }
  else if (!is.null(metric_options)) {
    spec_df %>%
      group_by(type, rank_level, measure, metric) %>%
      summarise(Median = median(data)) %>%
      as.data.frame() -> spec_df_median

    g1 <- ggplot(spec_df_median, aes(rank_level, Median, color = type)) +
      geom_point(aes(shape = metric)) +
      facet_wrap(~ measure, scales = "free") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
      )
  }
  return(g1)
}

# ########### optional ??? #####################
# ##optional breakdown of MDM ###
# get_MDM_breakdown_plots <- function(orig_phylo){
#create melted dataframe of taxa relative abundance
#' @export
make_mdt <- function(biom_file){
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

# create barplot breakdown plot - # Unassigned, Ambiguous, Uncultured #####
# get_MDM_breakdown_plot <- function(mdt){
#   mdt %>% group_by(SampleID) %>% tally() -> mdt_full_counts_df
#     mdt <- as.data.frame(mdt)
#     test_num_MDM_across_all_tax <- lapply(colnames(mdt)[2:7], function(tax_level){
#       mdt %>% group_by(SampleID) %>% tally() -> mdt_full_counts_df
#       mdt_amb <- mdt[grep("Ambiguous", mdt[,tax_level]),]
#       mdt_unc <- mdt[grep("uncultured", mdt[,tax_level]),]
#       mdt_NA <- mdt[which(is.na(mdt[,tax_level]) == TRUE),]
#       newl <- lapply(list(mdt_amb, mdt_unc, mdt_NA), function(x_df) {
#         if(dim(x_df) != 0){
#           x_df %>% group_by(SampleID) %>% tally() -> x_df2
#         }
#         else{
#           x_df2 = data.frame(SampleID = c("HS", "HY", "DS", "ARC"),
#                              n = rep(0, 4))
#         }
#         return(x_df2)
#       })
#       names(newl) = c("Ambiguous", "Uncultured", "Unassigned")
#       full_MDM_df <- do.call(rbind, newl)
#       full_MDM_df$MDM_Type = rep(c("Ambiguous", "Uncultured", "Unassigned"),each=4)
#       full_MDM_df$taxclasslevel = tax_level
#       #full_MDM_df %>% mutate(sum_all_MDM = Ambiguous.n + Uncultured.n + Unassigned.n) -> full_MDM_df1
#       #full_MDM_df2 <- full_MDM_df1[,c("Ambiguous.SampleID", "sum_all_MDM")]
#       # full_and_MDM_df <- merge(mdt_full_counts_df, full_MDM_df2, by.x="SampleID", by.y="Ambiguous.SampleID", all.x=TRUE)
#       # full_and_MDM_df %>% mutate(K_counts = n - sum_all_MDM) -> full_and_MDM_df2
#       # # full_and_MDM_df2$taxclasslevel = tax_level
#       # colnames(full_and_MDM_df2)[2:4] = c("TotalCounts", "MDM", "K")
#       # print(head(full_and_MDM_df2))
#       return(full_MDM_df)
#     })
#
#     test_num_MDM_all_tax_df <- do.call(rbind, test_num_MDM_across_all_tax)
#
#     test_num_df_combo <- merge(test_num_MDM_all_tax_df, mdt_full_counts_df, by.x = "SampleID", by.y="SampleID", all.x=TRUE)
#     test_num_df_combo_m <- (melt(test_num_df_combo, id.vars = c("SampleID", "MDM_Type", "taxclasslevel", "n.y")))
#
#
#     taxlevelnames = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
#     test_num_df_combo_m$taxclasslevel <- plyr::revalue(test_num_df_combo_m$taxclasslevel, c
#                                                  ("Rank2" = taxlevelnames[1],
#                                                    "Rank3" = taxlevelnames[2],
#                                                    "Rank4" = taxlevelnames[3],
#                                                    "Rank5" = taxlevelnames[4],
#                                                    "Rank6" = taxlevelnames[5],
#                                                    "Rank7" = taxlevelnames[6]))
#     MDMbreakdownplots_all_levels <- lapply(taxlevelnames, function(taxlevel){
#       g <- ggplot(test_num_df_combo_m[test_num_df_combo_m$taxclasslevel == taxlevel,],
#                   aes(SampleID, value, fill=MDM_Type)) + ggtitle(paste(taxlevel, "MDM Breakdown", sep=" ")) +
#         geom_bar(stat="identity", position="stack") + theme_classic() + theme(axis.title.x = element_blank()) +
#         scale_fill_manual(values = c("gold", "navy", "darkred")) + ylab("Number of OTUs")
#       print(g)
#     })
#     return(MDMbreakdownplots_all_levels)
# }
#
#
# ### create lineplot per environment
#first make prevalence df including all taxa information
#' @export
make_prev_df <- function(biom_file){
  phylo <- import_biom(biom_file)
  prevdf = apply(X = otu_table(phylo),MARGIN = ifelse(taxa_are_rows(ps), yes=1, no=2), FUN=function(x){sum(x>0)}) #find prevalence of each taxa per sample
  #bind together taxonomic information and prevalence and total abundance of each taxa into one dataframe
  prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(phylo), tax_table(phylo))
  prevdf <- prevdf[prevdf$Prevalence > 0,] #remove taxa not present in samples
  return(prevdf)
}

# #find prevalence at genus level and create lineplots
#' @export
get_back_counts_for_line_plotf <- function(prev_df, met_name) {
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
####################################################################


### heatmap function for OTU prevalence ####
#' @export
create_abund_heatmap <- function(orig_phylo_w_map, met_type_title, filepathname) {
  manualcolors <- c(
    "black", "forestgreen", "red2", "orange", "cornflowerblue",
    "magenta", "darkolivegreen4",
    "indianred1", "tan4", "darkblue",
    "mediumorchid1", "firebrick4", "yellowgreen", "lightsalmon", "tan3",
    "tan1", "darkgray", "wheat4", "#DDAD4B", "chartreuse", "seagreen1",
    "moccasin", "mediumvioletred", "seagreen", "cadetblue1",
    "darkolivegreen1", "tan2", "tomato3", "#7CE3D8", "gainsboro"
  )
  met_otu_table_mat <- as(otu_table(orig_phylo_w_map)[-nrow(otu_table(orig_phylo_w_map)), ], "matrix")
  met_otu_table_mat_clr <- clr(met_otu_table_mat) # function clr - centered-log-ratio transformation of data - comes from SpiecEasi package
  # change order of otus and samples to match
  dendo_otunames2 <- hclust(dist(met_otu_table_mat_clr)) # hclust(dist) is same function as will be used in heatmap
  dendo_samplenames2 <- hclust(dist(t(met_otu_table_mat_clr)))

  # change otus
  otu_label2 <- as.data.frame(dendo_otunames2$labels)
  met_tax <- as.data.frame(tax_table(orig_phylo_w_map))
  met_tax$ID <- rownames(met_tax)
  otu_label2_df2 <- merge(otu_label2, met_tax, by.x = "dendo_otunames2$labels", by.y = "ID", all.x = TRUE)

  otu_label2_df2[, ] <- lapply(otu_label2_df2[, ], as.character)
  otu_label2_df2[is.na(otu_label2_df2)] <- "Unknown"
  new_otu_colors_clr <- with(otu_label2_df2, data.frame(otu_Class_label = levels(factor(Class)), color = rainbow(length(unique(Class)))))
  otu_label2_df2_clr <- merge(otu_label2_df2, new_otu_colors_clr, by.x = "Class", by.y = "otu_Class_label", all.x = TRUE)
  colnames(otu_label2_df2_clr)[9] <- "otu_colors"
  otu_label2_df2_clr$otu_colors <- as.character(otu_label2_df2_clr$otu_colors)
  otu_label2_df2_clr[otu_label2_df2_clr$Class == "Unknown", ]$otu_colors <- "black"
  otu_label2_df2_clr$org_type <- "white"
  otu_label2_df2_clr[otu_label2_df2_clr$Class == "Unknown", ]$org_type <- "black"
  # otu_sample_col_clr <- t(as.matrix(otu_label2_df2_clr$otu_colors))
  new_otu_sample_col_clr <- as.matrix(rbind(otu_label2_df2_clr$otu_colors, otu_label2_df2_clr$org_type)) # use new_otu_sample_col_clr matrix to assign colors to OTUs (columns)
  rownames(new_otu_sample_col_clr) <- c("Class", "Category") # label colors in heatmap

  # change samples
  sample_label_clr <- as.data.frame(dendo_samplenames2$labels)
  sample_data_from_map <- as.data.frame(as.matrix(sample_data(orig_phylo_w_map)))


  sample_label_df2_clr <- merge(sample_label_clr, sample_data_from_map, by.x = "dendo_samplenames2$labels", by.y = "SampleID", all.x = TRUE)
  num_of_proj <- length(unique(sample_label_df2_clr$ProjectName))
  if (num_of_proj > 11) {
    new_proj_colors_clr <- new_proj_colors_clr <- with(sample_label_df2_clr, data.frame(Project = levels(factor(ProjectName)), color = sample(manualcolors, num_of_proj, replace = FALSE)))
  } else {
    new_proj_colors_clr <- with(sample_label_df2_clr, data.frame(Project = levels(factor(ProjectName)), color = brewer.pal(num_of_proj, "Spectral")))
  }
  new_sample_label_df2_clr <- merge(sample_label_df2_clr, new_proj_colors_clr, by.x = "ProjectName", by.y = "Project", all.x = TRUE)
  colnames(new_sample_label_df2_clr)[16] <- "colors_of_samples"
  new_sample_label_df2_clr$Location <- as.character(new_sample_label_df2_clr$Location)
  hscountrynames <- gsub("\\:.*", "", new_sample_label_df2_clr$Location) # remove everything after : (sublocation) so only Country_name is present
  hscountrynames2 <- gsub("\\_.*", "", hscountrynames) # remove everything after _ so only country is present
  new_sample_label_df2_clr$Country <- hscountrynames2

  num_to_change <- length(unique(new_sample_label_df2_clr$Country))
  new_country_colors_clr <- with(new_sample_label_df2_clr, data.frame(Country = levels(factor(Country)), color = brewer.pal(num_to_change, "Paired")))
  new_sample_label_df2_clr2 <- merge(new_sample_label_df2_clr, new_country_colors_clr, by.x = "Country", by.y = "Country", all.x = TRUE)
  sample_label_colors2_clr <- as.matrix(cbind(as.character(new_sample_label_df2_clr2$colors_of_samples), as.character(new_sample_label_df2_clr2$color))) # use this matrix to assign color to Rows (Samples)
  colnames(sample_label_colors2_clr) <- c("Project", "Country") # label colors of project, country in heatmap

  # now create heatmap with corresponding legends
  par(mar = c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
  #rich8equal <- c("#000041", "#0000CB", "#0081FF", "#02DA81", "#80FE1A", "#FDEE02", "#FFAB00", "#FF3300")
  #png(paste0("~/Tatyana_MDM/", met_type_title, "_", "clr_transformed", ".png"), width = 13, height = 11, units = "in", res = 300)
  pdf(paste0(filepathname, met_type_title, "_", "clr_transformed", ".pdf"))
  #heatmap.3(x = met_otu_table_mat_clr, col = rich8equal, scale = "none", ColSideColors = sample_label_colors2_clr, ColSideColorsSize = 2, RowSideColors = new_otu_sample_col_clr, RowSideColorsSize = 2, keysize = 1, KeyValueName = "OTU Count (after clr)", labCol = FALSE, labRow = FALSE, xlab = "Samples", ylab = "OTUs", main = met_type_title)
  heatmap.3(x = met_otu_table_mat_clr, col = viridis::viridis, scale = "none", ColSideColors = sample_label_colors2_clr, ColSideColorsSize = 2, RowSideColors = new_otu_sample_col_clr, RowSideColorsSize = 2, keysize = 1, KeyValueName = "OTU Count (after clr)", labCol = FALSE, labRow = FALSE, xlab = "Samples", ylab = "OTUs", main = met_type_title)
  # legend("topright", legend = c(unique(as.character(new_sample_label_df2_clr2$ProjectName)), unique(new_sample_label_df2_clr2$Country)), fill =  c(unique(new_sample_label_df2_clr2$colors_of_samples), unique(as.character(new_sample_label_df2_clr2$color))),bty = "n", y.intersp = 0.7, cex=0.7, ncol=2, yjust=0, xjust = 1.5, xpd = TRUE)
  legend("topright", legend = c(unique(new_sample_label_df2_clr2$Country)), fill = c(unique(as.character(new_sample_label_df2_clr2$color))), bty = "n", y.intersp = 0.7, cex = 0.7, ncol = 2, yjust = 0, xjust = 1.5, xpd = TRUE, inset = c(0, -0.02))
  legend("bottomleft", legend = c("Known", "Unknown"), fill = c("White", "Black"), border = TRUE, bty = "n", y.intersp = 0.7, cex = 0.7, xjust = 0)
  dev.off()
}

#### validation checks #####
#### sparcc function and nets (using provided sparcc function from SpiecEasi library)
#' @export
get_sparcc_info <- function(orig_phylo, met_name) {
  # put otu_table of phyloseq obj into community count matrix format for sparcc
  otu_matrix <- as(otu_table(orig_phylo), "matrix")
  print(dim(otu_matrix))
  # transform matrix so that rows are samples, columns are OTUs for sparcc format
  otu_matrix_t <- t(otu_matrix)
  print(dim(otu_matrix_t))
  orig_sparcc <- sparcc(otu_matrix_t)
  # orig_sparcc should now be a list of 2- a symmetric Cov matrix, a symmetric Cor matrix
  # change rownames/colnames of matrix to OTU names
  rownames(orig_sparcc$Cor) <- colnames(orig_sparcc$Cor) <- taxa_names(otu_table(orig_phylo))



  # corr_val_graphy = list()
  # corr_threshold_val <- seq(0.3,0.8,0.1)
  # for(i in 1:length(corr_threshold_val)){
  val <- 0.3
  new_sparcc.graph <- abs(orig_sparcc$Cor) >= val
  diag(new_sparcc.graph) <- 0 # make diagonal of symmetric correlation matrix =0
  new_sparcc.graph <- Matrix(new_sparcc.graph, sparse = TRUE) # put into compatible matrix format for igraph
  new_sparcc_igraph <- adj2igraph(new_sparcc.graph) # convert to igraph
  print(new_sparcc_igraph)
  V(new_sparcc_igraph)$name <- rownames(orig_sparcc$Cor) # change vertex names of graph to OTU names
  print(V(new_sparcc_igraph)$name)

  plot_a <- plot_network(new_sparcc_igraph, orig_phylo, type = "taxa", color = "Class", label = NULL, title = paste(met_name, "SparCC network", sep = " ")) + theme(legend.position = "none")
  plot_b <- plot_network(new_sparcc_igraph, orig_phylo, type = "taxa", color = "Genus", label = NULL, title = paste(met_name, "SparCC hub network", sep = " "), point_size = hub_score(new_sparcc_igraph)$vector * 10) + theme(legend.position = "none")
  all_plots <- list(plot_a, plot_b)
  all_info <- list(new_sparcc_igraph, all_plots)
  return(all_info) # will return list object of 2 - 1st object will be cclasso igraph for further analyses, 2nd obj = list of 2 plots, 1st = cclasso network plot (at class level), 2nd = hub network plot (at genus level)
}

### cclasso function and nets #####
#### use CCLasso function from CCLasso.R ###
#' @export
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3,
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {
  n <- nrow(x)
  p <- ncol(x)

  if (counts) {
    x <- x + pseudo
    x <- x / rowSums(x)
  }
  x <- log(x)
  vx2 <- var(x)

  # Diagonal weight for loss function
  rmean_vx2 <- rowMeans(vx2)
  wd <- 1 / diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2))
  wd2 <- sqrt(wd)

  # Some global parameters for optimization with single lambda
  rho <- 1
  u_f <- eigen(diag(p) - 1 / p)$vectors
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p]
  wd_u_eig <- eigen(wd_u)
  d0_wd <- 1 / ((rep(wd_u_eig$values, each = p - 1) + wd_u_eig$values) /
    (2 * rho) + 1)
  u0_wd <- wd_u_eig$vectors

  # Golden section method for the selection of lambda (log10 scale)
  sigma <- vx2
  lam_int2 <- log10(range(lam_int))
  a1 <- lam_int2[1]
  b1 <- lam_int2[2]
  # Store lambda and corresponding cross validation's loss
  lams <- NULL
  fvals <- NULL
  # Two trial points in first
  a2 <- a1 + 0.382 * (b1 - a1)
  b2 <- a1 + 0.618 * (b1 - a1)
  fb2 <- cv_loss_cclasso(
    lambda2 = 10^b2 / rho, x = x, k_cv = k_cv,
    sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
    wd2 = wd2
  )
  lams <- c(lams, b2)
  fvals <- c(fvals, fb2$cv_loss)
  fa2 <- cv_loss_cclasso(
    lambda2 = 10^a2 / rho, x = x, k_cv = k_cv,
    sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
    wd2 = wd2
  )
  lams <- c(lams, a2)
  fvals <- c(fvals, fa2$cv_loss)
  # Error tolerance for convergence
  err_lam2 <- 1e-1 * max(1, lam_int2)
  err_fval <- 1e-4

  err <- b1 - a1
  k <- 0
  while (err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss)

    if (fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2
      a2 <- b2
      fa2 <- fb2
      b2 <- a1 + 0.618 * (b1 - a1)
      fb2 <- cv_loss_cclasso(
        lambda2 = 10^b2 / rho, x = x, k_cv = k_cv,
        sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
        wd2 = wd2
      )

      lams <- c(lams, b2)
      fvals <- c(fvals, fb2$cv_loss)
    } else {
      b1 <- b2
      b2 <- a2
      fb2 <- fa2
      a2 <- a1 + 0.382 * (b1 - a1)
      fa2 <- cv_loss_cclasso(
        lambda2 = 10^a2 / rho, x = x, k_cv = k_cv,
        sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
        wd2 = wd2
      )

      lams <- c(lams, a2)
      fvals <- c(fvals, fa2$cv_loss)
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss)

    k <- k + 1
    err <- b1 - a1
    if (abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break
    }
  }
  info_cv <- list(
    lams = lams, fvals = fvals, k = k + 2,
    lam_int = 10^c(a1, b1)
  )
  if (a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",",
      10^b1, "])\n",
      sep = ""
    )
  }

  lambda <- 10^((a2 + b2) / 2)
  # Bootstrap for cclasso
  lambda2 <- lambda / rho
  info_boot <- boot_cclasso(
    x = x, sigma = fb2$sigma, lambda2 = lambda2,
    n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd
  )

  return(list(
    var_w = info_boot$var_w, cor_w = info_boot$cor_w,
    p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv
  ))
}
# Bootstrap for cclasso
#' @export
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20,
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x)
  p <- ncol(x)

  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1) / 2, ncol = n_boot + 1)
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1)
  cors_mat <- matrix(0, p, p)
  ind_low <- lower.tri(cors_mat)

  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T),
    ncol = n_boot
  )
  for (k in 1:n_boot) {
    ind_samp <- sam_boot[, k]
    sigma2 <- cclasso_sub(
      sigma = sigma, vx = var(x[ind_samp, ]),
      lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd
    )

    vars_boot[, k] <- diag(sigma2)
    Is <- 1 / sqrt(vars_boot[, k])
    cors_mat <- Is * sigma2 * rep(Is, each = p)
    cors_boot[, k] <- cors_mat[ind_low]
  }
  Is <- 1 / sqrt(diag(sigma))
  cors_mat <- sigma * Is * rep(Is, each = p)
  cors_boot[, n_boot + 1] <- cors_mat[ind_low]
  vars_boot[, n_boot + 1] <- diag(sigma)

  #----------------------------------------
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot)
  #----------------------------------------
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n)
  cors_art0 <- cor(sam_art0)[ind_low]
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)))
  sigma_art <- cclasso_sub(
    sigma = sigma, vx = var(sam_art),
    lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd
  )
  Is <- 1 / sqrt(diag(sigma_art))
  cors_mat <- Is * sigma_art * rep(Is, each = p)
  cors_art2 <- cors_mat[ind_low]
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log(((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) *
    (1 - cors_art0)))
  tmp <- abs(cors_art2) >= tol_cor
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0)
  # Modification of estimation for cclasso
  cors2 <- log((1 + cors_boot) / (1 - cors_boot))
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) +
    (cors_boot <= -tol_cor) * (cors2 - bias02)
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1))
  cors2_mat <- diag(p)
  cors2_mat[ind_low] <- cors2mod
  cors2_mat <- t(cors2_mat)
  cors2_mat[ind_low] <- cors2mod
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2)
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals)
  pval_mat <- diag(p)
  pval_mat[ind_low] <- p_vals
  pval_mat <- t(pval_mat)
  pval_mat[ind_low] <- p_vals
  #----------------------------------------

  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat))
}
#-------------------------------------------------------------------------------
# cross validation's loss of cclasso for single lambda
#' @export
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma,
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x)
  p <- ncol(x)

  n_b <- floor(n / k_cv)
  cv_loss <- 0
  for (k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k)
    vxk <- var(x[itest, ])
    vx2k <- var(x[-itest, ])

    sigma <- cclasso_sub(
      sigma = sigma, vx = vx2k, lambda2 = lambda2,
      wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd
    )

    dsig <- sigma - vxk
    tmp <- rowMeans(dsig)
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp)
    cv_loss <- cv_loss + norm(wd2 * dsig, "F")^2
  }

  return(list(cv_loss = cv_loss, sigma = sigma))
}
#-------------------------------------------------------------------------------
# cclasso for single lambda
#' @export
cclasso_sub <- function(sigma, vx, lambda2,
                        wd, u_f, u0_wd, d0_wd,
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma)
  sigma2 <- sigma
  LAMBDA <- matrix(0, p, p)
  lambda2 <- matrix(lambda2, p, p)
  diag(lambda2) <- 0

  k <- 0
  err <- 1
  while (err > x_tol && k < k_max) {
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f
    x_sigma[-p, -p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) *
      d0_wd) %*% t(u0_wd)
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f)
    # Update sigma2
    A <- LAMBDA + sigma_new
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) *
      (A + lambda2)
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new)

    err <- max(
      abs(sigma_new - sigma) / (abs(sigma) + 1),
      abs(sigma2_new - sigma2) / (abs(sigma2) + 1)
    )
    k <- k + 1
    sigma <- sigma_new
    sigma2 <- sigma2_new
  }

  if (k >= k_max) {
    cat(
      "WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max,
      "&& Relative error:", err, "!\n"
    )
  }

  return(sigma)
}

#######
#' @export
get_cclasso_info <- function(orig_phylo, met_name) {
  hy_cclasso_m <- as.matrix(t(otu_table(orig_phylo)))
  hy_test_cclasso_f <- cclasso(hy_cclasso_m, counts = TRUE)
  hy_cclasso_corr_m <- as.matrix(hy_test_cclasso_f$cor_w)
  rownames(hy_cclasso_corr_m) <- colnames(hy_cclasso_corr_m) <- rownames(otu_table(orig_phylo))
  hy_cclasso_pval_m <- as.matrix(hy_test_cclasso_f$p_vals)
  dim(hy_cclasso_pval_m)
  rownames(hy_cclasso_pval_m) <- colnames(hy_cclasso_pval_m) <- rownames(otu_table(orig_phylo))

  hy_cclasso_cormat1 <- ifelse(hy_cclasso_pval_m <= 0.05, hy_cclasso_corr_m, 0)
  diag(hy_cclasso_cormat1) <- 0
  hy_cclasso_cormat1_corr_thresh_0.6 <- abs(hy_cclasso_cormat1) >= 0.6
  hy_cclasso_cormat1_adjmat <- Matrix(hy_cclasso_cormat1_corr_thresh_0.6, sparse = TRUE)
  library(SpiecEasi)
  hy_cclasso_igraph <- adj2igraph(hy_cclasso_cormat1_adjmat, vertex.attr = list(name = taxa_names(orig_phylo))) # 291 nodes, 477 edges = check
  met_cclasso_graph <- hy_cclasso_igraph
  # V(hy_cclasso_igraph)$name <- rownames(otu_table(new_hy_phyloseqobj_final))

  plot_a <- plot_network(hy_cclasso_igraph, orig_phylo, type = "taxa", color = "Class", label = NULL, title = paste(met_name, "CCLASSO network", sep = " ")) + theme(legend.position = "none")
  plot_b <- plot_network(hy_cclasso_igraph, orig_phylo, type = "taxa", color = "Genus", label = NULL, title = paste(met_name, "CCLASSO hub network", sep = " "), point_size = hub_score(hy_cclasso_igraph)$vector * 10) + theme(legend.position = "none")
  all_plots <- list(plot_a, plot_b)
  all_info <- list(met_cclasso_graph, all_plots)
  return(all_info) # will return list object of 2 - 1st object will be cclasso igraph for further analyses, 2nd obj = list of 2 plots, 1st = cclasso network plot (at class level), 2nd = hub network plot (at genus level)
}

# hs_cclasso_info <- get_cclasso_info(new_hs_phyloseqobj_final, "Hot Springs")
# hy_cclasso_info <- get_cclasso_info(new_hy_phyloseqobj_final, "Hypersaline")
# ds_cclasso_info <- get_cclasso_info(new_ds_phyloseqobj_final, "Deep Sea")
# arc_cclasso_info <- get_cclasso_info(new_arc_phyloseqobj_final, "Polar")

##### pearson function and nets #####
# function to get back pearson network and pearson hub network plots for each environment
#' @export
get_pearson_info <- function(orig_phylo, met_name) {
  print(met_name)
  hy_test_pearson <- rcorr(as.matrix(t(otu_table(orig_phylo))), type = "pearson")
  hy_test_pearson_cormat <- as.matrix(hy_test_pearson$r)
  hy_test_pearson_pmat <- as.matrix(hy_test_pearson$P)
  hy_test_pearson_pmat[is.na(hy_test_pearson_pmat)] <- 1
  hy_test_pearson_cormat1 <- ifelse(hy_test_pearson_pmat <= 0.05, hy_test_pearson_cormat, 0)


  hy_test_pearson_cormat2 <- abs(hy_test_pearson_cormat1) >= 0.6
  hy_test_pearson_cormat2[1:10, 1:10]
  hy_test_pearson_adjmat <- Matrix(hy_test_pearson_cormat2, sparse = TRUE)


  # library(SpiecEasi)
  hy_test_pearson_igraph <- adj2igraph(hy_test_pearson_adjmat, vertex.attr = list(name = taxa_names(orig_phylo))) # 291 nodes, 477 edges = check
  # V(hy_cclasso_igraph)$name <- rownames(otu_table(new_hy_phyloseqobj_final))
  print(hy_test_pearson_igraph)
  plot_a <- plot_network(hy_test_pearson_igraph, orig_phylo, type = "taxa", color = "Class", label = NULL, title = paste(met_name, "Pearson network", sep = " ")) + theme(legend.position = "none")
  plot_b <- plot_network(hy_test_pearson_igraph, orig_phylo, type = "taxa", color = "Genus", label = NULL, title = paste(met_name, "Pearson hub network", sep = " "), point_size = hub_score(hy_test_pearson_igraph)$vector * 10) + theme(legend.position = "none")
  all_plots <- list(plot_a, plot_b)
  all_info <- list(hy_test_pearson_igraph, all_plots)
  return(all_info)
}

# hs_pearson_info <- get_pearson_info(new_hs_phyloseqobj_final, "Hot Springs")
# hy_pearson_info <- get_pearson_info(new_hy_phyloseqobj_final, "Hypersaline")
# ds_pearson_info <- get_pearson_info(new_ds_phyloseqobj_final, "Deep Sea")
# arc_pearson_info <- get_pearson_info(new_arc_phyloseqobj_final, "Polar")

### use pearson_info[[1]] to access and plot igraph object ####
#### function to plot all networks in same shape/manner as original (SpiecEasi) ###
## modify plot_network function from phyloseq to change layout to reference network style ###
#' @export
plotnet_comp_v2 <- function(g, g_orig_for_layout, physeq = NULL, type = "samples", color = NULL, shape = NULL,
                            point_size = 4, alpha = 1, label = "value", hjust = 1.35,
                            line_weight = 0.5, line_color = color, line_alpha = 0.4,
                            layout.method = layout.fruchterman.reingold, title = NULL) {
  if (vcount(g) < 2) {
    stop("The graph you provided, `g`, has too few vertices. \n         Check your graph, or the output of `make_network` and try again.")
  }
  if (type %in% c("taxa", "species", "OTUs", "otus", "otu")) {
    type <- "taxa"
  }
  edgeDF <- data.frame(get.edgelist(g))
  edgeDF$id <- 1:length(edgeDF[, 1])
  vertDF <- layout.method(g_orig_for_layout)
  colnames(vertDF) <- c("x", "y")
  vertDF <- data.frame(
    value = get.vertex.attribute(g, "name"),
    vertDF
  )
  if (!is.null(physeq)) {
    extraData <- NULL
    if (type == "samples" & !is.null(sample_data(
      physeq,
      FALSE
    ))) {
      extraData <- data.frame(sample_data(physeq))[as.character(vertDF$value), ,
        drop = FALSE
      ]
    }
    else if (type == "taxa" & !is.null(tax_table(
      physeq,
      FALSE
    ))) {
      extraData <- data.frame(tax_table(physeq))[as.character(vertDF$value), ,
        drop = FALSE
      ]
    }
    if (!is.null(extraData)) {
      vertDF <- data.frame(vertDF, extraData)
    }
  }
  graphDF <- merge(melt(edgeDF, id = "id"), vertDF,
    by = "value"
  )
  p <- ggplot(vertDF, aes(x, y))
  p <- p + theme_bw() + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(), axis.ticks = element_blank(),
    panel.border = element_blank()
  )
  p <- p + geom_point(aes_string(color = color, shape = shape),
    size = point_size, na.rm = TRUE
  )
  if (!is.null(label)) {
    p <- p + geom_text(aes_string(label = label),
      size = 2,
      hjust = hjust, na.rm = TRUE
    )
  }
  p <- p + geom_line(aes_string(group = "id", color = line_color),
    graphDF,
    size = line_weight, alpha = line_alpha, na.rm = TRUE
  )
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

##heatmap.3 function for heatmaps
#' @export
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){

  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }

  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)

    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }

    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }

  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }

  if (!missing(ColSideColors)) {

    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }

  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }

    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

# #all hs plots
# hs_graph_plot_newcomp <- plotnet_comp_v2(hs_graph, new_hs_phsloseqobj_final, g_orig_for_layout = hs_graph, type='taxa', color="Class", label = NULL, title="HS SpiecEasi")
# hs_sparcc_plot_newcomp <- plotnet_comp_v2(hs_sparcc_graph[[1]], new_hs_phyloseqobj_final, g_orig_for_layout = hs_graph, type='taxa', color="Class", label = NULL, title="HS SparCC")
# hs_pearson_plot_newcomp <- plotnet_comp_v2(hs_pearson_info[[1]], new_hs_phyloseqobj_final, g_orig_for_layout = hs_graph, type='taxa', color="Class", label = NULL, title="HS Pearson")
# hs_cclasso_plot_newcomp <- plotnet_comp_v2(hs_cclasso_info[[1]], new_hs_phyloseqobj_final, g_orig_for_layout = hs_graph, type='taxa', color="Class", label = NULL, title="HS CCLasso")

# package.skeleton(name="MDMAnalyzerR_v0", path= ".")

# usethis::use_package("phyloseq")
# usethis::use_package("igraph")
# usethis::use_package("ggplot2")
# usethis::use_package("dplyr")
# usethis::use_package("data.table")
# usethis::use_package("ggpubr")
# usethis::use_package("Hmisc")
# usethis::use_package("seqtime")

# usethis::use_package("Matrix")
# usethis::use_package("magrittr")
# use_data(unassignedOTUs_MDM, internal=TRUE)
# use_data(tax_rank_names, internal=TRUE, overwrite = TRUE)
# use_data(hs_graph, internal=TRUE, overwrite = TRUE)
# use_data(new_hs_phyloseqobj_final, internal=TRUE, overwrite = TRUE)
# use_data(new_hs_phyloseqobj_final)
# use_data(hs_graph)
#use_data(hs_phylo_w_map)

#### #' @importFrom seqtime filterTaxonMatrix ### -removed line
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
