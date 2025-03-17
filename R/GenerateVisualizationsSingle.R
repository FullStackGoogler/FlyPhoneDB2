#' Heatmap generation for Multi Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param DEG The file path for the differentially expressed genes file
#' @param pathwayObj The core component file used in the analysis
#' @param pathways The pathways of interest. Leave this unspecified if you want to generate heatmaps for all pathways in pathwayObj.
#'
#' @export
HeatmapMulti <- function(output_dir = NULL, DEG, pathwayObj, pathways = NULL) {
  DEG_OBJ <- read.csv(DEG, row.names = 1)
  interested_pathways <- if(!is.null(pathways)) unique(pathways) else unique(pathwayObj$pathway)

  for(pathway in interested_pathways) {
    if(!pathway %in% pathwayObj$pathway) {
      stop(paste0("Error: '", pathway, "' is not in pathwayObj!"))
    }
  }

  # Reshape data: Spread 'cell_type' into columns with 'avg_log2FC' as values
  heatmap_data <- DEG_OBJ %>%
    select(gene, cell_type, avg_log2FC) %>%
    spread(key = cell_type, value = avg_log2FC)

  # Convert to matrix, with genes as rownames
  rownames(heatmap_data) <- heatmap_data$gene
  heatmap_data$gene <- NULL  # Remove gene column after setting as rownames
  heatmap_matrix <- as.matrix(heatmap_data)

  # Define the color palette mapping
  palette_length <- 50  # Number of colors in the palette
  my_breaks <- seq(-3, 3, length.out = palette_length + 1)  # Breaks from -4 to 4

  # Create a color palette
  my_colors <- colorRampPalette(c("blue", "white", "red"))(palette_length)

  for(signaling_path in interested_pathways) {
    signaling_subset <- pathwayObj %>% filter(pathway == signaling_path) %>% select(2,4)
    signaling_subset$role <- as.factor(signaling_subset$role)
    signaling_subset$role <- factor(signaling_subset$role, levels = c("ligand", "receptor", "reporter", "TF", "other"))
    signaling_subset <- signaling_subset %>% arrange(role)
    heatmap_matrix_subset <- heatmap_data %>% filter(rownames(heatmap_data) %in% unique(signaling_subset$gene))
    heatmap_matrix_subset$gene <- rownames(heatmap_matrix_subset)
    heatmap_matrix_subset <- as.matrix(heatmap_matrix_subset)
    heatmap_matrix_subset[is.na(heatmap_matrix_subset)] <- 0
    df <- as.data.frame(as.table(heatmap_matrix_subset))
    colnames(df) <- c("gene", "cell_type", "expression")
    df$expression <- as.numeric(as.character(df$expression))
    merged_df <- df %>%
      left_join(signaling_subset, by = "gene") %>% filter(! cell_type == "gene")

    # Plotting the heatmap with ggplot2
    p <- ggplot(merged_df, aes(x = cell_type, y = gene, fill = expression)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey50"
      ) +
      #theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing.y = unit(2, "lines"),  # Increase space between panels
        strip.placement = "outside",  # Moves strip labels outside of panel
        strip.text.y.left = element_text(angle = 0), # Role label on the right
        strip.text.y.right = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", color = NA)
      ) +
      labs(title = signaling_path, x = "Cell Type", y = "Gene", fill = "Log2FC") +
      facet_grid(rows = vars(role), scales = "free_y", space = "free_y", switch = "y")

    if(!is.null(output_dir)) {
      ggsave(p, file = paste0(output_dir,"/", gsub("[_/, ]", "-", signaling_path), ".png"),
             width = 12, # The width of the plot in inches
             height = 8)
    } else {
      print(p)
    }
  }
}

#' Heatmap generation for Single Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param seuratObj A Seurat Object where counts = raw count data, and meta.data = metadata
#' @param pathwayObj The core component file used in the analysis
#' @param pathways The pathways of interest. Leave this unspecified if you want to generate heatmaps for all pathways in pathwayObj.
#'
#' @export
HeatmapSingle <- function(output_dir = NULL, seuratObj, pathwayObj, pathways = NULL) {
  interested_pathways <- if(!is.null(pathways)) unique(pathways) else unique(pathwayObj$pathway)

  for(pathway in interested_pathways) {
    if(!pathway %in% pathwayObj$pathway) {
      stop(paste0("Error: '", pathway, "' is not in pathwayObj!"))
    }
  }

  seuratObj <- NormalizeData(seuratObj)

  # Set colors of clusters
  color_palette <- NULL

  celltype_count <- length(sort(unique(seuratObj$cluster)))

  if(celltype_count < 35) { # Use polychrome palette; up to 34 colors
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome") # Save palette information
    color_palette <- color_palette[3:36]
  } else { # Too many celltypes, use varibow instead
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = cell_type_count, palette = "varibow", shuffle_pal = TRUE) # Save palette information
  }

  names(color_palette) = sort(unique(seuratObj$cluster))

  # Get average expression of each cluster for all genes
  avgexp <- AverageExpression(seuratObj, assay = "RNA", return.seurat = T, group.by = c("cluster"))

  factor_levels <- avgexp@active.ident

  for(i in interested_pathways) {
    cat("\n")
    cat("## ", i, " {.tabset} \n")
    df <- subset(pathwayObj, pathway == i)
    genes <- df$gene
    cat("\n")
    p <- DoHeatmap(avgexp,
                   features = genes,
                   label = TRUE,
                   draw.lines = FALSE,
                   raster = FALSE,
                   angle = 90,
                   size = 3,
                   slot = "scale.data",
                   group.bar = TRUE,
                   group.colors = color_palette) +
      scale_fill_gradientn(colors = c("#313695","#FFFFBF","#A50026")) +
      guides(color = "none") +
      theme(
        plot.margin = margin(t = 55, r = 10, b = 10, l = 10), # Adjust margins
      )

    if(!is.null(output_dir)) {
      ggsave(p, file = paste0(output_dir,"/", gsub("[_/, ]", "-", i), ".png"),
             width = 12, # The width of the plot in inches
             height = 8)
    } else {
      print(p)
    }
  }
}

#' Dot plot generation for Single Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param seuratObj A Seurat Object where counts = raw count data, and meta.data = metadata
#' @param pathwayObj The core component file used in the analysis
#' @param pathways The pathways of interest. Leave this unspecified if you want to generate heatmaps for all pathways in pathwayObj.
#'
#' @export
DotPlotSingle <- function(output_dir = NULL, seuratObj, pathwayObj, pathways = NULL) {
  interested_pathways <- if(!is.null(pathways)) unique(pathways) else unique(pathwayObj$pathway)

  for(pathway in interested_pathways) {
    if(!pathway %in% pathwayObj$pathway) {
      stop(paste0("Error: '", pathway, "' is not in pathwayObj!"))
    }
  }

  seuratObj <- NormalizeData(seuratObj)

  # Set colors of clusters
  ncelltype <- length(unique(seuratObj$cluster))
  polychrome_pal <- scCustomize::DiscretePalette_scCustomize(num_colors = (ncelltype+2), palette = "polychrome")
  polychrome_pal <- polychrome_pal[3:(ncelltype+2)]
  names(polychrome_pal) = sort((unique(seuratObj$cluster)))

  for(signaling_path in interested_pathways) {
    signaling_subset <- pathwayObj %>% filter(pathway == signaling_path) %>% select(2,4)
    signaling_subset$role <- as.factor(signaling_subset$role)
    signaling_subset$role <- factor(signaling_subset$role, levels = c("ligand", "receptor", "reporter", "TF", "other"))
    signaling_subset <- signaling_subset %>% arrange(role)
    a <- DotPlot(seuratObj, features = unique(signaling_subset$gene), group.by = "cluster")[['data']]
    a <- a %>% tibble::remove_rownames() %>% rename(gene = features.plot)
    a <- a %>% full_join(signaling_subset)
    a <- a %>% select(id, gene, role, avg.exp.scaled, pct.exp)
    a <- na.omit(a)
    # Create a dot plot
    p <- ggplot(a, aes(x = id, y = gene)) +
      geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      facet_grid(rows = vars(role), scales = "free_y", space = "free_y") +
      #theme_minimal() +
      labs(color = "Average Expression\n(Scaled)", size = "Percent Expressed") +
      xlab("Celltype") + ylab("Gene") + ggtitle(signaling_path) +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(angle = 0, size = 10, face = "bold"),
            strip.placement = "outside",  # Moves strip labels outside of panel
            strip.text.y.left = element_text(angle = 0), # Role label on the right
            strip.text.y.right = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "lightgrey", color = NA),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black")
      )

    if(!is.null(output_dir)) {
      ggsave(p, file = paste0(output_dir,"/", gsub("[_/, ]", "-", signaling_path), ".png"),
             width = 12,
             height = 12,
             dpi = 300)
    } else {
      print(p)
    }
  }
}

#' Circle plot generation for Multi Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param results The file path to the results Excel file generated by RunMultiAnalysis()
#' @param pathwayObj The core component file used in the analysis
#' @param pathways The pathways of interest. Leave this unspecified if you want to generate heatmaps for all pathways in pathwayObj.
#' @param cellTypes The celltypes of interest. Leave this unspecified if you want to explore all cell types.
#'
#' @export
CirclePlotMulti <- function(output_dir = NULL, results, pathwayObj, pathways = NULL, cellTypes = NULL) {
  interested_pathways <- if(!is.null(pathways)) unique(pathways) else unique(pathwayObj$pathway)

  for(pathway in interested_pathways) {
    if(!pathway %in% pathwayObj$pathway) {
      stop(paste0("Error: '", pathway, "' is not in pathwayObj!"))
    }
  }

  df1 <- read_excel(results, sheet = 1)
  df2 <- read_excel(results, sheet = 2)
  df3 <- read_excel(results, sheet = 3)

  df2$score_control <- 0
  df3$score_mutant <- 0

  interaction_df <- rbind(df1, df2, df3)
  interaction_df$secretor <- gsub("/", "_", interaction_df$secretor) #FIXME? Mainly "ISC/EB" is problem
  interaction_df$receptor <- gsub("/", "_", interaction_df$receptor) #FIXME? Mainly "ISC/EB" is problem

  interaction_list_control <- interaction_df %>% select(c(Gene_secreted, Gene_receptor, pathway_receptor, secretor, receptor, score_control, pval_control))
  colnames(interaction_list_control)[6] <- "control_score"
  interaction_list_mutant <- interaction_df %>% select(c(Gene_secreted, Gene_receptor, pathway_receptor, secretor, receptor, score_mutant, pval_mutant))
  colnames(interaction_list_mutant)[6] <- "mutant_score"

  # Assign color to each cell type
  paste0(unique(append(interaction_df$secretor, interaction_df$receptor)))

  # Set colors of clusters
  color_palette <- NULL

  celltype_count <- length(sort(unique(append(interaction_df$secretor, interaction_df$receptor))))

  if(celltype_count < 35) { # Use polychrome palette; up to 34 colors
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome") # Save palette information
    color_palette <- color_palette[3:36]
  } else { # Too many celltypes, use varibow instead
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = cell_type_count, palette = "varibow", shuffle_pal = TRUE) # Save palette information
  }
  names(color_palette) <- sort(unique(append(interaction_df$secretor, interaction_df$receptor)))
  celltypes <- names(color_palette)
  celltypes <- na.omit(names(color_palette))
  celltypes <- gsub("/", "_", celltypes) #FIXME? Mainly "ISC/EB" is problem

  for(celltype in cellTypes) {
    if(!celltype %in% celltypes) {
      stop(paste0("Error: '", celltype, "' not found in your results file!"))
    }
  }

  celltypes_to_check <- if(!is.null(cellTypes)) unique(cellTypes) else celltypes

  for(celltype in celltypes_to_check) {
    # Create a new workbook for the summed interaction scores
    wb <- createWorkbook()

    # Loop through each pathway
    for(pathway in interested_pathways) {
      cat(paste0(pathway, "\n"))

      # Process mutant interactions
      interaction_pathway_long_mutant <- subset(interaction_list_mutant, pathway_receptor == pathway) %>%
        filter(secretor == celltype & receptor != celltype) %>%
        filter(mutant_score > 0) %>%
        mutate(variable = paste0(secretor, ".", receptor)) %>%
        group_by(variable) %>%
        summarize(sum_score = sum(mutant_score))

      # Process control interactions
      interaction_pathway_long_control <- subset(interaction_list_control, pathway_receptor == pathway) %>%
        filter(secretor == celltype & receptor != celltype) %>%
        filter(control_score > 0) %>%
        mutate(variable = paste0(secretor, ".", receptor)) %>%
        group_by(variable) %>%
        summarize(sum_score = sum(control_score))

      combined_df <- full_join(interaction_pathway_long_control,
                               interaction_pathway_long_mutant,
                               by = "variable",
                               suffix = c("_control", "_mutant"))

      combined_df <- combined_df %>%
        mutate(
          sum_score_control = coalesce(sum_score_control, 0),
          sum_score_mutant = coalesce(sum_score_mutant, 0)
        )

      net <- combined_df %>%
        mutate(score_difference = sum_score_mutant - sum_score_control) %>%
        select(variable, score_difference)

      # Convert to sender/receiver format
      net$variable <- gsub("\\.", ">", net$variable)
      net <- net %>%
        tidyr::separate(variable, c("sender", "receiver"), ">") %>%
        mutate(pair = if_else(sender < receiver, paste0(sender, "/", receiver), paste0(receiver, "/", sender))) %>%
        group_by(pair) %>%
        filter(score_difference == max(score_difference, na.rm = TRUE)) %>% #There was 1 warning in `filter()`. In ℹ In argument: `score_difference == max(score_difference)`. Caused by warning in `max()`: ! no non-missing arguments to max; returning -Inf
        ungroup() %>%
        select(-pair) %>%
        filter(score_difference != 0)

      # Complete graph setup
      empty_celltype <- setdiff(celltypes, unique(c(net$sender, net$receiver)))
      for (ct in empty_celltype) {
        line <- c(ct, ct, 0)
        net <- rbind(net, line)
      }

      colnames(net) <- c("sender", "receiver", "n")
      net$n <- as.numeric(net$n)

      if(sum(net$n) == 0){
        next
      }

      # Create graph
      g <- graph_from_data_frame(net, directed = TRUE)
      x <- as_adjacency_matrix(g, attr = "n", sparse = FALSE)
      x <- x[celltypes, celltypes]
      g <- graph_from_adjacency_matrix(x, mode = "directed", weighted = TRUE)

      # Define color mapping for edge coloring
      E(g)$color <- ifelse(E(g)$weight > 0, "red", "blue")

      # Determine visual properties
      edge.start <- ends(g, es = E(g), names = FALSE)
      coords <- layout_(g, in_circle())
      coords_scale <- scale(coords)

      V(g)$size <- 20
      V(g)$color <- color_palette[V(g)$name]  # Use the corrected color palette
      V(g)$label.color <- "black"
      V(g)$label.cex <- 0.8
      if (max(E(g)$weight) == min(E(g)$weight)) {
        E(g)$width <- 1
      } else {
        E(g)$width <- 0.5 + abs(E(g)$weight) / max(abs(E(g)$weight)) * 4
      }
      E(g)$arrow.width <- 3
      E(g)$label.color <- 'black'

      # Plot setup and save
      if(!is.null(output_dir)) {
        png(file = paste0(output_dir, "/", gsub("[_/, ]", "-", celltype), "_", gsub("[_/, ]", "-", pathway), ".png"),
            width = 10,
            height = 10,
            units = "in",
            res = 300)

        plot(g, edge.curved = 0.2, vertex.shape = 'circle',
             layout = coords_scale, margin = 0.2, edge.arrow.size = 0.5,
             vertex.frame.color = "white", label = FALSE)

        dev.off()

        # Assign pathway to summed interaction score table
        net$pathway <- pathway

        # Add a worksheet to the workbook
        addWorksheet(wb, pathway)

        # Write the dataframe to the sheet
        writeData(wb, sheet = pathway, net)

        cat(paste0("Done! Moving on... \n"))
      } else {
        plot_title <- paste("Pathway:", pathway, "- Cell type:", celltype)

        plot(g, edge.curved = 0.2, vertex.shape = 'circle',
             layout = coords_scale, margin = 0.2, edge.arrow.size = 0.5,
             vertex.frame.color = "white", label = FALSE, main = plot_title)
        cat(paste0("Plot displayed for ", pathway, " in cell type ", celltype, "\n"))
      }
    }

    if (!is.null(output_dir)) {
      saveWorkbook(wb, paste0(output_dir, "/", gsub("[_/, ]", "-", celltype), "_summed-interaction-scores.xlsx"), overwrite = TRUE)
    }
  }
}

#' Circle plot generation for Single Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param results The file path to the results file (long or long_filtered (recommended)) generated by RunSingleAnalysis()
#' @param pathwayObj The core component file used in the analysis
#' @param pathways The pathways of interest. Leave this unspecified if you want to generate heatmaps for all pathways in pathwayObj.
#' @param cellTypes The celltypes of interest. Leave this unspecified if you want to explore all cell types.
#'
#' @export
CirclePlotSingle <- function(output_dir = NULL, results, pathwayObj, pathways = NULL, cellTypes = NULL) {
  interested_pathways <- if(!is.null(pathways)) unique(pathways) else unique(pathwayObj$pathway)

  for(pathway in interested_pathways) {
    if(!pathway %in% pathwayObj$pathway) {
      stop(paste0("Error: '", pathway, "' is not in pathwayObj!"))
    }
  }

  interaction_list <- read.csv(results)
  filtered = grepl("filtered", results)
  interaction_list$secretor <- gsub("/", "_", interaction_list$secretor) #FIXME? Mainly "ISC/EB" is problem
  interaction_list$receptor <- gsub("/", "_", interaction_list$receptor) #FIXME? Mainly "ISC/EB" is problem

  sample_name <- if(filtered) gsub("^.*interaction-long-filtered_(.*)\\.csv$", "\\1", results) else gsub("^.*interaction-long_(.*)\\.csv$", "\\1", results)
  sample_name <- gsub("[_/, ]", "-", sample_name)
  name_extension <- if(!sample_name == "") paste0(sample_name, "_") else ""

  # Print cell types to std out, assign color to each cell type
  paste0(unique(append(interaction_list$secretor, interaction_list$receptor)))

  # Set colors of clusters
  color_palette <- NULL

  celltype_count <- length(sort(unique(append(interaction_list$secretor, interaction_list$receptor))))

  if(celltype_count < 35) { # Use polychrome palette; up to 34 colors
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome") # Save palette information
    color_palette <- color_palette[3:36]
  } else { # Too many celltypes, use varibow instead
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = cell_type_count, palette = "varibow", shuffle_pal = TRUE) # Save palette information
  }
  names(color_palette) <- sort(unique(append(interaction_list$secretor, interaction_list$receptor)))
  celltypes <- na.omit(names(color_palette))
  celltypes <- gsub("/", "_", celltypes) #FIXME? Mainly "ISC/EB" is problem

  for(celltype in cellTypes) {
    if(!celltype %in% celltypes) {
      stop(paste0("Error: '", celltype, "' not found in your results file!"))
    }
  }

  celltypes_to_check <- if(!is.null(cellTypes)) unique(cellTypes) else celltypes

  for(celltype in celltypes_to_check) {
    wb <- createWorkbook()

    for(pathway in interested_pathways) {

      cat(paste0(pathway, "\n"))

      interaction_pathway_long <- subset(interaction_list, pathway_receptor == pathway) %>%
        filter(secretor == celltype & receptor != celltype) %>%
        filter(score > 0) %>%
        mutate(variable = paste0(secretor, ".", receptor)) %>%
        select(4,5,6,8)

      data <- interaction_pathway_long

      data$variable <- gsub(".", ">", data$variable, fixed = TRUE)

      col <- color_palette
      label=FALSE
      edge.curved=0.5
      shape='circle'
      layout=in_circle()
      vertex.size=20
      margin=0.2
      vertex.label.cex=0.8
      vertex.label.color='black'
      arrow.width=3
      edge.label.color='black'
      edge.label.cex=1
      edge.max.width=4 # the maximum thickness of the line is 4

      colnames(data)[3] <- "score"

      net <- data %>% group_by(variable) %>% summarize(sum_score = sum(score)) # calculate the average score group by cell type

      net$variable <- gsub(".", ">", net$variable, fixed = TRUE)

      net <- net %>%
        tidyr::separate(variable, c("sender", "receiver"), ">")

      net <- net %>%
        mutate(pair = if_else(sender < receiver, paste0(sender, "/", receiver), paste0(receiver, "/", sender)))

      net <- net %>%
        group_by(pair) %>%
        filter(sum_score == max(sum_score)) %>%
        ungroup() %>%
        select(-pair)

      net <- net %>%
        filter(sum_score > 0.0)

      empty_celltype <- setdiff(celltypes, unique(c(net$sender, net$receiver)))

      for(ct in empty_celltype) {
        line <- c(ct, ct, 0)
        net <- rbind(net, line)
      }

      colnames(net) <- c("sender", "receiver", "n")
      net$n <- as.numeric(net$n)

      # This code chunk may be changed to only show the top 15 signals (uncomment slice_head(n=15) %>%)
      top_val <- net %>%
        arrange(desc(n)) %>%
        slice_head(n = 50) %>%
        pull(n)

      net <- net %>%
        mutate(n = ifelse(n %in% top_val, n, 0))

      net<-as.data.frame(net,stringsAsFactors=FALSE)
      g<-graph_from_data_frame(net,directed=TRUE)
      x <- as_adjacency_matrix(g, attr="n", sparse=FALSE)
      x <- x[celltypes, celltypes]
      g <- graph_from_adjacency_matrix(x, mode = "directed", weighted = T)
      edge.start <- ends(g, es=E(g), names=FALSE)
      coords<-layout_(g,layout)

      if(sum(net$n) == 0){
        next
      }

      if(nrow(coords)!=1){
        coords_scale=scale(coords)
      }else{
        coords_scale<-coords
      }

      loop.angle<-ifelse(coords_scale[V(g),1]>0,-atan(coords_scale[V(g),2]/coords_scale[V(g),1]),pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
      V(g)$size<-vertex.size
      V(g)$color<-col[V(g)]
      V(g)$label.color<-vertex.label.color
      V(g)$label.cex<-vertex.label.cex
      if(label){
        E(g)$label<-E(g)$n
      }

      if(max(E(g)$weight)==min(E(g)$weight)){
        E(g)$width<-1 # if all the average scores are the same, set all the line thickness to 1
      }else{
        E(g)$width <- 0.5 + E(g)$weight/max(E(g)$weight)*edge.max.width # otherwise, set the line thickness linearly related to the average score
      }
      E(g)$arrow.width<-arrow.width
      E(g)$label.color<-edge.label.color
      E(g)$color<-V(g)$color[edge.start[,1]]

      if(!is.null(output_dir)) {
        png(file=paste0(output_dir, "/", name_extension, gsub("[_/, ]", "-", celltype), "_", gsub(" ","-",pathway), ".png"),
            width = 10,
            height = 10,
            units = "in",
            res = 300)

        plot(g,edge.curved=0.2,vertex.shape=shape,
             layout=coords_scale,margin=margin,edge.arrow.size=0.5, vertex.frame.color="white"
             , label=FALSE)

        dev.off()

        # Assign pathway to summed interaction score table
        net$pathway <- pathway

        # Add a worksheet to the workbook
        addWorksheet(wb, pathway)

        # Write the dataframe to the sheet
        writeData(wb, sheet = pathway, net)

        cat(paste0("Done! Moving on... \n"))
      } else {
        plot_title <- paste("Pathway:", pathway, "- Cell type:", celltype)

        plot(g,edge.curved=0.2,vertex.shape=shape,
             layout=coords_scale,margin=margin,edge.arrow.size=0.5, vertex.frame.color="white"
             , label=FALSE, main = plot_title)

        cat(paste0("Plot displayed for ", pathway, " in cell type ", celltype, "\n"))
      }
    }

    if(!is.null(output_dir)) {
      # Save the workbook
      saveWorkbook(wb, paste0(output_dir, "/", gsub("[_/, ]", "-", celltype),"_summed-interaction-scores.xlsx"), overwrite = TRUE)
    }
  }
}

#' Chord diagram generation for Multi Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param results The file path to the results Excel file generated by RunMultiAnalysis()
#' @param cellTypes The celltypes of interest. Leave this unspecified if you want to explore all cell types.
#'
#' @export
ChordDiagramMulti <- function(output_dir = NULL, results, cellTypes = NULL) {
  df1 <- read_excel(results, sheet = 1)
  df2 <- read_excel(results, sheet = 2)
  df3 <- read_excel(results, sheet = 3)

  df2$score_control <- 0
  df3$score_mutant <- 0

  interaction_df <- rbind(df1, df2, df3)
  interaction_df$secretor <- gsub("/", "_", interaction_df$secretor) #FIXME? Mainly "ISC/EB" is problem
  interaction_df$receptor <- gsub("/", "_", interaction_df$receptor) #FIXME? Mainly "ISC/EB" is problem

  interaction_list <- interaction_df %>% select(c(Gene_secreted, Gene_receptor, pathway_receptor, secretor, receptor, score_control, pval_control))
  interaction_list2 <- interaction_df %>% select(c(Gene_secreted, Gene_receptor, pathway_receptor, secretor, receptor, score_mutant, pval_mutant))

  paste0(unique(append(interaction_list$secretor, interaction_list$receptor)))

  # Set colors of clusters
  color_palette <- NULL

  celltype_count <- length(sort(unique(append(interaction_list$secretor, interaction_list$receptor))))

  if(celltype_count < 35) { # Use polychrome palette; up to 34 colors
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome") # Save palette information
    color_palette <- color_palette[3:36]
  } else { # Too many celltypes, use varibow instead
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = cell_type_count, palette = "varibow", shuffle_pal = TRUE) # Save palette information
  }
  names(color_palette) <- sort(unique(append(interaction_list$secretor, interaction_list$receptor)))

  celltypes <- unique(interaction_list$secretor)

  for(celltype in cellTypes) {
    if(!celltype %in% celltypes) {
      stop(paste0("Error: '", celltype, "' not found in your results file!"))
    }
  }

  sheet_names <- excel_sheets(results)
  mutant_name <- gsub("Significant_In_", "", sheet_names[2])
  control_name <- gsub("Significant_In_", "", sheet_names[3])

  celltypes_to_check <- if(!is.null(cellTypes)) unique(cellTypes) else celltypes

  for(celltype in celltypes_to_check) {
    interaction_list_1 <- interaction_list
    interaction_list_2 <- interaction_list2

    colnames(interaction_list_1)[6:7] <- c("score", "pval")
    colnames(interaction_list_2)[6:7] <- c("score", "pval")

    # Process data for the specific celltype
    interaction_list_subset <- interaction_list_1 %>%
      filter(secretor == celltype & receptor != celltype)

    if(nrow(interaction_list_subset) == 0){
      return(NULL)
    }

    # Get summarized interaction scores
    sorting_order <- interaction_list_subset %>%
      group_by(receptor) %>%
      summarize(total_score = sum(score)) %>%
      arrange(desc(total_score)) %>%
      ungroup()

    # Sort dataframe for chord diagram
    sorted_chord_df <- interaction_list_subset
    sorted_chord_df$receptor <- factor(sorted_chord_df$receptor, levels = sorting_order$receptor)
    sorted_chord_df <- sorted_chord_df %>% left_join(sorting_order) %>% arrange(desc(total_score)) %>% select(-total_score)
    sorted_chord_df <- sorted_chord_df[,c(5,4,6)]

    # Process data for the specific celltype
    interaction_list_subset_2 <- interaction_list_2 %>%
      filter(secretor == celltype & receptor != celltype)

    if(nrow(interaction_list_subset_2) == 0){
      return(NULL)
    }

    # Get summarized interaction scores
    sorting_order_2 <- interaction_list_subset_2 %>%
      group_by(receptor) %>%
      summarize(total_score = sum(score)) %>%
      arrange(desc(total_score)) %>%
      ungroup()

    # Sort dataframe for chord diagram
    sorted_chord_df_2 <- interaction_list_subset_2
    sorted_chord_df_2$receptor <- factor(sorted_chord_df_2$receptor, levels = sorting_order_2$receptor)
    sorted_chord_df_2 <- sorted_chord_df_2 %>% left_join(sorting_order_2) %>% arrange(desc(total_score)) %>% select(-total_score)
    sorted_chord_df_2 <- sorted_chord_df_2[,c(5,4,6)]
    ###

    if(sum(sorted_chord_df$score) > sum(sorted_chord_df_2$score)){
      gap = calc_gap(sorted_chord_df, sorted_chord_df_2, big.gap = 10, small.gap = 0.1)
    }else{
      gap = calc_gap(sorted_chord_df_2, sorted_chord_df, big.gap = 10, small.gap = 0.05)
    }

    #   case_when(sum(sorted_chord_df$score) > sum(sorted_chord_df_2$score) ~ calc_gap(sorted_chord_df, sorted_chord_df_2, big.gap = 30, small.gap = 3),
    #             sum(sorted_chord_df$score) < sum(sorted_chord_df_2$score) ~ calc_gap(sorted_chord_df_2, sorted_chord_df, big.gap = 30, small.gap = 3))
    # # gap = calc_gap(sorted_chord_df, sorted_chord_df_2, big.gap = 30, small.gap = 3)

    # Clear previous circos plots
    circos.clear()
    circos.par(circle.margin = c(0.5, 0.5, 0.1, 0.1))

    # Open a temporary graphics device
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, width=3000, height=3000, res=300)

    # Draw the chord diagram
    chordDiagramFromDataFrame(
      sorted_chord_df,
      grid.col = color_palette,
      big.gap = gap,
      small.gap = 1,
      annotationTrack = c("grid"),
      preAllocateTracks = list(track.height = 0.4)
    )

    circos.track(track.index = 1,
                 panel.fun = function(x, y) {
                   circos.text(CELL_META$xcenter, CELL_META$ylim[1] + cm_h(0.75), CELL_META$sector.index,
                               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                 }, bg.border = NA)
    title(main = control_name)
    dev.off()  # Close the graphics device

    # Convert the PNG to a rasterGrob
    image <- png::readPNG(temp_file)
    unlink(temp_file)  # Remove temporary file
    f1 <- rasterGrob(image)

    # SECOND PLOT #

    # Clear previous circos plots
    circos.clear()
    circos.par(circle.margin = c(0.5, 0.5, 0.1, 0.1))

    # Open a temporary graphics device
    temp_file_2 <- tempfile(fileext = ".png")
    png(temp_file_2, width=3000, height=3000, res=300)

    # Draw the chord diagram
    chordDiagramFromDataFrame(
      sorted_chord_df_2,
      grid.col = color_palette,
      big.gap = gap,
      small.gap = 1,
      annotationTrack = c("grid"),
      preAllocateTracks = list(track.height = 0.4)
    )

    circos.track(track.index = 1,
                 panel.fun = function(x, y) {
                   circos.text(CELL_META$xcenter, CELL_META$ylim[1] + cm_h(0.75), CELL_META$sector.index,
                               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                 }, bg.border = NA)
    title(main = mutant_name)
    dev.off()  # Close the graphics device

    # Convert the PNG to a rasterGrob
    image_2 <- png::readPNG(temp_file_2)
    unlink(temp_file_2)  # Remove temporary file
    f2 <- rasterGrob(image_2)
    p <- cowplot::plot_grid(f1, f2, nrow = 1)  # Arrange side by side, left is always control

    if(!is.null(output_dir)) {
      ggsave(filename = paste0(output_dir, "/", gsub("[_/, ]", "-", celltype), "_chord-diagram.png"), width = 6000, height = 3000, units = "px", )
    } else {
      print(p)
    }
  }
}

#' Chord diagram generation for Single Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param results The file path to the results file (long or long_filtered (recommended)) generated by RunSingleAnalysis()
#' @param cellTypes The celltypes of interest. Leave this unspecified if you want to explore all cell types.
#'
#' @export
ChordDiagramSingle <- function(output_dir = NULL, results, cellTypes = NULL) {
  interaction_list <- read.csv(results)
  filtered = grepl("filtered", results)
  interaction_list$secretor <- gsub("/", "_", interaction_list$secretor) #FIXME? Mainly "ISC/EB" is problem
  interaction_list$receptor <- gsub("/", "_", interaction_list$receptor) #FIXME? Mainly "ISC/EB" is problem

  sample_name <- if(filtered) gsub("^.*interaction-long-filtered_(.*)\\.csv$", "\\1", results) else gsub("^.*interaction-long_(.*)\\.csv$", "\\1", results)
  sample_name <- gsub("[_/, ]", "-", sample_name)
  name_extension <- if(!sample_name == "") paste0(sample_name, "_") else ""

  # Print cell types to std out, assign color to each cell type
  paste0(unique(append(interaction_list$secretor, interaction_list$receptor)))

  # Set colors of clusters
  color_palette <- NULL

  celltype_count <- length(sort(unique(append(interaction_list$secretor, interaction_list$receptor))))

  if(celltype_count < 35) { # Use polychrome palette; up to 34 colors
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome") # Save palette information
    color_palette <- color_palette[3:36]
  } else { # Too many celltypes, use varibow instead
    color_palette <- scCustomize::DiscretePalette_scCustomize(num_colors = cell_type_count, palette = "varibow", shuffle_pal = TRUE) # Save palette information
  }
  names(color_palette) <- sort(unique(append(interaction_list$secretor, interaction_list$receptor)))
  celltypes <- na.omit(names(color_palette))
  celltypes <- gsub("/", "_", celltypes) #FIXME? Mainly "ISC/EB" is problem

  for(celltype in cellTypes) {
    if(!celltype %in% celltypes) {
      stop(paste0("Error: '", celltype, "' not found in your results file!"))
    }
  }

  celltypes_to_check <- if(!is.null(cellTypes)) unique(cellTypes) else celltypes

  # Find the maximum width required for cell types
  longest_name_length <- max(nchar(celltypes))
  # Compute the width and height multiplier based on the longest name
  width_height_multiplier <- longest_name_length / 20

  for(celltype in celltypes_to_check) {
    print(paste0("Making diagram for celltype: ", celltype))

    # Process only celltype of interest
    interaction_list_cell <- interaction_list %>% filter(secretor == celltype & receptor != celltype)

    # Get summarized interaction scores
    sorting_order <- interaction_list_cell %>%
      group_by(receptor) %>%
      summarize(total_score = sum(score)) %>%
      arrange(desc(total_score)) %>%
      ungroup()

    # Sort dataframe in descending order and prepare dataframe for chord diagram function
    sorted_chord_df <- interaction_list_cell
    sorted_chord_df$receptor <- factor(sorted_chord_df$receptor, levels = sorting_order$receptor)
    sorted_chord_df <- sorted_chord_df %>% left_join(sorting_order) %>% arrange(desc(total_score)) %>% select(-total_score)
    sorted_chord_df <- sorted_chord_df[,c(5,4,6)]

    # Check if output_dir is not NULL
    if(!is.null(output_dir)) {
      png(paste0(output_dir, "/", name_extension, gsub("[_/, ]", "-", celltype), "_chord-diagram.png"),
          width = 10,
          height = 10,
          units = "in",
          res = 300)

      # Clear any previous circos plots
      circos.clear()
      circos.par(circle.margin=c(0.5, 0.5, 0.1, 0.1))

      # Define a reasonable estimate for the track height
      label_width <- 0.3  # Fixed width, adjust based on visual inspection and need

      # Draw the chord diagram
      chordDiagramFromDataFrame(
        sorted_chord_df,
        grid.col = color_palette,
        big.gap = 30,
        small.gap = 4,
        annotationTrack = c("grid"),
        preAllocateTracks = list(track.height = label_width)  # Use the fixed label width
      )

      # Add text within track
      circos.track(track.index = 1,
                   panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ylim[1] + cm_h(0.75), CELL_META$sector.index,
                                 facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                   }, bg.border = NA)

      dev.off()
    } else {
      # Clear any previous circos plots
      circos.clear()
      circos.par(circle.margin=c(0.5, 0.5, 0.1, 0.1))

      # Define a reasonable estimate for the track height
      label_width <- 0.3  # Fixed width, adjust based on visual inspection and need

      # Draw the chord diagram
      chordDiagramFromDataFrame(
        sorted_chord_df,
        grid.col = color_palette,
        big.gap = 30,
        small.gap = 4,
        annotationTrack = c("grid"),
        preAllocateTracks = list(track.height = label_width)  # Use the fixed label width
      )

      # Add text within track
      circos.track(track.index = 1,
                   panel.fun = function(x, y) {
                     circos.text(CELL_META$xcenter, CELL_META$ylim[1] + cm_h(0.75), CELL_META$sector.index,
                                 facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                   }, bg.border = NA)
    }
  }
}

#' Interaction Strength scatterplot generation for Multi Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param results The file path to the results Excel file generated by RunMultiAnalysis()
#'
#' @export
InteractionStrengthMulti <- function(output_dir = NULL, results) {
  df1 <- read_excel(results, sheet = 1)
  df2 <- read_excel(results, sheet = 2)
  df3 <- read_excel(results, sheet = 3)

  df2$score_control <- 0
  df3$score_mutant <- 0

  interaction_df <- rbind(df1, df2, df3)
  interaction_df$secretor <- gsub("/", "_", interaction_df$secretor)
  interaction_df$receptor <- gsub("/", "_", interaction_df$receptor)

  sheet_names <- excel_sheets(results)
  mutant_name <- gsub("Significant_In_", "", sheet_names[2])
  control_name <- gsub("Significant_In_", "", sheet_names[3])

  control_df <- interaction_df %>% select(c(Gene_secreted, Gene_receptor, pathway_receptor, secretor, receptor, score_control, pval_control))
  mutant_df <- interaction_df %>% select(c(Gene_secreted, Gene_receptor, pathway_receptor, secretor, receptor, score_mutant, pval_mutant))

  # Calculate scores for control
  outgoing_control <- control_df %>%
    group_by(secretor) %>%
    summarize(outgoing_score = sum(score_control, na.rm = TRUE))

  incoming_control <- control_df %>%
    group_by(receptor) %>%
    summarize(incoming_score = sum(score_control, na.rm = TRUE))

  scores_control <- full_join(outgoing_control, incoming_control, by = c("secretor" = "receptor")) %>%
    mutate(cell_type = secretor)
  scores_control[is.na(scores_control)] <- 0

  # Calculate scores for mutant
  outgoing_mutant <- mutant_df %>%
    group_by(secretor) %>%
    summarize(outgoing_score = sum(score_mutant, na.rm = TRUE))

  incoming_mutant <- mutant_df %>%
    group_by(receptor) %>%
    summarize(incoming_score = sum(score_mutant, na.rm = TRUE))

  scores_mutant <- full_join(outgoing_mutant, incoming_mutant, by = c("secretor" = "receptor")) %>%
    mutate(cell_type = secretor)
  scores_mutant[is.na(scores_mutant)] <- 0

  # Determine combined axis limits
  x_range_combined <- range(c(scores_control$incoming_score, scores_mutant$incoming_score), na.rm = TRUE)
  y_range_combined <- range(c(scores_control$outgoing_score, scores_mutant$outgoing_score), na.rm = TRUE)

  x_limits <- c(floor(x_range_combined[1] * 0.9), ceiling(x_range_combined[2] * 1.1))
  y_limits <- c(floor(y_range_combined[1] * 0.9), ceiling(y_range_combined[2] * 1.1))

  # Define custom color palette
  n_colors <- max(length(unique(scores_control$cell_type)), length(unique(scores_mutant$cell_type)))
  palette_colors <- NULL

  if(n_colors < 35) { # Use polychrome palette; up to 34 colors
    palette_colors <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome") # Save palette information
    palette_colors <- palette_colors[3:36]
  } else { # Too many celltypes, use varibow instead
    palette_colors <- scCustomize::DiscretePalette_scCustomize(num_colors = n_colors, palette = "varibow", shuffle_pal = TRUE) # Save palette information
  }

  # Plot for control
  a1 <- ggplot(scores_control, aes(x = incoming_score, y = outgoing_score)) +
    geom_point(aes(color = cell_type), size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = cell_type, color = cell_type), size = 3) +
    scale_color_manual(values = palette_colors) +
    labs(title = paste0("Interaction Strength: ", control_name), x = "Incoming Score", y = "Outgoing Score", color = "Cell Type") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    ggpubr::theme_pubr() +
    theme(legend.position = "right")

  # Plot for mutant
  a2 <- ggplot(scores_mutant, aes(x = incoming_score, y = outgoing_score)) +
    geom_point(aes(color = cell_type), size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = cell_type, color = cell_type), size = 3) +
    scale_color_manual(values = palette_colors) +
    labs(title = paste0("Interaction Strength: ", mutant_name), x = "Incoming Score", y = "Outgoing Score", color = "Cell Type") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    ggpubr::theme_pubr() +
    theme(legend.position = "right")

  # Combine the plots using patchwork
  combined_plot <- a1 / a2 + plot_layout(guides = 'collect')

  if(!is.null(output_dir)) {
    ggsave(paste0(output_dir, "/Interaction_Strength_DEG_filtered.png"), combined_plot, width = 12, height = 10, dpi = 300)
  } else {
    print(combined_plot)
  }
}

#' Interaction Strength scatterplot generation for Single Analysis data
#'
#' @param output_dir Directory to save generated visualizations to. If unspecified, visualizations will not be saved and instead will simply be printed out.
#' @param resultFileList The file path to the results file (long or long_filtered (recommended)) generated by RunSingleAnalysis(). If multiple files are inputted, InteractionStrengthSingle() will scale all axes to be the same.
#'
#' @export
InteractionStrengthSingle <- function(output_dir = NULL, resultFileList) {
  all_outgoing_scores_list <- list()
  all_incoming_scores_list <- list()

  for (i in seq_along(resultFileList)) {
    file <- unlist(resultFileList[i])
    interaction_df <- read.csv(file)

    outgoing <- interaction_df %>%
      group_by(secretor) %>%
      summarize(outgoing_score = sum(score, na.rm = TRUE))

    incoming <- interaction_df %>%
      group_by(receptor) %>%
      summarize(incoming_score = sum(score, na.rm = TRUE))

    all_outgoing_scores_list[[i]] <- outgoing
    all_incoming_scores_list[[i]] <- incoming
  }

  combined_outgoing <- bind_rows(all_outgoing_scores_list)
  combined_incoming <- bind_rows(all_incoming_scores_list)

  # Calculate the ranges for the axis limits across all dataframes
  x_range <- range(combined_incoming$incoming_score, na.rm = TRUE)
  y_range <- range(combined_outgoing$outgoing_score, na.rm = TRUE)

  # Add padding to the limits for better visualization
  x_limits <- c(floor(x_range[1] * 0.9), ceiling(x_range[2] * 1.1))
  y_limits <- c(floor(y_range[1] * 0.9), ceiling(y_range[2] * 1.1))

  # Plot each dataset with common scales
  for (i in seq_along(resultFileList)) {
    file <- unlist(resultFileList[i])
    filtered = grepl("filtered", file)
    print(filtered)
    print(file)
    sample_name <- if(filtered) gsub(".*/interaction-long-filtered_(.*)\\.csv$", "\\1", file) else gsub(".*/interaction-long_(.*)\\.csv$", "\\1", file)
    print(sample_name)
    sample_name <- gsub("[_/, ]", "-", sample_name)
    scores <- full_join(all_outgoing_scores_list[[i]], all_incoming_scores_list[[i]], by = c("secretor" = "receptor"))
    scores[is.na(scores)] <- 0
    scores <- scores %>%
      mutate(cell_type = secretor)

    n_colors <- length(unique(scores$cell_type))
    palette_colors <- scales::hue_pal()(n_colors)

    plot_title <- paste("Interaction Strength:", sample_name)

    p <- ggplot(scores, aes(x = incoming_score, y = outgoing_score)) +
      geom_point(aes(color = cell_type), size = 3, alpha = 0.7) +
      geom_text_repel(aes(label = cell_type, color = cell_type),
                      size = 3, box.padding = 0.3, point.padding = 0.2, max.overlaps = 10) +
      theme_minimal(base_size = 18) +
      scale_color_manual(values = palette_colors) +
      labs(title = plot_title,
           x = "Incoming Score",
           y = "Outgoing Score",
           color = "Cell Type") +
      scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
      scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
      ggpubr::theme_pubr() +
      ggpubr::labs_pubr() +
      theme(legend.position = "right")

    if(!is.null(output_dir)) {
      output_filename <- paste0(output_dir, "/", sample_name, "_interaction-strength.png")
      ggsave(output_filename, plot = p, width = 12, height = 10, dpi = 300)
    } else {
      print(p)
    }
  }
}







