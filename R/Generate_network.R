#' Generate KEGG Pathway Network Visualization
#'
#' Creates an interactive Plotly network graph showing relationships between
#' significantly changed metabolites and their associated KEGG pathways, modules,
#' enzymes, and reactions.
#'
#' @param i Integer. Index of the contrast to visualize.
#' @param mydata List. Network data structure containing graph objects.
#' @param dir Character. Direction indicator ("up" or "down").
#' @param dir_word Character. Human-readable direction ("more" or "less").
#' @param contrasts Character vector. List of contrasts from statistical analysis.
#' @param test_type Character. Statistical test type ("t.test", "anova", etc.).
#' @param results_df Data frame. Results table containing KEGG IDs and confidence levels.
#'
#' @return Plotly network visualization object.
#'
#' @details
#' The network includes:
#' - **Pathways** (nodes colored by type)
#' - **Significantly changed compounds** (highlighted in red, shown as squares)
#' - **Related enzymes, reactions, and modules** (colored by category)
#'
#' Node shapes:
#' - Circle: Non-input metabolites
#' - Square: Input metabolites (significantly changed)
#'
#' Node colors:
#' - Dark2 palette: Pathway, module, enzyme, reaction, compound
#' - Red: Significantly changed input compounds
#'
#' @export

generate_network <- function(i, mydata, dir, dir_word, 
                              contrasts, test_type, results_df) {
  
  # Required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network generation.")
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots.")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for color palettes.")
  }
  if (!requireNamespace("stringi", quietly = TRUE)) {
    stop("Package 'stringi' is required for string manipulation.")
  }
  
  library(igraph)
  library(plotly)
  library(RColorBrewer)
  library(stringi)
  
  # ===========================
  # Extract Graph Components
  # ===========================
  
  # Get node list
  mydata[[paste("vs", dir, sep = ".")]][[i]] <- 
    V(mydata[[paste(dir, "graph", sep = ".")]][[i]])
  
  # Get edge list
  mydata[[paste("es", dir, sep = ".")]][[i]] <- 
    as.data.frame(get.edgelist(mydata[[paste(dir, "graph", sep = ".")]][[i]]))
  
  # Get node attribute dataframe
  mydata[[paste("node.data", dir, sep = ".")]][[i]] <- 
    get.data.frame(mydata[[paste(dir, "graph", sep = ".")]][[i]], what = "vertices")
  
  # ===========================
  # Process Node Attributes
  # ===========================
  
  node_data <- mydata[[paste("node.data", dir, sep = ".")]][[i]] %>%
    mutate(input = as.character(input)) %>%
    mutate(com = case_when(
      endsWith(input, "TRUE") ~ "6",
      TRUE ~ as.character(com)
    )) %>%
    mutate(input = case_when(
      endsWith(input, "TRUE") ~ "Input",
      endsWith(input, "FALSE") ~ "",
    ))
  
  # Recode node categories
  if (test_type == "nostats") {
    node_data <- node_data %>%
      mutate_at(vars(com), ~ stri_replace_all_regex(
        .,
        pattern = as.character(1:6),
        replacement = c("pathway", "module", "enzyme", "reaction", 
                       "compound", ">1.5 log2FC changed compound"),
        vectorize = FALSE
      ))
  } else {
    node_data <- node_data %>%
      mutate_at(vars(com), ~ stri_replace_all_regex(
        .,
        pattern = as.character(1:6),
        replacement = c("pathway", "module", "enzyme", "reaction", 
                       "compound", "Significantly changed compound"),
        vectorize = FALSE
      ))
  }
  
  mydata[[paste("node.data", dir, sep = ".")]][[i]] <- node_data
  
  # Count nodes and edges
  mydata[[paste("Nv", dir, sep = ".")]][[i]] <- 
    length(mydata[[paste("vs", dir, sep = ".")]][[i]])
  mydata[[paste("Ne", dir, sep = ".")]][[i]] <- 
    length(mydata[[paste("es", dir, sep = ".")]][[i]][[1]])
  
  # ===========================
  # Calculate Layout
  # ===========================
  
  mydata[[paste("L", dir, sep = ".")]][[i]] <- 
    layout.fruchterman.reingold(mydata[[paste(dir, "graph", sep = ".")]][[i]])
  
  mydata[[paste("Xn", dir, sep = ".")]][[i]] <- 
    mydata[[paste("L", dir, sep = ".")]][[i]][, 1]
  mydata[[paste("Yn", dir, sep = ".")]][[i]] <- 
    mydata[[paste("L", dir, sep = ".")]][[i]][, 2]
  
  # ===========================
  # Load KEGG Mapping Data
  # ===========================
  
  # FIXED: Use relative path from project root
  kegg_map_path <- file.path("reference", "KEGG_map.df.RData")
  
  if (!file.exists(kegg_map_path)) {
    stop(sprintf("KEGG mapping file not found: %s", kegg_map_path))
  }
  
  KEGG_map <- readRDS(kegg_map_path)
  
  # ===========================
  # Merge Node Names with KEGG Data
  # ===========================
  
  mydata[[paste("vs.names", dir, ".")]][[i]] <- 
    data.frame(names(mydata[[paste("vs", dir, sep = ".")]][[i]])) %>%
    mutate(num = as.numeric(row.names(.)))
  
  colnames(mydata[[paste("vs.names", dir, ".")]][[i]]) <- c("vs.name", "num")
  
  mydata[[paste("vs.names.KEGG", dir, ".")]][[i]] <- 
    merge(mydata[[paste("vs.names", dir, ".")]][[i]], KEGG_map, 
          by.x = "vs.name", by.y = "V1", all.x = TRUE) %>%
    arrange(num)
  
  # ===========================
  # Add Confidence Levels
  # ===========================
  
  mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]] <- 
    results_df %>%
    filter(contrast == contrasts[[i]]) %>%
    select(KEGG, ID_confidence) %>%
    filter(KEGG != "") %>%
    distinct()
  
  mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]] <- 
    merge(
      mydata[[paste("vs.names.KEGG", dir, ".")]][[i]],
      mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]],
      by.x = "vs.name", by.y = "KEGG"
    ) %>%
    mutate(ID_confidence = case_when(
      ID_confidence == "High" ~ "High Confidence ID",
      ID_confidence == "Low" ~ "Low Confidence ID",
      .default = as.character(ID_confidence)
    )) %>%
    mutate(report.results.Level = paste(V2, ID_confidence, sep = ";"))
  
  # Handle duplicates: prioritize high confidence
  mydata[[paste("vs.names.KEGG.confidence.2", dir, ".")]][[i]] <- 
    mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]] %>%
    filter(!vs.name %in% .[duplicated(.$vs.name), ]$vs.name)
  
  mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]] <- rbind(
    mydata[[paste("vs.names.KEGG.confidence.2", dir, ".")]][[i]],
    mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]] %>%
      filter(vs.name %in% .[duplicated(.$vs.name), ]$vs.name) %>%
      filter(ID_confidence == "High Confidence ID")
  )
  
  mydata[[paste("vs.names.KEGG", dir, ".")]][[i]] <- 
    merge(
      mydata[[paste("vs.names.KEGG", dir, ".")]][[i]],
      mydata[[paste("vs.names.KEGG.confidence", dir, ".")]][[i]],
      all.x = TRUE
    ) %>%
    mutate(V2 = case_when(
      !is.na(report.results.Level) ~ report.results.Level,
      .default = V2
    )) %>%
    select(-report.results.Level) %>%
    arrange(num)
  
  # Add line breaks for hover text
  mydata[[paste("vs.names.KEGG", dir, ".")]][[i]]$V2 <- 
    gsub(";", "\n", mydata[[paste("vs.names.KEGG", dir, ".")]][[i]]$V2)
  
  # ===========================
  # Create Network Plot
  # ===========================
  
  # Adjust node sizes (input nodes are larger)
  node_data$size <- ifelse(node_data$input == "Input", 10, 5)
  
  # Initialize Plotly figure
  network <- plot_ly(evaluate = TRUE)
  
  # Add nodes
  network <- add_trace(
    network,
    x = ~mydata[[paste("Xn", dir, sep = ".")]][[i]],
    y = ~mydata[[paste("Yn", dir, sep = ".")]][[i]],
    mode = "markers",
    text = mydata[[paste("vs.names.KEGG", dir, ".")]][[i]]$vs.name,
    hoverinfo = "text",
    hovertext = mydata[[paste("vs.names.KEGG", dir, ".")]][[i]]$V2,
    color = as.factor(node_data$com),
    colors = c(brewer.pal(5, "Dark2"), "red"),
    symbol = ~node_data$input,
    symbols = c("circle", "square"),
    size = ~node_data$size
  )
  
  # Add node labels
  network <- network %>%
    add_text(
      x = ~mydata[[paste("Xn", dir, sep = ".")]][[i]],
      y = ~mydata[[paste("Yn", dir, sep = ".")]][[i]],
      text = mydata[[paste("vs.names.KEGG", dir, ".")]][[i]]$vs.name,
      evaluate = TRUE
    )
  
  # ===========================
  # Create Edges
  # ===========================
  
  names(mydata[[paste("Xn", dir, sep = ".")]][[i]]) <- 
    names(mydata[[paste("vs", dir, sep = ".")]][[i]])
  names(mydata[[paste("Yn", dir, sep = ".")]][[i]]) <- 
    names(mydata[[paste("vs", dir, sep = ".")]][[i]])
  
  # Generate edge shapes
  calculate_edge_shape <- function(j) {
    v0 <- as.character(mydata[[paste("es", dir, sep = ".")]][[i]][j, ]$V1)
    v1 <- as.character(mydata[[paste("es", dir, sep = ".")]][[i]][j, ]$V2)
    
    list(
      type = "line",
      line = list(color = "red", width = 0.3),
      x0 = mydata[[paste("Xn", dir, sep = ".")]][[i]][v0],
      y0 = mydata[[paste("Yn", dir, sep = ".")]][[i]][v0],
      x1 = mydata[[paste("Xn", dir, sep = ".")]][[i]][v1],
      y1 = mydata[[paste("Yn", dir, sep = ".")]][[i]][v1]
    )
  }
  
  edge_shapes_list <- lapply(
    1:mydata[[paste("Ne", dir, sep = ".")]][[i]],
    calculate_edge_shape
  )
  
  mydata[[paste("edge_shapes", dir, sep = ".")]][[i]] <- edge_shapes_list
  
  # Define axis settings
  mydata[[paste("axis", dir, sep = ".")]][[i]] <- list(
    title = "",
    showgrid = FALSE,
    showticklabels = FALSE,
    zeroline = FALSE
  )
  
  # ===========================
  # Finalize Layout
  # ===========================
  
  if (test_type == "nostats") {
    title_text <- sprintf(
      "KEGG subnetwork for known metabolites %s abundant in\n%s relative to\n%s",
      dir_word,
      str_split(contrasts[[i]], pattern = "-")[[1]][1],
      str_split(contrasts[[i]], pattern = "-")[[1]][2]
    )
  } else {
    title_text <- sprintf(
      "KEGG subnetwork for known metabolites significantly %s abundant in\n%s relative to\n%s",
      dir_word,
      str_split(contrasts[[i]], pattern = "-")[[1]][1],
      str_split(contrasts[[i]], pattern = "-")[[1]][2]
    )
  }
  
  network <- layout(
    network,
    title = title_text,
    shapes = mydata[[paste("edge_shapes", dir, sep = ".")]][[i]],
    xaxis = mydata[[paste("axis", dir, sep = ".")]][[i]],
    yaxis = mydata[[paste("axis", dir, sep = ".")]][[i]],
    width = 800,
    height = 800,
    showlegend = TRUE,
    evaluate = TRUE
  )
  
  return(network)
}