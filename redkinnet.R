# to be documented through roxygen2
redkinnet <- function(
    kinship_matrix,
    kinship_coefficient_cutoff = "2nd_degree_cousins",
    max_links = 25
    ){
  
  # check kinship_matrix
  
  if(!is.matrix(kinship_matrix)){
    error("The kinship_matrix argument should provide a matrix.")
  }
  else if(dim(kinship_matrix)[1] != dim(kinship_matrix)[1]){
    error("The kinship matrix should have the same number of rows and columns.")
  } else if(is.null(colnames(kinship_matrix)) | is.null(rownames(kinship_matrix))){
    error("The kinship matrix should be provided with columns and row names.")
  } else {
    kin_mtx <- kinship_matrix
  }
  
  # check kinship coefficient cutoff
  
  if(is.numeric(kinship_coefficient_cutoff)){
    k_cutoff <- kinship_coefficient_cutoff
  } else if(kinship_coefficient_cutoff == "2nd_degree_cousins"){
    k_cutoff <- 0.025
  } # else if(){} # add other known kinship coefficient cutoffs (ask Pio and check in literature) 
  else {error("You must specify a valid kinship_coefficient_cutoff.")}
  cat(paste0("The selected kinship coefficient cutoff considered is ", k_cutoff))
  
  # global filter
  kin_mtx[kin_mtx <= k_cutoff] <- 0
  kin_mtx[kin_mtx > k_cutoff] <- 1
  diag(kin_mtx) <- 0
  
  # create network
  kin_net <- graph_from_adjacency_matrix(kin_mtx, mode = "undirected")
  
  # delete nodes with more than max_links
  kin_deg <- degree(kin_net)
  kin_net <- delete_vertices(kin_net, names(kin_deg)[kin_deg > max_links])
  cat(paste0("\nThe following subjects were excluded from the analysis because they have more than ", max_links, " (argument max_links) links: "))
  cat(names(kin_deg)[kin_deg > max_links])
  cat(paste0("\n", length(names(kin_deg)[kin_deg > max_links])," subjects excluded"))
  
  # get connected components and create df
  kin_comp <- components(kin_net)
  kin_data <- data.frame(comp_id = kin_comp$membership) # component id
  kin_data$name <- rownames(kin_data) # node name
  
  # component size 
  kin_data <- merge(kin_data, 
                    cbind(comp_id = unique(kin_comp$membership), 
                          comp_size = kin_comp$csize)
  )
  
  # node degree
  kin_data <- merge(kin_data, 
                    cbind(name = names(degree(kin_net)), 
                          degree = degree(kin_net))
  ) 
  
  # get isomorphisms class
  kin_data$iso_id <- rep(NA, nrow(kin_data))
  
  # assign iso class to each node
  iso3_id <- which(motifs(subgraph(kin_net, kin_data[kin_data$comp_size == 3,]$name), 3) > 0) - 1
  
  for(iso_id in iso3_id){
    names <- subgraph_isomorphisms(graph_from_isomorphism_class(3, iso_id, directed = FALSE),  
                                   subgraph(kin_net, kin_data[kin_data$comp_size == 3,]$name))
    
    kin_data[kin_data$name %in% names(unlist(names)), ]$iso_id <- iso_id
  }
  
  iso4_id <- which(motifs(subgraph(kin_net, kin_data[kin_data$comp_size == 4,]$name), 4) > 0) - 1
  
  for(iso_id in iso4_id){
    names <- subgraph_isomorphisms(graph_from_isomorphism_class(4, iso_id, directed = FALSE),  
                                   subgraph(kin_net, kin_data[kin_data$comp_size == 4,]$name))
    
    kin_data[kin_data$name %in% names(unlist(names)), ]$iso_id <- iso_id
  }
  
  iso5_id <- which(motifs(subgraph(kin_net, kin_data[kin_data$comp_size == 5,]$name), 5) > 0) - 1
  
  for(iso_id in iso5_id){
    names <- subgraph_isomorphisms(graph_from_isomorphism_class(5, iso_id, directed = FALSE),  
                                   subgraph(kin_net, kin_data[kin_data$comp_size == 5,]$name))
    
    kin_data[kin_data$name %in% names(unlist(names)), ]$iso_id <- iso_id
  }
  
  iso6_id <- which(motifs(subgraph(kin_net, kin_data[kin_data$comp_size == 6,]$name), 6) > 0) - 1
  
  for(iso_id in iso6_id){
    names <- subgraph_isomorphisms(graph_from_isomorphism_class(6, iso_id, directed = FALSE),  
                                   subgraph(kin_net, kin_data[kin_data$comp_size == 6,]$name))
    
    kin_data[kin_data$name %in% names(unlist(names)), ]$iso_id <- iso_id
  }
  
  # max path length (computed only for comp size 5 and 6)
  g <- subgraph(kin_net, kin_data[kin_data$comp_size %in% 5:6,]$name)
  max_path_lengths <- sapply(V(g), function(v) {
    max(distances(g, v = v, to = V(g))[is.finite(distances(g, v = v, to = V(g)))])
  })
  
  kin_data <- merge(kin_data, 
                    cbind(name = names(max_path_lengths), 
                          max_path_length = max_path_lengths),
                    all.x = TRUE)
  
  kin_data$degree <- as.integer(kin_data$degree)
  kin_data <- dplyr::arrange(kin_data, comp_id)
  
  # return dataset
  return(kin_data)
}
