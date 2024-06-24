# Load necessary libraries
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(ape)
library(tidyverse)

# Function to read data and plot the dendrogram
main <- function() {
  # Define file paths
  smiles_combinations_file <- "~/Git/mia-playground/skeletons_for_mces.csv"
  # smiles_combinations_file <- "~/Git/mia-playground/skeletons_distances_mhfp.csv"
  skeleton_names_file <- "~/Git/mia-playground/src/miaplayground/data/skeletons.csv"
  skeleton_names_new_file <- "~/Git/mia-playground/src/miaplayground/data/skeletons_new.csv"
  distances_mces_file <- "~/Git/mia-playground/skeletons_distances_mces.csv"
  distances_mhfp_file <- "~/Git/mia-playground/skeletons_distances_mhfp.csv"
  
  # Load data
  smiles_combinations <- read_csv(smiles_combinations_file,
                                  col_names = c("index", "smiles_1", "smiles_2"))
  skeleton_names <- read_csv(skeleton_names_file) |>
    # rename(name = ATTRIBUTE_Skeletons, smarts = SMARTS,smiles = Skeleton_SMILES) |>
    rename(name = Skeletons, smiles = SMILES) |>
    # distinct(name, smarts,smiles) |>
    distinct(name, smiles)
  distances <- read_csv(distances_mces_file,
                        col_names = c("index", "time", "distance", "type"))
  
  # Merge data based on index
  merged_data <- left_join(smiles_combinations, distances, by = "index") |>
    # merged_data <- smiles_combinations |>
    #   rename(distance=X4) |>
    left_join(skeleton_names, by = c("smiles_1" = "smiles")) |>
    left_join(skeleton_names, by = c("smiles_2" = "smiles")) |>
    distinct()
  # merged_data <- smiles_combinations |>
  #   rename(distance=X4) |>
  #   left_join(skeleton_names, by = c("smiles_1" = "smarts")) |>
  #   left_join(skeleton_names, by = c("smiles_2" = "smarts")) |>
  #   distinct()
  #   # mutate(distance=as.integer(-distance))
  # merged_data <- smiles_combinations |>
  #   rename(distance=X5) |>
  #   left_join(skeleton_names, by = c("smiles_1" = "smiles")) |>
  #   left_join(skeleton_names, by = c("smiles_2" = "smiles")) |>
  #   distinct()
  
  ## TODO
  merged_data |>
    write_csv("~/Git/mia-playground/temp.csv")
  
  # Create a distance matrix
  mat <- merged_data |>
    distinct(skeleton = name.x, skeleton_2 = name.y, distance) |>
    pivot_wider(names_from = skeleton_2, values_from = distance)
  name_1 <- mat$skeleton[1]
  name_2 <- mat$skeleton[2]
  name_last <- colnames(mat)[length(colnames(mat))]
  mat <- mat |>
    mutate(!!as.name(name_1) := NA, ) |>
    relocate(!!as.name(name_1), .before = !!as.name(name_2)) |>
    add_row() |>
    mutate(skeleton = ifelse(
      test = is.na(skeleton),
      yes = name_last,
      no = skeleton
    )) |>
    column_to_rownames("skeleton") |>
    t()
  
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  
  # mat[mat = 0] <- 1
  # mat = as.dist(mat,upper = TRUE)
  
  # Convert distance matrix to a dist object
  distance_matrix <- as.dist(1 - mat)
  distance_matrix <- as.dist(1 / mat)
  distance_matrix <- mat
  
  pheatmap::pheatmap(
    mat = mat,
    display_numbers = TRUE,
    fontsize_number = 6,
    color = khroma::color(palette = "batlowW", reverse = TRUE)(length(unique(
      distances$distance
    )))
  )
  
  
  distance_matrix[is.na(distance_matrix)] <- 0
  
  # Build the tree
  hc <- hclust(as.dist(distance_matrix))
  plot(hc)
  clus.den <- as.dendrogram(object = hc)
  
  # Redefine skeletons with <= 3 distance
  obj <- cut(clus.den, h = 3)
  plot(obj$lower[[1]])
  test <- obj$lower[lapply(obj$lower, length) > 1]
  labels <- lapply(test, labels)
  labels_pasted <- lapply(labels, paste, collapse = "+")
  skeleton_names_new <- lapply(1:length(labels), function(x) {
    skeleton_names_new <- skeleton_names |>
      mutate(name_new = ifelse(
        test = name %in% labels[[x]],
        yes = labels_pasted[[x]],
        no = name
      )) |>
      filter(name_new != name)
  }) |> bind_rows() |>
    right_join(skeleton_names) |>
    mutate(name_new = ifelse(
      test = is.na(name_new),
      yes = name,
      no = name_new
    ))
  
  skeleton_names_new |>
    write_csv(file = skeleton_names_new_file)
  
  #
  #   tree <- hc |>
  #     ggtree::ggtree()
  #   images <- paste0(
  #     "https://www.simolecule.com/cdkdepict/depict/bot/svg?smi=",
  #     unique(merged_data$smiles_1)
  #   )
  #
  #   tree_styled <- tree +
  #     ggtree::geom_tiplab(
  #       mapping = ggtree::aes(
  #         image = images,
  #         size = ggplot2::rel(0.025),
  #         offset = ggplot2::rel(0.05),
  #         geom = "image"
  #       )
  #     )
  #
  #   # Plot the dendrogram
  #   plot(
  #     as.phylo(hc),
  #     tip.label = rownames(mat),
  #     cex = 0.6,
  #     no.margin = TRUE
  #   )
}

# Run the main function
main()
