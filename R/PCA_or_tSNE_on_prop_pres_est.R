PCA_or_tSNE_on_pro_pres_est <- function(tSNE = TRUE, PCA = FALSE, list_theta_path, var_exp_stars = NULL, display_plot = TRUE, save_plot = FALSE){
  #' Run PCA or tSNE on pixels with as coordinates, probabilities of presence of species. Plot Axis 1 vs axis 2, axis 2 vs axis 3 and dendrogram. Possibility to hide plot and save it in .png.
  #'
  #' @param tSNE boolean. compute tSNE reduction of dimensions algorithm. If TRUE, `PCA` must be FALSE, default is TRUE.
  #' @param PCA boolean. compute PCA reduction of dimensions algorithm. If TRUE, `tSNE` must be FALSE, default is FALSE.
  #' @param list_theta_path character vector.  with full path of Tiff files who contain probabilities of presence of each species to consider.
  #' @param var_exp_stars multilayer stars object. optional, only consider if `PCA` = TRUE. Plot as supplementary variables, default is NULL.
  #' @param display_plot boolean. Show plot, default is TRUE.
  #' @param save_plot boolean. Save plots in plot folder with png format.
  #' @return float matrix. with as much row as number of pixels and 3 columns for tSNE algorithm and as columns as number of axis with more than 1% of explained variance for PCA.
  #'
  #' @import Rtsne
  #' @import stars
  #' @import ggplot2
  #' @import ade4
  #' @import factoextra
  #' @export

  if (tSNE * PCA){
    print("Chose only one among tSNE and PCA")
    stop()
  }
  cell <- which(!is.na(split(read_stars(list_theta_path[1]))[[1]]), arr.ind = TRUE)[, 1:2] # keep only x & y axis
  matrix_values <- NULL
  for (i in list_theta_path) {
    # create matrix with probabilities for all species on each pixels
    matrix_values <- cbind(matrix_values, values(rast(i), na.rm = TRUE))
  }
  if (tSNE){
    tSNE_fit <- Rtsne(scale(matrix_values), dims = 3, max_iter = 5000, num_threads = 4,
                      verbose = TRUE)
    pixel_df <- as.data.frame(tSNE_fit$Y)
    dist_pixel <- dist(pixel_df)
    axis1_2 <- ggplot(data = pixel_df, aes(x = V1, y = V2)) +
      geom_point() +
      theme(legend.position = "bottom")
    axis2_3 <- ggplot(data = pixel_df, aes(x = V2, y = V3)) +
      geom_point() +
      theme(legend.position = "bottom")
    if (display_plot){
      print(axis1_2)
      print(axis2_3)
    }
    if (save_plot){
      dir.create(here("plot"), showWarnings = FALSE)
      ggsave(here("plot", "tSNE_axis1_axis2.png"), plot = axis1_2)
      ggsave(here("plot", "tSNE_axis2_axis3.png"), plot = axis2_3)
    }
  }

  if (PCA){
    if (!is.null(var_exp_stars)){
      if (length(var_exp_stars) != 1){
        var_exp_stars <- merge(var_exp_stars)
      }
      matrix_var_exp <- matrix(var_exp_stars[[1]][!is.na(var_exp_stars[[1]])], ncol = length(split(var_exp_stars)))
      pca_with_env <- PCA(cbind(matrix_values, matrix_var_exp), graph = FALSE, ncp = 10,
                          quanti.sup = dim(matrix_values)[2] + 1:length(split(var_exp_stars)))
    }else{
      pca_with_env <- PCA(matrix_values, graph = FALSE, ncp = 10)
    }
    dist_pixel <- dist(get_pca_ind(pca_with_env)$coord, method = "euclidian")
    pixel_df <- get_pca_ind(pca_with_env)$coord[, 1: which.min(fviz_eig(pca_with_env)$data[,2][fviz_eig(pca_with_env)$data[,2] > 1])]
    exp_variances <- fviz_eig(pca_with_env)
    PCA_1_2 <- fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), fill.ind = "orange", repel = TRUE,
                               col.quanti.sup = "black", invisible = "var", pointshape = 21, labelsize = 7, axes = c(1, 2),
                               col.ind = "NA", alpha = 0.8) +
      ggtitle("PCA on pixels with explanatories variables in sup") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    PCA_2_3 <- fviz_pca_biplot(pca_with_env, choix = "ind", label = c("quali", "quanti.sup"), fill.ind = "orange", repel = TRUE,
                               col.quanti.sup = "black", invisible = "var", pointshape = 21, labelsize = 7, axes = c(2, 3),
                               col.ind = "NA", alpha = 0.8) +
      ggtitle("PCA on pixels with explanatories variables in sup") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

    if (display_plot){
      print(exp_variances)
      print(paste0("Explained variances by first ", which.min(fviz_eig(pca_with_env)$data[,2][fviz_eig(pca_with_env)$data[,2] > 1]),
                   " axis : ",
                   round(sum(fviz_eig(pca_with_env)$data[,2][fviz_eig(pca_with_env)$data[,2] > 1])),
                   "%"))
      print(PCA_1_2)
      print(PCA_2_3)
    }
    if (save_plot){
      ggsave(here("plot", "PCA_explained_variances.png"), plot = exp_variances)
      ggsave(here("plot", "PCA_with_environ_sup_axis_1_2.png"), plot = PCA_1_2)
      ggsave(here("plot", "PCA_with_environ_sup_axis_1_2.png"), plot = PCA_2_3)
    }
  }

  # Plot CAH for choosen method
  if (display_plot){
    plot(hclust(dist_pixel), labels = FALSE, main = "", xlab = "", ylab = "", sub = "", ylim = "none")
  }
  if (save_plot){
    method = c(PCA, tSNE)
    names(method) <- c("PCA", "tSNE")
    png(here("plot", paste0("dendro_", names(method)[max(method)], ".png")))
    plot(hclust(dist_pixel), labels = FALSE, main = "", xlab = "", ylab = "", sub = "", ylim = "none")
    dev.off()
  }
  return(pixel_df)
}
