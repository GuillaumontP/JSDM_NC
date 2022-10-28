plot_HCA_EM_Kmeans <- function(method = "KM", pixel_df, nb_group, display_plot_3d = FALSE, display_plot_2d = TRUE, save_plot_2d = FALSE){
  #' Display plots in 2 & 3D with number of groups set and a choosen method. Possibility to hide 2D plots and save it.
  #'
  #' @param method character. clustering method to chose among c("HCA", EM, "KM"). See `details` for more explications. Default is KM.
  #' @param pixel_df dataframe. output of "dist" function in "stats" library or "PCA_or_tSNE_on_pro_pres_est" function.
  #' @param nb_group int. number of groups for clustering. Please set nb_group to high than 1.
  #' @param display_plot_3d boolean. Display new window with 3d plot, default is FALSE.
  #' @param display_plot_2d boolean. Display 2d plots, default is TRUE.
  #' @param save_plot_2d boolean. Save in plot folder as .png file all 2d plots, default is FALSE.
  #' @return int vector. Group number for each pixel.
  #' @details `method` values are HCA : Hierarchical Cluster Analysis; EM : Expectation Maximisation; KM : K-means algorithm
  #'
  #' @import EMCluster
  #' @import ggplot2
  #' @import viridis
  #' @import rgl
  #' @import stats
  #' @import factoextra
  #'
  #' @export

  if (nb_group == 1){
    print("Please select more than 1 group")
    stop()
  }
  if (method == "KM"){
    KM_class <- kmeans(pixel_df, centers = nb_group, iter.max = 20, nstart = 1000)
    pixel_group <- KM_class$cluster
  }
  if (method == "EM"){
    init_EM <- rand.EM(pixel_df, nclass = nb_group, min.n = 10)
    EM_class <- assign.class(pixel_df, emcluster(pixel_df, init_EM))
    pixel_group <- EM_class$class
  }
  if (method == "CAH"){
    dist_pixel <- dist(pixel_df)
    pixel_group <- cutree(hclust(dist_pixel), k = nb_group)
  }
  colnames(pixel_df) <- paste0("V", 1:dim(pixel_df)[2])
  plot1_2 <- ggplot(data = data.frame(pixel_df, pixel_group),
                    aes(x = V1, y = V2, color = as.factor(pixel_group))) +
    scale_color_manual(values = viridis(nb_group)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
  plot2_3 <- ggplot(data = data.frame(pixel_df, pixel_group),
                    aes(x = V2, y = V3, color = as.factor(pixel_group))) +
    scale_color_manual(values = viridis(nb_group)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
  if (display_plot_2d){
    print(plot1_2)
    print(plot2_3)
  }
  if (save_plot_2d){
    dir.create(here("plot"), showWarnings = FALSE)
    ggsave(here("plot", paste0("axis_1_2_", method, "_groups.png")), plot = plot1_2)
    ggsave(here("plot", paste0("axis_2_3_", method, "_groups.png")), plot = plot2_3)
  }
  if (display_plot_3d){
    col_3d <- viridis(nb_group)[pixel_group]
    pixel_df <- data.frame(pixel_df)
    par3d(windowRect = c(20, 30, 800, 800))
    plot3d(pixel_df$V1, pixel_df$V2, pixel_df$V3, col = col_3d, xlab = "x", ylab = "y", zlab = "z", size = 5, pch = 16)
    legend3d("bottomright", legend = paste0("Group_", 1:length(unique(pixel_group))), pch = 16, col = viridis(nb_group),
             cex = 2, inset = c(0.02))
  }
  return(pixel_group)
}
