## Functions for srIBI by Sizun Jiang
plot.umap = function(x, labels,
                     main="UMAP visualization",
                     colors=c("#ff7f00", "#e377c2", "#17becf"),
                     pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=1) {
  
  layout = x
  if (class(x)=="umap") {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

# From:
# https://blog.datascienceheroes.com/playing-with-dimensions-from-clustering-pca-t-sne-to-carl-sagan/
plot_cluster=function(data, var_cluster, palette)
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}

plot_cluster_multicolor=function(data, var_cluster, mycolor)
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_manual(values = mycolor) 
}

plot_gradient=function(data, var_cluster)
{
  data[[var_cluster]] = as.numeric(data[[var_cluster]])
  mid = mean(data[[var_cluster]])
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.25) +
    # guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    #ggtitle(as.character(var_cluster)) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.key.size = unit(3,"line")) + 
    scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                          high="red", space ="Lab") 
    #scale_color_gradient2(palette) 
}

get_cluster_mask = function(celldf, xdim, ydim, cluster_col_name, clusterID){
  # CellDF should have a column called ID, and column containing cluster
  # Returns a mask with same dimensions
    overlay = matrix(0,xdim, ydim)
    w = which(celldf[[cluster_col_name]] == clusterID)
  for(i in w){
    overlay[(celldf$x_start[i]*step_xy):(celldf$x_start[i]*step_xy + celldf$x_size[i] - 1),
            (celldf$y_start[i]*step_xy):(celldf$y_start[i]*step_xy+ celldf$y_size[i] - 1)] = 1
  }
    return(overlay)
}

find_knn=function(cluster_list, cluster_k, NNs, random_num)
  {
  # Set color palette
  mycolors = c(brewer.pal(name = "Set1", 9), brewer.pal(name="Set3", 12))
  
  # Version that removes cluster_k
  mycolors = mycolors[-cluster_k]
  
  # Add a randomized control with normal distribution
  set.seed(123456)
  random_x = runif(n = random_num, min = min(cluster_list[[cluster_k]]$x_start), max = max(cluster_list[[cluster_k]]$x_start))/1000
  random_y = runif(n = random_num, min = min(cluster_list[[cluster_k]]$y_start), max = max(cluster_list[[cluster_k]]$y_start))/1000
  random_z = runif(n = random_num, min = min(cluster_list[[cluster_k]]$z_start), max = max(cluster_list[[cluster_k]]$z_start))/1000
  
  random_pp3 = pp3(random_x,
                   random_y,
                   random_z,
                   box3(c(0,1)))
  # Requires a list per cluster, containing a df of the voxels in the cluster, make sure the cluster info and name are included Eg
  # Cycle through cellID
  knn_df_cells = do.call(rbind, lapply(unique(cluster_list[[cluster_k]]$cellID), function(cell){
    cluster_k_cell = (cluster_list[[cluster_k]])[which(cluster_list[[cluster_k]]$cellID == cell),]
    cluster_k_pp3 = pp3(cluster_k_cell$x_start/1000, 
                        cluster_k_cell$y_start/1000, 
                        cluster_k_cell$z_start/1000, box3(c(0,1)))
    knn_list = lapply(cluster_list[-which(names(cluster_list) == cluster_k)], function(df){
      # Get X,Y,Z coordinates and turn them into a box
      cluster_test_pp3 = pp3(df$x_start/1000, 
                             df$y_start/1000, 
                             df$z_start/1000, box3(c(0,1)))
      NN_df = nncross(cluster_k_pp3, cluster_test_pp3, k=1:NNs, what = c("dist", "which"))
      return(NN_df)
    })
    knn_list$random = nncross(cluster_k_pp3, random_pp3, k=1:NNs, what = c("dist", "which"))
    knn_df = do.call(rbind, lapply(names(knn_list), function(i){
      df = knn_list[[i]]
      new_df = data.frame(unlist(df[,1:NNs]), rep(i, length(unlist(df[,1:NNs]))))
      # rownames(new_df) = unlist(df[,(NNs + 1):(NNs*2)])
      colnames(new_df) = c("Dist", "Cluster")
      return(new_df)
    }))
  }))
  
  # wilcoxon_results = unlist(lapply(names(knn_list), function(i){
  #   rand = knn_list$random
  #   df = knn_list[[i]]
  #   wil = wilcox.test(unlist(df[,1:NNs]), unlist(rand[,1:NNs]))
  #   return(wil$p.value)
  # }))
  # names(wilcoxon_results) = names(knn_list)
  # # ttest
  # ttest_results = unlist(lapply(names(knn_list), function(i){
  #   rand = knn_list$random
  #   df = knn_list[[i]]
  #   wil = t.test(unlist(df[,1:NNs]), unlist(rand[,1:NNs]))
  #   return(wil$p.value)
  # }))
  # names(ttest_results) = names(knn_list)
  # ks test for distribution
  # ks_results = unlist(lapply(names(knn_list), function(i){
  #   rand = knn_list$random
  #   df = knn_list[[i]]
  #   wil = ks.test(unlist(df[,1:NNs]), unlist(rand[,1:NNs]))
  #   return(wil$p.value)
  # }))
  # names(ks_results) = names(knn_list)
  
  # Plot density plot
  p = ggplot(knn_df_cells, aes(x = Dist, fill = Cluster)) +
    #geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 100)
    geom_density(alpha = 0.65) + scale_fill_manual(values = mycolors) + 
    ggtitle(paste0("Distance of ", NNs, " Nearest Neighbors around Cluster ", cluster_k)) +
    xlab("Distance (*1e9 Pixels)") +
    ylab("Density")
  # Plot violin plot
  p2 = ggplot(knn_df_cells, aes(x = Cluster, y = Dist, fill = Cluster)) +
    #geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 100)
    geom_violin(alpha = 0.65, trim = TRUE) + scale_fill_manual(values = mycolors) +
    geom_boxplot(width = 0.1, outlier.shape = NA) + ggtitle(paste0("Distance of ", NNs, " Nearest Neighbors around Cluster ", cluster_k)) +
    xlab("Cluster Number") +
    ylab(paste("Distance (*1e9 Pixels)"))
  ret= list(knn_df_cells, p, p2)
  names(ret) = c("knn_df", "plotDensity", "plotViolin")
  return(ret)
}

plot_marker_dist=function(cluster_list, channel_name, whisker_cutoff){
  mycolors = c(brewer.pal(name = "Set1", 9), brewer.pal(name="Set3", 12))
  channel_df = do.call(rbind, lapply(names(cluster_list), function(i){
    df = cluster_list[[i]]
    z = (df[[channel_name]]-mean(df[[channel_name]]))/(sd(df[[channel_name]])+0.0001)
    z = z[is.finite(z)]
    new_df = data.frame(z, rep(i, length(df[[channel_name]])))
    names(new_df) = c("channel", "cluster")
    return(new_df)
  }))
  p0 = ggplot(channel_df, aes(x = cluster, y = channel, fill = cluster)) +
    #geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 100)
    geom_violin(alpha = 0.65, trim = TRUE) + scale_fill_manual(values = mycolors) +
    geom_boxplot(width = 0.1, outlier.shape = NA) + ggtitle(paste0("Count distribution of ", channel_name," per cluster")) +
    xlab("Cluster Number") +
    ylab(paste("Z-score distribution of counts"))
    # Test limit
  if (missing(whisker_cutoff)){
    return(p0)
  } else{
  ylim1 = boxplot.stats(channel_df$channel)$stats[c(1, 5)]
  p1 = p0 + coord_cartesian(ylim = ylim1*whisker_cutoff)
  return(p1)
  }
}

# Extract features based on a voxel size
extract_features = function(cell_list, x, y, z, step_xy, step_z, discard){
  # From a cell_list (list of tiffs), apply a window size of:
  # x * y * z pixels, stepping step_xy and z pixel sized steps
  # discard/ignore pixels around each border with discard
  ext_features = lapply(names(cell_list), function(tiff){
  # Shave
  tiff_small = cell_list[[tiff]][(discard+1):(dim(cell_list[[tiff]])[1]-discard),
                                      (discard+1):(dim(cell_list[[tiff]])[2]-discard),
                                      1:dim(cell_list[[tiff]])[3]]
  x_start = 1:floor((dim(tiff_small)[1]-x)/step_xy)
  y_start = 1:floor((dim(tiff_small)[2]-y)/step_xy)
  z_start = 1:floor((dim(tiff_small)[3]-z)/step_z)
  x_count = numeric(length(x_start) * length(y_start) * length(z_start))
  y_count = numeric(length(x_count))
  z_count = numeric(length(x_count))
  ave_vec = numeric(length(x_start))
  id_vec = 1:(length(x_count))
  id = 1
  # Step x, y, z
  for(i in x_start) {
    for(j in y_start) {
      for(k in z_start) {
        new_i = (i-1)*step_xy+1
        new_j = (j-1)*step_xy+1
        new_k = (k-1)*step_z+1
        ave_feat = mean(tiff_small[new_i:(new_i+x-1),
                                   new_j:(new_j+y-1),
                                   new_k:(new_k+z-1)])
        ## Vectorized the code to be much faster! (no incrementally increases in size within the loop!)
        ave_vec[id] = ave_feat
        x_count[id] = i
        y_count[id] = j
        z_count[id] = k
        id = id+1
      }
    }
  }
  df = data.frame(id_vec, ave_vec, x_count, y_count, z_count,
                  rep(x, length(x_count)), rep(y, length(y_count)), rep(z, length(z_count)))
  colnames(df) = c("ID", tiff, "x_start", "y_start", "z_start", "x_size", "y_size", "z_size")
  return(df)
})
names(ext_features) = names(cell_list)
return(ext_features)
}

# Sample cells and return a df of cells with a cellID
sample_cells = function(cells_feat, sample_num, cellID){
  subset_cells = cells_feat[sample(cells_feat$ID, sample_num, replace = F),]
  subset_cells$cellID = cellID
  return(subset_cells)
}

# Permutation test
# subset_cells has a DF containing xyz positions, cellID
perm_AB = function(subset_cells, namA, namB, radius, n, mc.cores){
  cluster_list = lapply(levels(subset_cells$h_cluster), function(x){
    w = which(subset_cells$h_cluster == x)
    return(subset_cells[w,])
  })
  
  num_AB = unlist(lapply(unique(subset_cells$cellID), function(cell){
    all_cluster_cell = do.call(rbind, lapply(cluster_list, function(df){
      df_cell = df[which(df$cellID == cell),]
    }))
    # Takes in a DF consisting of all clusters with a single cell ID, together with xyz start positions
    # Make a 3D box
    cluster_A = all_cluster_cell[which(all_cluster_cell$h_cluster == namA),]
    cluster_B = all_cluster_cell[which(all_cluster_cell$h_cluster == namB),]
    cluster_A_pp3 = pp3(cluster_A$x_start/1000, 
                        cluster_A$y_start/1000, 
                        cluster_A$z_start/1000, box3(c(0,1)))
    cluster_B_pp3 = pp3(cluster_B$x_start/1000, 
                        cluster_B$y_start/1000, 
                        cluster_B$z_start/1000, box3(c(0,1)))
    # Get distance matrix between A and B
    dist_matrix = crossdist.pp3(cluster_A_pp3, cluster_B_pp3, periodic = F, squared = T)
    num_in_radius = sum(dist_matrix < ((radius/1000)^2))
    # Number of AB within radius
    num_AB = num_in_radius/dim(cluster_A)[1]
  }))
  names(num_AB) = unique(subset_cells$cellID)
  
  # num_AB_perm = vector(,n)
  # shuffle clusters and try for n times
  num_AB_perm = lapply(unique(subset_cells$cellID), function(cell){
    all_cluster_cell = do.call(rbind, lapply(cluster_list, function(df){
      df_cell = df[which(df$cellID == cell),]
    }))
    num_AB_perm_cell = unlist(mclapply(1:n, function(i){
      rand_h_cluster = sample(all_cluster_cell$h_cluster)
      rand_A = all_cluster_cell[which(rand_h_cluster == namA),]
      rand_B = all_cluster_cell[which(rand_h_cluster == namB),]
      rand_A_pp3 = pp3(rand_A$x_start/1000, 
                       rand_A$y_start/1000, 
                       rand_A$z_start/1000, box3(c(0,1)))
      rand_B_pp3 = pp3(rand_B$x_start/1000, 
                       rand_B$y_start/1000, 
                       rand_B$z_start/1000, box3(c(0,1)))
      dist_matrix = crossdist.pp3(rand_A_pp3, rand_B_pp3, periodic = F, squared = T)
      num_in_radius = sum(dist_matrix < ((radius/1000)^2))
      return(num_in_radius/dim(rand_A)[1])
    }, mc.cores = mc.cores))
    return(num_AB_perm_cell)
  })
  names(num_AB_perm) = unique(subset_cells$cellID)
  
  plots = lapply(unique(subset_cells$cellID), function(cell){
    p = ggplot(as.data.frame(num_AB_perm[[cell]]), aes(x = num_AB_perm[[cell]])) + geom_density() + 
      geom_vline(xintercept = num_AB[[cell]], linetype = "dotted", color = "red", size = 1) +
      ggtitle(paste0(cell, " Interactions between Clusters ", namA, " and ", namB, ". ", n, " Permutations per FOV")) +
      xlab(paste0("Interaction Frequency within Radius ", radius)) +
      ylab(paste("Frequency"))
    return(p)
  })
  names(plots) = unique(subset_cells$cellID)
  # Adjusted Pvalues (There is a +1 so that pval is never 0!)
  adj.pval = p.adjust(unlist(lapply(unique(subset_cells$cellID), function(cell){
    if(mean(num_AB_perm[[cell]]) > num_AB[[cell]]){
      pval = (sum(num_AB_perm[[cell]] < num_AB[[cell]] )+ 1)/(length(num_AB_perm[[cell]]) + 1)
    } else{
      pval = (sum(num_AB_perm[[cell]] > num_AB[[cell]] )+ 1)/(length(num_AB_perm[[cell]]) + 1)
    }
    return(pval)
  })), method = "BH")
  
  # log2 fold Enrichment of interaction over permutation mean
  log2FC = log2(mean(num_AB)/mean(unlist(num_AB_perm)))
  
  ret = list(plots, adj.pval, num_AB, num_AB_perm, log2FC)
  names(ret) = c("plots", "adj.pval", "AB", "AB_perm", "log2FC")
  return(ret)
}