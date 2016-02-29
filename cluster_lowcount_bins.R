#This function will cluster bins from `X_bin` when they have low count.
#It will cluster according to its closest neighboor.
#
#cluster_lowcount_bins: It clusters 1 dimensional data.
#cluster_lowcount_bins_2D: It clusters 2 dimensional data.
#
#THIS IS REMOVED:
#cluster_lowcount_bins_return_break_points: It clusters 1 dimensional data. It returns the breaking points for clustering the data
#
#
#
cluster_lowcount_bins <- function(X_bin, n_thres=10, do_edges_only=T){
#This function will cluster data points from the edge in case
  a=aggregate(data.frame(count=X_bin), by=list(X_bin=X_bin),length);
  a <- a[order(a$X_bin),];
  
  ind_lt_thres_ordered <- order(a$count)[a$count[order(a$count)]<n_thres];
  
  n_bins <- length(a$X_bin);
  
  if(length(ind_lt_thres_ordered)==0){
    return(X_bin)
  }else{
    if(ind_lt_thres_ordered[1]==1){
      X_bin[X_bin==a$X_bin[1]] = a$X_bin[2];
      cluster_lowcount_bins(X_bin, n_thres, do_edges_only);
    }else if (ind_lt_thres_ordered[1]==n_bins){
      X_bin[X_bin==a$X_bin[n_bins]] = a$X_bin[n_bins-1];
      cluster_lowcount_bins(X_bin, n_thres, do_edges_only);
    }else if(do_edges_only){
      stop("This script was set to cluster only at the edges!");
    }else{
      #It will cluster with neighboor that contains least counts.
      if(a$count[ind_lt_thres_ordered[1]-1] < a$count[ind_lt_thres_ordered[1]+1]){
        choose_neighboor=ind_lt_thres_ordered[1]-1;
      }else{
        choose_neighboor=ind_lt_thres_ordered[1]+1;
      }
      X_bin[X_bin==a$X_bin[ind_lt_thres_ordered[1]]] = a$X_bin[choose_neighboor];
      cluster_lowcount_bins(X_bin, n_thres, do_edges_only);
    }
  }
}

cluster_lowcount_bins_with_label <- function(X_bin, X_label, n_thres=10, do_edges_only=T){
  
  X_label_unique <- unique(X_label);
  n_labels <- length(X_label_unique);
  
  X_bin_update <- X_bin;
  for(l in 1:n_labels){
    
    cur_label = X_label_unique[l];
    ind_cur_label <- X_label==cur_label;
    
    X_bin_update[ind_cur_label] <- cluster_lowcount_bins(X_bin = X_bin[ind_cur_label],n_thres = n_thres, do_edges_only = do_edges_only);
  }
  return(X_bin_update)
}
  
cluster_lowcount_bins_2D <- function(X_bin, Y_bin, n_thres=10, do_edges_only=T,
                                     X_bin_size, Y_bin_size,nloop=0){
  #This function will cluster data points from the edge in case
  a=aggregate(data.frame(count=X_bin), by=list(X_bin=X_bin, Y_bin=Y_bin),length);
  a <- a[order(a$X_bin, a$Y_bin),]; #The aggregated data is ordered according to X_bin and Y_bin;
  
  ind_lt_thres_ordered <- order(a$count)[a$count[order(a$count)]<n_thres];
  
  #It returns the final data when all bins have at least `n_thres` counts.
  if(length(ind_lt_thres_ordered)==0){
    return( data.frame(X_bin=X_bin, Y_bin=Y_bin, nloop=nloop+1));
  }
  #
  #n_bins <- length(a$X_bin);
  cur_lowcount_X_bin = a$X_bin[ind_lt_thres_ordered[1]];
  cur_lowcount_Y_bin = a$Y_bin[ind_lt_thres_ordered[1]];
  if(cur_lowcount_X_bin!=a$X_bin[1] & 
       cur_lowcount_X_bin!=max(a$X_bin) & 
       cur_lowcount_Y_bin!=a$Y_bin[1] &
       cur_lowcount_Y_bin!=max(a$Y_bin)){
    if(do_edges_only){
      stop("This script was set to cluster only at the edges!");
    }
  }
  
  #Choose only points with closest distance to X and Y to aggregate;
  dist_to_lowest_X_bin <- ((a$X_bin - cur_lowcount_X_bin)/X_bin_size)^2 + ((a$Y_bin - cur_lowcount_Y_bin)/Y_bin_size)^2;
  min_dist <- min(dist_to_lowest_X_bin[dist_to_lowest_X_bin!=0]);
  ind_min_dist = which(dist_to_lowest_X_bin==min_dist);
  
  #Get best neighboor;
  ind_best_neighboor= ind_min_dist[order(a$count[ind_min_dist])][1];
  X_bin_neighboor= a$X_bin[ind_best_neighboor];
  Y_bin_neighboor= a$Y_bin[ind_best_neighboor];
  
  
  X_bin[X_bin==cur_lowcount_X_bin & Y_bin==cur_lowcount_Y_bin] = X_bin_neighboor;
  Y_bin[X_bin==cur_lowcount_X_bin & Y_bin==cur_lowcount_Y_bin] = Y_bin_neighboor;
  
  cluster_lowcount_bins_2D(X_bin, Y_bin, n_thres, do_edges_only,
                           X_bin_size, Y_bin_size, nloop=nloop+1);
  
}