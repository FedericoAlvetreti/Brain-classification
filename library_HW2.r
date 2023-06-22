# Libraries ---------------------------------------------------------------
library(igraph)
library(readr)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(e1071)
library(tidyverse)
library(MNITemplate)
library(rgl)
library(vars)
library(caret)
library(graphkernels)
library(Hmisc)
library(sgof)
library(fclust)
library(cluster)
library(stepPlr)
library(tseries)
library(stats)
# Preprocessing -----------------------------------------------------------

extract_roi <- function(df) {
  
  # Get roi's columns names
  data_cols <- names(df[1,])[4:(ncol(df[1,])-2)]
  
  # Select only roi's columns and converting to tibble
  df_tibble <- as_tibble(df[1,][,4:(ncol(df[1,])-2)])
  
  # Convert columns into rows
  df_long <- df_tibble %>%
    pivot_longer(cols = all_of(data_cols), 
                 names_to = "column_name", 
                 values_to = "value") %>%
    separate(column_name, 
             into = c("time_instant", "ROI_name"), sep = "_")
  
  # Convert rows into columns
  df_wide <- df_long %>%
    pivot_wider(names_from = ROI_name, values_from = value)
  
  return(df_wide[,2:ncol(df_wide)])
}

find_null_series <- function(df) {
  
  # Method specifically meant for the train-set
  
  # Initialize a vector to store the rows with at least one zero signal
  
  supp <- c()
  
  for (i in 1:nrow(df)) {
    
    data <- extract_roi(df[i,])
    if (ncol(data %>% select(where(~ all(. == 0)))) > 0) {
      supp <- append(supp, i)
    }
  }
  return(df[-supp,])
}

correct_null_series <- function(df) {
  
  # Method specifically meant for the test-set
  for (i in 1:nrow(df)) {
    
    # Detect whether there's a null signal in the row 
    
    if (sum(df[i,]==0) > 115) {
      
      for (j in 1:ncol(df)-114) {
        
        # Find all the null signals in the row and correct them
        
        if (all(df[i,j:(j+114)] == c(rep(0,115)))) {
          
          df[i,j:(j+114)] <- colMeans(data.frame(df[,(j+114)]))
          
        }
      }
    }
  }
  return(df)
}

min_max_rescaling <- function(array) {
  
  # This method rescales the data to the range [0,1]
  
  return((array-min(array))/(max(array)-min(array)))
}

differentiate <- function(series) {
  
  # Apply first order differencing to the series
  
  for (i in 2:length(series)) {
    series[i] <- series[i] - series[i-1]
  }
  return(series)
}


make_stationary <- function(series) {
  # Calculate the first difference
  diff_series <- differentiate(series)
  
  # Check if differencing is sufficient
  adf_test <- adf.test(diff_series)
  p_value <- adf_test$p.value
  
  if (p_value < 0.05) {
    
    # Differencing is sufficient, return the differenced series
    diff_series[1] <- 0
    return(diff_series)
    
  } else {
    
    # Differencing is not sufficient, go recursive
    return(make_stationary(diff_series))
  }
}


fft_denoiser <- function(x, n_components=0.005) {
  n <- length(x)
  
  # Compute the FFT
  fft_result <- fft(x)
  
  # Compute power spectral density (PSD)
  # Squared magnitude of each FFT coefficient
  PSD <- fft_result * Conj(fft_result) / n
  
  # Keep high frequencies
  mask <- Re(PSD) > n_components
  fft_result <- mask * fft_result
  
  # Inverse Fourier transform
  clean_data <- Re(fft(fft_result, inverse = TRUE))
  
  return(clean_data)
}


# Feature engineering -----------------------------------------------------

best_k_distro <- function(df){
  get_distances <- function(ts){
    library(Hmisc)
    P_mat <- rcorr(as.matrix(ts),type = "spearman")$P
    n <- dim(P_mat)[1]
    P_mat[is.na(P_mat)] <- runif(sum(is.na(P_mat)))
    v <- c(P_mat)
    x <- data.frame(row.names = c(1:length(v)))
    x$P <- v
    x$idx <- c(1:length(v))
    x <- x[order(x$P,na.last = T),]
    x$P <- sort(v) * length(v) / c(1:length(v))
    x <- x[order(x$idx,na.last = T),]
    v <- x$P
    P_mat <- matrix(v,n,n)
    diag(P_mat) <- rep(0,n)
    return(P_mat)
  }
  get_best_k <- function(ts,k_max=20,p=2){
    library(cluster)
    W <- get_distances(ts)
    A <- 1-W
    D <- diag(colSums(A,na.rm = T))
    W[is.na(W)] <- 0
    L <- D-W
    sl_score <- c()
    slope_score <- c()
    for(k in c(2:k_max)){
      U <- as.matrix(eigen(L)$vectors)[,c(116:(116-k+1))]
      km <- kmeans(U, centers = k, nstart=25)
      ss <- silhouette(km$cluster, dist(U))
      sl_score <- c(sl_score,mean(ss[, 3]))
    }
    for(i in c(1:(length(sl_score)-1))){
      slope_score <- c(slope_score,((sl_score[i]-sl_score[i+1])*sl_score[i]^p))
    }
    return(which.max(slope_score)+1)
  }
  extract_roi <- function(df) {
    library(tidyverse)
    df_tibble <- as_tibble(df)
    df_long <- df_tibble %>%
      pivot_longer(cols = names(df), 
                   names_to = "column_name", 
                   values_to = "value") %>%
      separate(column_name, 
               into = c("time_instant", "ROI_name"), sep = "_")
    df_wide <- df_long %>%
      pivot_wider(names_from = ROI_name, values_from = value)
    
    return(data.frame(df_wide[-1]))
    
  }
  return_features <- function(df_row){
    
    clean_row <- df_row[,!names(df_row) %in% c("id","male","female","age", "y")]
    
    fRMI <- extract_roi(clean_row)
    
    return(get_best_k(fRMI))
    
  }
  result <- foreach(i = c(1:dim(df)[1]), .combine = rbind) %dopar% 
    return_features(df[i,])
  return(result)
}

get_patients_features <- function(df, k_opt){
  
  # Get all feature engineering stuff
  get_features <- function(ts, k_opt){
    library(igraph)
    library(graphkernels)
    library(fclust)
    library(tseries)
    
    # Get dissimilarity matrix obtained by the p value corrected 
    # for the false discovery rate FDR corresponding to the Spearman correlation
    get_distances <- function(ts){
      library(Hmisc)
      
      P_mat <- rcorr(as.matrix(ts),type = "spearman")$P
      n <- dim(P_mat)[1]
      diag(P_mat) <- 0
      v <- c(P_mat)
      x <- data.frame(row.names = c(1:length(v)))
      x$P <- v
      x$idx <- c(1:length(v))
      x <- x[order(x$P,na.last = T),]
      x$P <- c(rep(0,116),sort(v[117:length(v)]) * (length(v)-116) / c(1:(length(v)-116)))
      x <- x[order(x$idx,na.last = T),]
      v <- x$P
      P_mat <- matrix(v,n,n)
      
      
      return(P_mat)
    }
    
    # Get the optimal number of cluster for the patient
    get_best_k <- function(ts,k_max=20,p=2){
      library(cluster)
      W <- get_distances(ts)
      A <- 1-W
      D <- diag(colSums(A,na.rm = T))
      W[is.na(W)] <- 0
      L<-D-W
      sl_score <- c()
      slope_score <- c()
      for(k in c(2:k_max)){
        U <- as.matrix(eigen(L)$vectors)[,c(116:(116-k+1))]
        km <- kmeans(U, centers = k, nstart=25)
        ss <- silhouette(km$cluster, dist(U))
        sl_score <- c(sl_score,mean(ss[, 3]))
      }
      for(i in c(1:(length(sl_score)-1))){
        slope_score <- c(slope_score,((sl_score[i]-sl_score[i+1])*sl_score[i]^p))
      }
      return(which.max(slope_score)+1)
    }
    
    # Get Eigenvector centrality (how well is connected each node to other important nodes)
    EVC <- function(adj_mat){
      E <- eigen(adj_mat)
      idx <- which.max(abs(E$values))
      vec <- Re(E$vectors[,idx])
      return(vec)
    }
    
    # Get Shannon entropy (how dubious an assignment of an ROI into a cluster was compared with the others)
    entropy <- function(clu_mat){
      features <- c()
      for( i in c(1:dim(clu_mat)[1])){
        features <- c(features,(-sum(clu_mat[i,] * log(clu_mat[i,]))))
      }
      return(features)
    }
    
    # Get modularity (how well separated the ROIs were in terms of clustering)
    modularity <- function(m = 116^2, n = 116, adj_mat, degree_vec, clu){
      result <- matrix(0,116,116)
      for ( i in c(1:n)){
        for ( j in c(1:n)){
          result[i,j] <- (adj_mat[i,j]-degree_vec[i]*degree_vec[j]/(2*m))*(clu[i]==clu[j])
        }
      }
      result[is.na(result)] <- 0
      result <- sum(result)/(2*m)
      return(result)
    }
    
    # Compute the dissimilarity matrix
    W <- get_distances(ts)
    
    # Compute the adjacency matrix by calculating 1 – the p value corresponding to the Spearman correlation 
    A <- 1-W
    
    # Degrees vector 
    degrees <- colSums(A,na.rm = T)
    
    # Degree matrix 
    D <- diag(degrees)
    
    # Compute Laplacian matrix
    L<-D-W
    
    # Build U as the matrix of eigenvectors
    U <- as.data.frame((eigen(L)$vectors)[,c(116:(116-k_opt+1))])
    
    # Fuzzy clustering 
    C <- FKM(X = U, k = k_opt ,m = 1.4)$U
    
    # Get "hard clusters"
    w <- max.col(C)
    
    
    K <- get_best_k(ts) # 1 feature
    mod <- modularity(adj_mat = A, degree_vec = degrees, clu = w) # 1 feature
    entr <- entropy(C) # 116 features
    evc <- EVC(A) # 116 features
    
    # Set a threshold to keep only relevant edges
    treshold <- 1-(0.05/(0.5*116^2)) # 1-alpha corrected with Bonferroni 
    A[A < treshold] <- 0
    
    graph <- graph.adjacency(A, mode="undirected", weighted=TRUE)
    
    # Graph features
    G <- list(graph)
    graphlet <- CalculateConnectedGraphletKernel(G,3) # 1 feature 
    ass_coeff <- assortativity_degree(graph, directed = T) # 1 feature
    mean_distance <- mean_distance(graph,directed = FALSE) # 1 feature
    str <- quantile(strength(graph)) # 5 features
    clust_coeff <- transitivity(graph,type = 'weighted') # 116 features
    
    
    return(c(K,
             mod,
             graphlet,
             ass_coeff,
             mean_distance,
             str,
             entr,
             evc,
             clust_coeff))
  }
  
  # Extract a N x M dataframe where N='time instants', M = 'number of ROI'
  extract_roi <- function(df) {
    library(tidyverse)
    
    df_tibble <- as_tibble(df)
    df_long <- df_tibble %>%
      pivot_longer(cols = names(df), 
                   names_to = "column_name", 
                   values_to = "value") %>%
      separate(column_name, 
               into = c("time_instant", "ROI_name"), sep = "_")
    df_wide <- df_long %>%
      pivot_wider(names_from = ROI_name, values_from = value)
    
    return(data.frame(df_wide[-1]))
    
  }
  
  # Combine the initial features with the fRMI extracted features
  return_features <- function(df_row, k_opt){
    
    clean_row <- df_row[,!names(df_row) %in% c("id","male","female","age", "y")]
    
    fRMI <- extract_roi(clean_row)
    
    extracted_features <- get_features(fRMI,k_opt)
    
    initial_features <- c(df_row$male,df_row$female, df_row$age)
    
    features <- c(initial_features,extracted_features)
    features <- unlist(features)
    
    return(features)
    
  }
  
  # Parallelize the feature extraction for all patients in df
  result <- foreach(i = c(1:dim(df)[1]), .combine = rbind) %dopar% 
    return_features(df[i,],k_opt)
  
  # Formatting
  result <- as.data.frame(result)
  row.names(result) <- paste("patient",c(1:dim(df)[1]),sep=" ")
  names(result) <- c("male","female","age",paste("X",c(1:358),sep=""))
  
  return(result)
}

# LOCO --------------------------------------------------------------------

# New operator to assign multiple variables in the same line of code
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}
 
# Loss function based on cross entropy loss
accuracy <- function(model, df) {
  preds <- as.numeric(unlist(predict(model, df)))
  y <- as.numeric(unlist(df["y"]))
  return(c(abs(preds - y)))
}

# actual loco 
LOCO <- function(df, split_ratio = 0.8, B = 100, alpha = 0.05){
  
  LOCO_iteration <- function(D_1, D_2, feature_index, B, alpha){
    library(caret)
    accuracy <- function(model, df) {
      preds <- as.numeric(unlist(predict(model, df)))
      y <- as.numeric(unlist(df["y"]))
      return(c(abs(preds - y)))
    }
    
    sane_model <- glm(y~., data = D_1, family = "binomial")
    
    loco_model <- glm(y~., data = D_1[-feature_index], family = "binomial")
    
    bootstrap_values <- rep(NA,B)
    
    for (b in c(1:B)) {
      boots_idxs <- sample(1:dim(D_2)[1], replace=TRUE)
      D_boots <- D_2[boots_idxs,]
      bootstrap_values[b] <- median(accuracy(loco_model, D_boots) - accuracy(sane_model, D_boots))
    }
    
    bootstrap_values <- sort(bootstrap_values)
    alpha <- alpha / (dim(D_1)[2] - 1)
    
    point_estimate <- mean(bootstrap_values)
    q_1 <- bootstrap_values[ceiling(B*alpha/2)]
    q_2 <- bootstrap_values[ceiling(B*(1-(alpha/2)))]
    
    return(c("feature" = feature_index, "lower" = q_1, "feature_importance" = point_estimate, "upper" = q_2))
    
  }
  split_data <- function(df, split_ratio){
    
    N <- dim(df)[1]
    sample <- sample(1:N, round(split_ratio * dim(df)[1]))
    D_1  <- df[sample,]
    D_2  <- df[-sample,]
    return(list("D_1"=D_1, "D_2"=D_2))
  }
  
  c(D_1,D_2) := split_data(df, split_ratio)
  
  number_of_features <- dim(df)[2]-1
  
  result <- foreach(feature = c(1:number_of_features), .combine = rbind) %dopar% 
    LOCO_iteration(D_1, D_2, feature, B, alpha )
  
  
  return(data.frame(result))
}
