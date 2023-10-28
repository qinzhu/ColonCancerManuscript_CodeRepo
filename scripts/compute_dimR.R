

#' @export
filter_cds <- function(cds, min_detect=1, min_numc_expressed = 10, min_disp_ratio=1){
    numc_expressed <- Matrix::rowSums(exprs(cds) > min_detect)
    g_filter <- numc_expressed >= min_numc_expressed
    sum(g_filter)
    cds_oidx <- cds[g_filter,]
    #cds_oidx@normalized_data_projection <- matrix() # bug in monocle3 cellset dataset class definition
    cds_oidx <- estimateDispersions(cds_oidx)
    disp_subset = dispersionTable(cds_oidx)
    disp_list = subset(disp_subset,(dispersion_empirical/dispersion_fit)>min_disp_ratio)
    disp_genes = disp_list$gene_id
    message(paste0("Final VEG number:", length(disp_genes)))
    cds_oidx <- setOrderingFilter(cds_oidx, disp_genes)
    return(cds_oidx)
}


#' @export
compute_pca_cds <- function(cds, num_dim =100, scvis=NULL, use_order_gene = T, residualModelFormulaStr = NULL, return_type=c("irlba","scvis","proj"), use_raw = F) {
    if(use_raw) {
        FM <- log2(exprs(cds) + 1)
    } else {
        FM <- normalize_expr_data2(cds, "log", 1, use_order_gene = use_order_gene)
    }

    xm <- Matrix::rowMeans(FM)
    xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
    FM <- FM[xsd > 0, ]

    if (!is.null(residualModelFormulaStr)) {
        X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), data = pData(cds), drop.unused.levels = TRUE)
        fit <- limma::lmFit(FM, X.model_mat)
        beta <- fit$coefficients[, -1, drop = FALSE]
        beta[is.na(beta)] <- 0
        FM <- as.matrix(FM) - beta %*% t(as.matrix(X.model_mat[, -1]))
    }
    #print(class(FM))
    irlba_res <- prcomp_irlba2(t(as.matrix(FM)), n = min(num_dim, min(dim(FM)) - 1), center = TRUE, scale. = TRUE) # Change to recent sparse version if necessary
    irlba_sdev <- irlba_res$sdev
    names(irlba_sdev) <- paste0("PC",1:length(irlba_sdev))
    pca_proj <- as.data.frame(irlba_res$x)
    rownames(pca_proj) <- colnames(cds)

    if(return_type == "scvis") {
        scvis@pca <- pca_proj
        return(scvis)
    } else if(return_type == "irlba"){
        return(irlba_res)
    }else {
        return(pca_proj)
    }
}

#' @export
prcomp_irlba2 <- function (x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE,
                           ...)
{
    a <- names(as.list(match.call()))
    ans <- list(scale = scale.)
    if ("tol" %in% a)
        warning("The `tol` truncation argument from `prcomp` is not supported by\n`prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to\ncontrol that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
    if (!is.matrix(x))
        x <- as.matrix(x)
    args <- list(A = x, nv = n)
    if (is.logical(center)) {
        if (center)
            args$center <- colMeans(x)
    }
    else args$center <- center
    if (is.logical(scale.)) {
        if (is.numeric(args$center)) {
            f <- function(i) sqrt(sum((x[, i] - args$center[i])^2)/(nrow(x) -
                                                                        1L))
            scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE)
            ans$scale_val <- scale.
            if (ans$scale)
                ans$totalvar <- ncol(x)
            else ans$totalvar <- sum(scale.^2)
        }
        else {
            if (ans$scale) {
                scale. <- apply(x, 2L, function(v) sqrt(sum(v^2)/max(1,
                                                                     length(v) - 1L)))
                ans$scale_val <- scale.
                f <- function(i) sqrt(sum((x[, i]/scale.[i])^2)/(nrow(x) -
                                                                     1L))
                ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi,
                                           USE.NAMES = FALSE)^2)
            }
            else {
                f <- function(i) sum(x[, i]^2)/(nrow(x) - 1L)
                ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi,
                                           USE.NAMES = FALSE))
            }
        }
        if (ans$scale)
            args$scale <- scale.
    }
    else {
        args$scale <- scale.
        f <- function(i) sqrt(sum((x[, i]/scale.[i])^2)/(nrow(x) -
                                                             1L))
        ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE))
    }
    if (!missing(...))
        args <- c(args, list(...))
    s <- do.call(irlba, args = args)
    ans$sdev <- s$d/sqrt(max(1, nrow(x) - 1))
    ans$rotation <- s$v
    colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)),
                                    sep = "")
    ans$center <- args$center
    if (retx) {
        ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
        colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)),
                                 sep = "")
    }
    class(ans) <- c("irlba_prcomp", "prcomp")
    ans
}


#' @export
compute_umap_pca <- function(pca_proj, use_dim = 50, 
                             n_component=2, 
                             metric = "cosine",
                             min_dist = 0.1,
                             n_neighbors = 15L,
                             fast_sgd = FALSE,
                             nn_method = "annoy", 
                             cores=1,
                             verbose=T, ...) {
    umap_proj <- uwot::umap(as.matrix(pca_proj[, 1:use_dim]),
                            n_components = n_component,
                            metric = metric,
                            min_dist = min_dist,
                            n_neighbors = n_neighbors,
                            fast_sgd = fast_sgd,
                            n_threads=cores,
                            verbose=verbose,
                            nn_method = nn_method,
                            ...)
    colnames(umap_proj) <- paste0("UMAP_", 1:n_component)
    rownames(umap_proj) <- rownames(pca_proj)
    return(umap_proj)
}


#' @export
compute_tsne_pca <- function(pca_proj, use_dim = 30, n_component = 2, perplexity = 30, ...) {
  tsne_res <- Rtsne::Rtsne(as.matrix(pca_proj[, 1:use_dim]), 
                           dims = n_component, pca = F, perplexity = perplexity, 
                           ...)
  tsne_proj <- as.data.frame(tsne_res$Y[, 1:n_component])
  colnames(tsne_proj) <- paste0("TSNE_", 1:n_component)
  rownames(tsne_proj) <- rownames(pca_proj)
  return(tsne_proj)
}



louvain_clustering <- function(data, k = 20) {
    neighborMatrix <- RANN::nn2(data, data, k + 1, searchtype = "standard")[[1]][,-1]
    links <- monocle:::jaccard_coeff(neighborMatrix, F)
    links <- links[links[,1]>0, ]
    relations <- as.data.frame(links)
    colnames(relations)<- c("from","to","weight")
    g <- igraph::graph.data.frame(relations, directed=FALSE)
    Q <- igraph::cluster_louvain(g)
    return(factor(igraph::membership(Q)))
}


#' @export
monocle_UMAP <- function (X, python_home = system("which python", intern = TRUE),
          log = TRUE, n_neighbors = 15L, n_component = 2L, metric = "correlation",
          n_epochs = NULL, negative_sample_rate = 5L, learning_rate = 1,
          init = "spectral", min_dist = 0.1, spread = 1, set_op_mix_ratio = 1,
          local_connectivity = 1L, repulsion_strength = 1, a = NULL,
          b = NULL, random_state = 0L, metric_kwds = reticulate::dict(),
          angular_rp_forest = FALSE, target_n_neighbors = -1L, target_metric = "categorical",
          target_metric_kwds = reticulate::dict(), target_weight = 0.5,
          transform_seed = 42L, verbose = FALSE, return_all = FALSE)
{
    reticulate::use_python(python_home)
    tryCatch({
        reticulate::import("umap")
    }, warning = function(w) {
    }, error = function(e) {
        print(e)
        stop("please pass the python home directory where umap is installed with python_home argument!")
    }, finally = {
    })
    reticulate::source_python(paste(system.file(package = "monocle"),
                                    "umap.py", sep = "/"))
    if (length(grep("Matrix", class(X))) == 0) {
        X <- as(as.matrix(X), "TsparseMatrix")
    }
    else {
        X <- as(X, "TsparseMatrix")
    }
    i <- as.integer(X@i)
    j <- as.integer(X@j)
    if (log) {
        val <- log(X@x + 1)
    }
    else {
        val <- X@x
    }
    dim <- as.integer(X@Dim)
    if (is.null(n_epochs) == F) {
        n_epochs <- as.integer(n_epochs)
    }
    if (is.null(a) == F) {
        a <- as.numeric(a)
    }
    if (is.null(b) == F) {
        n_epochs <- as.numeric(b)
    }
    if (is.list(metric_kwds) == F) {
        metric_kwds <- reticulate::dict()
    }
    else {
        metric_kwds <- reticulate::dict(metric_kwds)
    }
    if (is.list(target_metric_kwds) == F) {
        target_metric_kwds <- reticulate::dict()
    }
    else {
        target_metric_kwds <- reticulate::dict(target_metric_kwds)
    }
    umap_res <- umap(i, j, val, dim, as.integer(n_neighbors),
                     as.integer(n_component), as.character(metric), n_epochs,
                     as.integer(negative_sample_rate), as.numeric(learning_rate),
                     as.character(init), as.numeric(min_dist), as.numeric(spread),
                     as.numeric(set_op_mix_ratio), as.integer(local_connectivity),
                     as.numeric(repulsion_strength), a, b, as.integer(random_state),
                     metric_kwds, as.logical(angular_rp_forest), as.integer(target_n_neighbors),
                     as.character(target_metric), target_metric_kwds, as.numeric(target_weight),
                     as.integer(transform_seed), as.logical(verbose))
    if (return_all) {
        return(umap_res)
    }
    else {
        umap_res$embedding_
    }
}

#' @export
normalize_expr_data2 <-function (cds, norm_method = c("log", "vstExprs", "none"), pseudo_expr = 1, use_order_gene = TRUE, relative_expr = TRUE)
{
    FM <- exprs(cds)
    if(!use_order_gene) {
        fData(cds)$use_for_ordering <- NULL
    }
    if (is.null(fData(cds)$use_for_ordering) == FALSE && nrow(subset(fData(cds),
                                                                     use_for_ordering == TRUE)) > 0) {
        FM <- FM[fData(cds)$use_for_ordering, ]
    }
    norm_method <- match.arg(norm_method)
    if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        if (is.null(pseudo_expr)) {
            if (norm_method == "log")
                pseudo_expr = 1
            else pseudo_expr = 0
        }
        checkSizeFactors(cds)
        if (norm_method == "vstExprs") {
            if (relative_expr == FALSE)
                message("Warning: relative_expr is ignored when using norm_method == 'vstExprs'")
            if (is.null(fData(cds)$use_for_ordering) == FALSE &&
                nrow(subset(fData(cds), use_for_ordering == TRUE)) >
                0) {
                VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering,
                                       ], round_vals = FALSE)
            }
            else {
                VST_FM <- vstExprs(cds, round_vals = FALSE)
            }
            if (is.null(VST_FM) == FALSE) {
                FM <- VST_FM
            }
            else {
                stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
            }
        }
        else if (norm_method == "log") {
            if (relative_expr)
                FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
            if (is.null(pseudo_expr))
                pseudo_expr <- 1
            FM <- FM + pseudo_expr
            FM <- log2(FM)
        }
        else if (norm_method == "none") {
            FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
            FM <- FM + pseudo_expr
        }
    }
    else if (cds@expressionFamily@vfamily == "binomialff") {
        if (norm_method == "none") {
            ncounts <- FM > 0
            ncounts[ncounts != 0] <- 1
            FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
        }
        else {
            stop("Error: the only normalization method supported with binomial data is 'none'")
        }
    }
    else if (cds@expressionFamily@vfamily == "Tobit") {
        FM <- FM + pseudo_expr
        if (norm_method == "none") {
        }
        else if (norm_method == "log") {
            FM <- log2(FM)
        }
        else {
            stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
        }
    }
    else if (cds@expressionFamily@vfamily == "gaussianff") {
        if (norm_method == "none") {
            FM <- FM + pseudo_expr
        }
        else {
            stop("Error: the only normalization method supported with gaussian data is 'none'")
        }
    }
    return(FM)
}

#' @export
checkSizeFactors <- function (cds)
{
    if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        if (is.null(sizeFactors(cds))) {
            stop("Error: you must call estimateSizeFactors() before calling this function.")
        }
        if (sum(is.na(sizeFactors(cds))) > 0) {
            stop("Error: one or more cells has a size factor of NA.")
        }
    }
}


louvain_R <- function(X, python_home = system('which python', intern = TRUE), 
                      partition_method = 'CPMVertexPartition', 
                      initial_membership = NULL, 
                      weights = NULL, 
                      res = 0.6, 
                      node_sizes = NULL, 
                      random_seed = 0L, 
                      verbose = FALSE,
                      return_all = FALSE,
                      louvain_path = NULL) {
    
    reticulate::use_python(python_home)
    
    tryCatch({
        reticulate::import("louvain")
    }, warning = function(w) {
    }, error = function(e) {
        print (e)
        stop('please pass the python home directory where louvain is installed with python_home argument!')
    }, finally = {
    })
    
    reticulate::source_python(louvain_path)
    # X <- Matrix::t(X)
    if(length(grep('Matrix', class(X))) == 0){
        X <- as(as.matrix(X), 'TsparseMatrix')
    } else {
        X <- as(X, 'TsparseMatrix')
    }
    
    i <- as.integer(X@i)
    j <- as.integer(X@j)
    val <- X@x
    
    dim <- as.integer(X@Dim)
    
    if(is.null(partition_method) == F) {
        partition_method <- as.character(partition_method)
    }
    if(!is.null(random_seed)) {
        random_seed <- as.integer(random_seed)
    }
    # if(is.null(initial_membership) == F) { #initial_membership (list of int) 
    #   a <- as.numeric(a)
    # }
    # if(is.null(weights) == F) { # weights (list of double, or edge attribute) 
    #   n_epochs <- as.numeric(b)
    # }
    # if(is.null(res) == F) { # node_sizes (list of int, or vertex attribute)
    #   metric_kwds <- reticulate::dict()
    # } 
    # if(is.null(node_sizes) == F) { # resolution_parameter (double)
    #   metric_kwds <- reticulate::dict(metric_kwds)
    # }
    
    louvain_res <- louvain(i, j, val, dim, 
                           as.character(partition_method), 
                           initial_membership, 
                           weights, 
                           as.numeric(res),
                           node_sizes,
                           random_seed,
                           as.logical(verbose))
    
    if(return_all) {
        return(louvain_res)
    } else {
        list(membership = louvain_res$membership + 1, modularity = louvain_res$modularity)
    }
}


#' @export
louvain_clus <- function (data, k = 20, weight = F, louvain_iter = 1, resolution = NULL, 
                          random_seed = 0L, verbose = F, ...) 
{
    extra_arguments <- list(...)
    cell_names <- row.names(data)
    if (is.data.frame(data)) 
        data <- as.matrix(data)
    if (!is.matrix(data)) 
        stop("Wrong input data, should be a data frame of matrix!")
    if (k < 1) {
        stop("k must be a positive integer!")
    }
    else if (k > nrow(data) - 2) {
        stop("RANN counts the point itself, k must be smaller than\nthe total number of points - 1 (all other points) - 1 (itself)!")
    }
    if (verbose) {
        message("Run kNN based graph clustering starts:", "\n", 
                "  -Input data of ", nrow(data), " rows and ", ncol(data), 
                " columns", "\n", "  -k is set to ", k)
    }
    if (verbose) {
        cat("  Finding nearest neighbors...")
    }
    t1 <- system.time(tmp <- RANN::nn2(data, data, k + 1, searchtype = "standard"))
    neighborMatrix <- tmp[[1]][, -1]
    distMatrix <- tmp[[2]][, -1]
    if (verbose) {
        cat("DONE ~", t1[3], "s\n", " Compute jaccard coefficient between nearest-neighbor sets ...")
    }
    t2 <- system.time(links <- monocle:::jaccard_coeff(neighborMatrix, 
                                                       weight))
    if (verbose) {
        cat("DONE ~", t2[3], "s\n", " Build undirected graph from the weighted links ...")
    }
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    relations$from <- cell_names[relations$from]
    relations$to <- cell_names[relations$to]
    t3 <- system.time(g <- igraph::graph.data.frame(relations, 
                                                    directed = FALSE))
    if (verbose) {
        cat("DONE ~", t3[3], "s\n", " Run louvain clustering on the graph ...\n")
    }
    t_start <- Sys.time()
    Qp <- -1
    optim_res <- NULL
    best_max_resolution <- "No resolution"
    if (louvain_iter >= 2) {
        random_seed <- NULL
    }
    for (iter in 1:louvain_iter) {
        if (verbose) {
            cat("Running louvain iteration ", iter, "...\n")
        }
        if (!is.null(resolution)) {
            for (i in 1:length(resolution)) {
                cur_resolution <- resolution[i]
                louvain_args <- c(list(X = igraph::get.adjacency(g), 
                                       res = as.numeric(cur_resolution), random_seed = random_seed, 
                                       verbose = verbose), extra_arguments[names(extra_arguments) %in% 
                                                                               c("python_home", "partition_method", "initial_membership", 
                                                                                 "weights", "node_sizes", "return_all", "louvain_path")])
                Q <- do.call(louvain_R, louvain_args)
                Qt <- max(Q$modularity)
                if (verbose) {
                    message("Current iteration is ", iter, "; current resolution is ", 
                            cur_resolution, "; Modularity is ", Qt, "; Number of clusters are ", 
                            max(Q$membership))
                }
                if (Qt > Qp) {
                    optim_res <- Q
                    Qp <- Qt
                    best_max_resolution <- cur_resolution
                }
            }
        }
        else {
            Q <- igraph::cluster_louvain(g)
        }
        if (is.null(optim_res)) {
            Qp <- max(Q$modularity)
            optim_res <- Q
        }
        else {
            Qt <- max(Q$modularity)
            if (Qt > Qp) {
                optim_res <- Q
                Qp <- Qt
            }
        }
    }
    if (verbose) 
        message("Maximal modularity is ", Qp, "; corresponding resolution is ", 
                best_max_resolution)
    t_end <- Sys.time()
    if (verbose) {
        message("\nRun kNN based graph clustering DONE, totally takes ", 
                t_end - t_start, " s.")
        cat("  -Number of clusters:", length(unique(igraph::membership(optim_res))), 
            "\n")
    }
    if (igraph::vcount(g) < 3000) {
        coord <- NULL
        edge_links <- NULL
    }
    else {
        coord <- NULL
        edge_links <- NULL
    }
    V(g)$names <- as.character(V(g))
    #return(list(g = g, relations = relations, distMatrix = distMatrix, 
    #            coord = coord, edge_links = edge_links, optim_res = optim_res))
    return(factor(igraph::membership(optim_res)))
}


dens_clus<-function(data, rho_threshold=5, delta_threshold=10, peaks=NULL){
    dataDist <- dist(data)
    dataClust <- densityClust::densityClust(dataDist, gaussian = T)
    dataClust <- densityClust::findClusters(dataClust, 
                                            rho = rho_threshold, delta = delta_threshold, 
                                            peaks = peaks)
    return(factor(dataClust$clusters))
}



#' @export
leiden_clus <- function (embedding, k = 30) {
  require(RANN)
  require(leiden)
  snn <- RANN::nn2(embedding, k=k)$nn.idx
  adjacency_matrix <- matrix(0L, nrow(embedding), nrow(embedding))
  rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- rownames(embedding)
  for(ii in 1:nrow(embedding)) {
    adjacency_matrix[ii,rownames(embedding)[snn[ii,]]] <- 1L
  }
  partition <- leiden(adjacency_matrix)
  return(partition)
}



