library(Seurat)
library(dplyr)
library(SeuratData)
library(Matrix)
library(CelliD)
library(dnet)
library(reshape2)
library(infotheo)
library(igraph)
library(anticlust)
# library(singleCellHaystack)
library(multimode)
library(RANN)
library(future)
library(future.apply)
library(philentropy)
library(MASS)
#####################################
# 1. compute MCA embeddings
#
compute.mca <- function(object, dims.use = 1:10, genes.use = rownames(object)){

        genes.use.1 <- intersect(genes.use, rownames(object))
        object <- object %>%
                NormalizeData() %>%
                # ScaleData(features = rownames(object)) %>%
                RunMCA(features = genes.use.1)
        genes.use.update <- intersect(genes.use.1, rownames(object@reductions$mca@feature.loadings))
        coembed <- rbind(object@reductions$mca@feature.loadings[genes.use.update, dims.use],
                 object@reductions$mca@cell.embeddings[, dims.use])
        return(coembed)
}

# 2. compute density of gene sets of interest
# 2.1 compute grid point coordinates
compute.grid.coords <- function(coembed, genes.use, n.grids = 100){
        coembed <- scale(coembed[genes.use, ])
        cl <- balanced_clustering(coembed, K = n.grids)
        centroid.coords <- aggregate(coembed, list(cl), mean)[, -1]
        return(centroid.coords)
}

# 2.2 compute KL-divergence
# some are adapted from https://github.com/alexisvdb/singleCellHaystack/
compute.kld <- function(coembed, genes.use, n.grids = 100, gene.set.list, gene.set.cutoff = 3,
                        n.times = 100){
        coembed.scale <- scale(coembed[genes.use, ])
        grid.co <- compute.grid.coords(coembed = coembed,
                                       genes.use = genes.use,
                                       n.grids = n.grids)
        # for background
        dist.to.grid <- vectorized_pdist(A = as.matrix(coembed.scale), B = as.matrix(grid.co))
        bandwidth <- median(apply(dist.to.grid,1,min))
        dist.to.grid.norm <- dist.to.grid / bandwidth
        density.contributions <-
                exp(-dist.to.grid.norm * dist.to.grid.norm / 2)
        Q <- compute.db(density.df = density.contributions)

        # compute kld for each gene set and randomly sampled gene sets of the same size
        gene.set.names <- rep(NA, length(gene.set.list))
        klds <- rep(0, length(gene.set.list))
        len.gl <- rep(0, length(gene.set.list))
        rklds.avg <- rep(0, length(gene.set.list))
        rklds.sd <- rep(0, length(gene.set.list))
        pvalues <- rep(0, length(gene.set.list))
        for (i in 1:length(gene.set.list)){
                gene.set.name <- names(gene.set.list)[i]
                gene.set <- intersect(genes.use, gene.set.list[[i]])
#                print(gene.set)
                len.gene.set <- length(gene.set)
                if (len.gene.set < gene.set.cutoff){
                        next
                }
                gene.set.names[i] <- gene.set.name
                len.gl[i] <- len.gene.set
                P <- compute.db(density.df = density.contributions[gene.set, ])
                kld <- sum(P * log(P/Q))
                klds[i] <- log(kld)
                # sample
                rkld <- sapply(1:n.times, function(x)sample.kld(density.df = density.contributions,
                                                    ref = Q,
                                                    len.gene.set = len.gene.set))
                rkld.avg <- mean(log(rkld))
                rkld.sd <- sd(log(rkld))
                rklds.avg[i] <- rkld.avg
                rklds.sd[i] <- rkld.sd
                pvalue <- pnorm(log(kld), rkld.avg, rkld.sd, lower.tail = FALSE)
                pvalues[i] <- pvalue
        }
        out <- as.data.frame(cbind(gene.set.names,
                                   klds,
                                   len.gl,
                                   rklds.avg,
                                   rklds.sd,
                                   pvalues))
        out <- out[complete.cases(out), ]
        colnames(out) <- c("gene.set",
                           "kld",
                           "gene.set.length",
                           "rkld.mean",
                           "rkld.sd",
                           "p.value")
        out$p.adj <- p.adjust(p = out$p.value, method = "fdr")

        return(out)
}
# from an excellent post: https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
# enhanced the speed
vectorized_pdist <- function(A,B){
    an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
    bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))

    m = nrow(A)
    n = nrow(B)

    tmp = matrix(rep(an, n), nrow=m)
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    sqrt( tmp - 2 * tcrossprod(A,B) )
}

compute.db <- function(density.df){
        Q <- apply(density.df, 2, sum)
        pseudo <- 1e-300
        Q <- Q + pseudo
        Q <- Q / sum(Q)
        return(Q)
}


sample.kld <- function(density.df, ref, len.gene.set){
        idx <- sample(nrow(density.df), len.gene.set, replace = F)
        rP <- compute.db(density.df = density.df[idx, ])
        rkld.vec <- sum(rP * log(rP/ref))
        return(rkld.vec)
}

# 3. compute nearest neighbor graph for genes and cells
el_nn_search <- function(nn2_out){
  n.df <- nn2_out$nn.idx
  n.df <- cbind(1:nrow(n.df), n.df)
  el <- cbind(n.df[, 1], c(n.df[, -1]))
  return(el)
}

compute.nn.edges <- function(coembed, nn.use = 300){
        # distmat <- dist(coembed)

        # fnn <- FindNeighbors(
        #                     object = distmat,
        #                     k.param = nn.use,
        #                     return.neighbor = FALSE,
        # )
        nbrs <- nn2(coembed, k = nn.use)
        el_idx <- el_nn_search(nn2_out = nbrs)
        el_nn <- cbind(rownames(coembed)[el_idx[, 1]], rownames(coembed)[el_idx[, 2]])
        # el_nn <- get.edgelist(graph.adjacency(fnn$nn))
        return(el_nn)
}


# 4. compute label propagation from gene set to cells
seed.mat <- function(gene_set, graph.use){
        gs <- intersect(gene_set, names(V(graph.use)))
        ss <- data.frame(n = names(V(graph.use)))
        rownames(ss) <- ss$n
        ss$weight <- ifelse(rownames(ss) %in% gs, 1, 0)
        ss$n <- NULL
        # ss <- as.matrix(ss)
        return(ss)
}

seed.mat.list <- function(gene_set_list, graph.use){
        sml <- future.apply::future_lapply(names(gene_set_list), 
                                           function(x) seed.mat(gene_set = gene_set_list[[x]], 
                                                                graph.use = graph.use))
        sm <- do.call(cbind, sml)
        colnames(sm) <- names(gene_set_list)
        sm <- as.matrix(sm)
        return(sm)
}

run.rwr <- function(el, gene_set, cells){
        g <- graph_from_edgelist(as.matrix(el), directed = F)
        g <- simplify(g)
        ss <- seed.mat(gene_set = gene_set, graph.use = g)
        rwr <- dRWR(g, setSeeds = ss, parallel = F, verbose = F)
        rownames(rwr) <- rownames(ss)
        cell_vec <- rwr[cells, 1]
        cell_vec <- cell_vec / sum(cell_vec)
        return(cell_vec)
}

run.rwr.list <- function(el, gene_set_list, cells){
        g <- graph_from_edgelist(as.matrix(el), directed = F)
        g <- simplify(g)
        ss <- seed.mat.list(gene_set_list = gene_set_list, graph.use = g)
        rwrl <- dRWR(g, setSeeds = ss, parallel = T, verbose = F)
        rownames(rwrl) <- rownames(ss)
        cell_rwr <- rwrl[cells, ]
        colnames(cell_rwr) <- colnames(ss)
        cell_rwr <- cell_rwr[, colSums(is.na(cell_rwr)) == 0]
        cell_rwr_norm <- cell_rwr %*% diag(1/colSums(cell_rwr))

        rownames(cell_rwr_norm) <- rownames(cell_rwr)
        colnames(cell_rwr_norm) <- colnames(cell_rwr)
        cell_rwr_norm <- as.data.frame(cell_rwr_norm)
        return(cell_rwr_norm)
}


compute.cell.label <- function(cell_vec){
        nm <- names(cell_vec)
        m <- locmodes(cell_vec, mod0 = 2)
        split.value <- m$locations[2]
        cell.label <- ifelse(cell_vec < split.value, "negative", "positive")
        return(cell.label)
#        mt <- modetest

}

compute.cell.label.df <- function(cell_df){
        cell.labels <- future_apply(cell_df, 
                                    MARGIN = 2, 
                                    function(x) {compute.cell.label(x)
        })
        return(cell.labels)

}



compute.multimode.p <- function(cell_vec){
        nm <- names(cell_vec)
        m <- modetest(cell_vec)
        p.value <- m$p.value
        return(p.value)

}

compute.multimode.p.df <- function(cell_df){
        p.values <- future_apply(cell_df, 
                                    MARGIN = 2, 
                                    function(x) {compute.multimode.p(x)
        })
        names(p.values) <- colnames(cell_df)
        return(p.values)
}

# 5. compute the specificity of gene set when cluster information is available

# inspired by https://github.com/FloWuenne/scFunctions/blob/0d9ea609fa72210a151f7270e61bdee008e8fc88/R/calculate_rrs.R
compute.jsd <- function(x, y){
        input_df <- rbind(x, y)
        jsd_divergence <- suppressMessages(philentropy::JSD(input_df))
        jsd_distance <- 1-sqrt(jsd_divergence)
        return(jsd_distance)
}

compute.spec.single <- function(vec, positive, cell_df){
        num <- ifelse(vec == positive, 1, 0)
        num <- num / sum(num)
        spec.single <- future_apply(cell_df, 
                                    MARGIN = 2, 
                                    function(x) {compute.jsd(x = x, y = num)})
        names(spec.single) <- colnames(cell_df)
        return(spec.single)
}


compute.spec <- function(cell_df, metadata, cell_group){
        gene.sets <- colnames(cell_df)
        cells <- rownames(cell_df)
        metadata <- metadata[cells, ]
        cell_groups <- unique(as.character(metadata[[cell_group]]))

        jsd <- future.apply::future_lapply(cell_groups, 
                                           function(x) 
                                                   {compute.spec.single(vec = metadata[[cell_group]], 
                                                                        positive = x, 
                                                                        cell_df = cell_df)})
        jsd.df <- do.call(cbind, jsd)
        colnames(jsd.df) <- cell_groups
        return(jsd.df)
}

# find spatially 

compute.spatial.kld <- function(spatial.coords, weight_vec, n = 10){
        bg.weight <- rep(1/nrow(spatial.coords), nrow(spatial.coords))
        bg.dens <- kde2d.weighted(x = spatial.coords[, 1], 
                                  y = spatial.coords[, 2], 
                                  w = bg.weight, 
                                  n = n)
        Q <- c(bg.dens$z) + 1e-300
        dens <- kde2d.weighted(x = spatial.coords[, 1], 
                               y = spatial.coords[, 2], 
                               w = weight_vec, 
                               n = n)
        P <- c(dens$z)
        spatial.kld <- sum(P * log(P/Q))
        return(spatial.kld)
}

compute.spatial.kld.df <- function(spatial.coords, weight_df, n = 10){
        weight_df <- weight_df[rownames(spatial.coords), ] # make sure the order is the same
        klds <- future_apply(weight_df, 
                             MARGIN = 2, 
                             function(x) {compute.spatial.kld(spatial.coords = spatial.coords, 
                                                              weight_vec = x, 
                                                              n = n)})
        names(klds) <- colnames(weight_df)
        return(klds)
}

# https://stat.ethz.ch/pipermail/r-help/2006-June/107405.html
kde2d.weighted <- function (x, y, w, h, n, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
      stop("data vectors must be the same length")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}

