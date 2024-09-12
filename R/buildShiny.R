#' This function builds the minimum files required for Shiny
#'
#' @param counts A count matrix in sparse format: dgCMatrix.
#' @param meta.data Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts. "clusters" must be provided (see celltypeColumn and notes).
#' @param feature.set Set of feature used to calculate dendrogram. Typically highly variable and/or marker genes.
#' @param umap.coords Dimensionality reduction coordiant data.frame with 2 columns. Rownames must be equal to colnames of counts.
#' @param taxonomyDir The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param taxonomyName The file name to assign for the Taxonomy h5ad.
#' @param celltypeColumn Column name corresponding to where the clusters are located (default="cluster")
#' @param cluster_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. cluster_colors can also be provided in the metadata (see notes)
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' @param subsample The number of cells to retain per cluster
#' @param dend Existing dendrogram associated with this taxonomy (e.g., one calculated elsewhere).  If NULL (default) a new dendrogram will be calculated based on the input `feature.set`
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' 
#' Notes: Precomputed clusters must be provided.  In the anndata object these will be stored using the term "cluster".  If celltypeColumn is anything other than cluster, then any existing "cluster" column will be overwritten by celltypeColumn.  Values can be provided without colors and ids (e.g., "cluster") or with them (e.g., "cluster_label" + "cluster_color" + "cluster_id").  In this case cluster_colors is ignored and colors are taken directly from the metadata.  Cluster_id's will be overwritten to match dendrogram order.
#' 
#' @import scrattch.hicat 
#' @import scrattch.io
#' @import feather
#' @import tibble
#' @import dplyr
#' @import Matrix
#' @import pvclust
#' @import anndata
#'
#' @return AIT anndata object in the specified format (only if return.anndata=TRUE)
#'
#' @export
buildShiny = function(AIT.anndata){

  ## Get the data and metadata matrices
  counts    = counts[,kpSub]
  meta.data = meta.data[kpSub,]
  umap.coords = umap.coords[kpSub,]

  ## Create log TPM matrix from counts
  tpm.matrix = scrattch.hicat::logCPM(counts)

  ## Next, generate and output the data matrices. Make sure matrix class is dgCMatrix and not dgTMatrix.  
  sample_id = colnames(tpm.matrix); meta.data$sample_id = sample_id
  gene      = rownames(tpm.matrix)
  cluster   = meta.data$cluster; names(cluster) = rownames(meta.data) ## For get_cl_medians

  ## Counts feather - Not used for shiny, but useful for saving raw data nonetheless
  print("===== Building counts.feather =====")
  counts.tibble = as_tibble(counts)
  counts.tibble = cbind(gene, counts.tibble)
  counts.tibble = as_tibble(counts.tibble)
  write_feather(counts.tibble, file.path(taxonomyDir, "counts.feather"))

  ## Data_t feather
  print("===== Building data_t.feather =====")
  data.tibble = as_tibble(tpm.matrix)
  data.tibble = cbind(gene, data.tibble)
  data.tibble = as_tibble(data.tibble)
  write_feather(data.tibble, file.path(taxonomyDir, "data_t.feather"))
  
  ## Data feather
  print("===== Building data.feather =====")
  norm.data.t = t(as.matrix(tpm.matrix))
  norm.data.t = as_tibble(norm.data.t)
  norm.data.t = cbind(sample_id, norm.data.t)
  norm.data.t = as_tibble(norm.data.t)
  write_feather(norm.data.t, file.path(taxonomyDir, "data.feather"))
  
  ## Output tree
  saveRDS(dend, file.path(taxonomyDir,"dend.RData"))

  ## Output tree order
  outDend = data.frame(cluster = labels(dend), order = 1:length(labels(dend)))
  write.csv(outDend, file.path(taxonomyDir, "ordered_clusters.csv"))
    
  ## Write the desc file.  
  anno_desc = create_desc(meta.data, use_label_columns = TRUE)
  # Subset the desc file to match metadata_names, if provided
  if(!is.null(metadata_names)){
    desc <- anno_desc[match(names(metadata_names), as.character(as.matrix(anno_desc[,1]))),]
    desc[,2] <- as.character(metadata_names[as.character(as.matrix(desc[,1]))])
    desc <- desc[!is.na(desc$base),]  # Remove missing values
    anno_desc <- desc
  }
  write_feather(anno_desc, file.path(taxonomyDir,"desc.feather"))

  ## Minor reformatting of metadata file, then write metadata file
  meta.data$cluster = meta.data$cluster_label; 
  colnames(meta.data)[colnames(meta.data)=="sample_name"] <- "sample_name_old" # Rename "sample_name" to avoid shiny crashing
  if(!is.element("sample_id", colnames(meta.data))){ meta.data$sample_id = meta.data$sample_name_old } ## Sanity check for sample_id
  meta.data$cluster_id <- as.numeric(factor(meta.data$cluster_label,levels=labels(dend))) # Reorder cluster ids to match dendrogram
  write_feather(meta.data, file.path(taxonomyDir,"anno.feather"))

  ## Write the UMAP coordinates.  
  print("===== Building umap/tsne feathers (precalculated) =====")
  tsne      = data.frame(sample_id = rownames(umap.coords),
                         all_x = umap.coords[,1],
                         all_y = umap.coords[,2])
  tsne      = tsne[match(meta.data$sample_id, tsne$sample_id),]
  tsne_desc = data.frame(base = "all",
                         name = "All Cells UMAP")

  write_feather(tsne, file.path(taxonomyDir,"tsne.feather"))
  write_feather(tsne_desc, file.path(taxonomyDir,"tsne_desc.feather"))

  ##
  print("===== Building count, median, sum feathers =====")
  all_clusters = unique(meta.data$cluster_label) ## cluster_id
  count_gt0 = matrix(0, ncol = length(all_clusters), nrow = nrow(data.tibble))
  count_gt1 = sums = medianmat = count_gt0

  ## Compute the number of genes with greater than 0,1 counts, gene sums and medians
  for (i in 1:length(all_clusters)) {
    cluster = all_clusters[i]
    cluster_samples = which(meta.data$cluster_label == cluster)
    cluster_data    = tpm.matrix[,cluster_samples,drop=F]
    cluster_counts  = counts[,colnames(cluster_data),drop=F]
    count_gt0[, i]  = Matrix::rowSums(cluster_counts > 0)
    count_gt1[, i]  = Matrix::rowSums(cluster_counts > 1)
    sums[, i]       = Matrix::rowSums(cluster_counts)
    medianmat[, i]  = apply(cluster_data, 1, median)
  }

  ##
  colnames(count_gt0) = colnames(count_gt1) = colnames(sums) = colnames(medianmat) = all_clusters ## allClust
  count_gt0 = cbind(gene = gene, as.data.frame(count_gt0))
  count_gt1 = cbind(gene = gene, as.data.frame(count_gt1))
  sums = cbind(gene = gene, as.data.frame(sums))
  medianmat = cbind(gene = gene, as.data.frame(medianmat))

  ##
  count_n = meta.data %>% 
              arrange(cluster_id) %>% 
              group_by(cluster_id) %>% 
              summarise(n_cells = n())

  ##
  write_feather(count_gt0, file.path(taxonomyDir,"count_gt0.feather"))
  write_feather(count_gt1, file.path(taxonomyDir,"count_gt1.feather"))
  write_feather(count_n,   file.path(taxonomyDir,"count_n.feather"))
  write_feather(medianmat, file.path(taxonomyDir,"medians.feather"))
  write_feather(sums,      file.path(taxonomyDir,"sums.feather"))
}
