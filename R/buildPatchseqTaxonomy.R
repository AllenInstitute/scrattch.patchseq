#' Save marker genes for patchSeqQC
#'
#' This function saves all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm) to AIT.anndata$uns. This is only used for patch-seq analysis.  Requirements for input include:
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param mode.name A name to identify the taxonomy mode to add QC metrics for (default is the current mode). If the mode does not yet exist, buildPatchseqTaxonomy will call buildTaxonomyMode using the sampled data from the `class.column`, `subclass.column`, and `off.target.types`.
#' @param subclass.column Column name corresponding to the moderate-resolution cell types used for the cell types of interest (default = "subclass_label").
#' @param class.column Column name corresponding to the low-resolution cell types used for the off-target cell types (default = "class_label").
#' @param off.target.types A character vector of off-target (also known as 'contamination') cell types.  This must include at least one of the cell types found in "class.column" and/or "subclass.column" (both columns are checked)

#' @param subclass.subsample The number of cells to retain for PatchseqQC contamination calculation (default = 100, probably no need to change).
#' @param num.markers The maximum number of markers to calculate per node per direction (default = 50)
#' @param taxonomyDir The location to save shiny output (default = current working directory).
#' @param add.dendrogram.markers If TRUE (default=TRUE), will also add dendrogram markers to prep the taxonomy for tree mapping. Default is TRUE because the membership values calculated here as well as the associated tree mapping capabilities is required for a subset of QL metrics, including KL divergence calculations.
#' @param ... Additional variables to be passed to `addDendrogramMarkers` or `buildTaxonomyMode`
#' 
#' The following variables are added to AIT.anndata$uns:  
#' $dend[[mode.name]]  
#' $filter[[mode.name]]  
#' $QC_markers[[mode.name]]  
#' ...$markers,  
#' ...$countsQC,  
#' ...$cpmQC,  
#' ...$classBr,  
#' ...$subclassF,  
#' ...$allMarkers  
#' $memb[[mode.name]]  # Only if add.dendrogram.markers=TRUE
#' ...$memb.ref,  
#' ...$map.df.ref  
#' 
#' @import patchseqtools
#' @import scrattch.hicat
#' @import reticulate
#'
#' @return AIT.anndata An updated AIT.anndata variable with the above content added to AIT.anndata$uns for the relevant mode.name.
#'
#' @export
addPatchseqQCMetrics = function(AIT.anndata,
                                mode.name = AIT.anndata$uns$mode, 
                                subclass.column = "subclass_label",
                                class.column = "class_label",
                                off.target.types = NULL,
                                subclass.subsample = 100,
                                num.markers = 50,
                                taxonomyDir = file.path(AIT.anndata$uns$taxonomyDir),
                                add.dendrogram.markers = TRUE,
                                ...
){
  
  ## FIRST DO SOME CHECKS OF THE INPUT VARIABLES

  ## Check mode is legal
  if(!is.character(mode.name)) stop("mode.name must be a character string.")
  mode.name <- mode.name[1]

  ## Check that required fields for patchSeqQC exist 
  if(!is.element(subclass.column, colnames(AIT.anndata$obs))){stop(paste(subclass.column,"is not a column in the metadata data frame."))}
  if(!is.element(class.column, colnames(AIT.anndata$obs))){stop(paste(class.column,"is not a column in the metadata data frame."))}
  if(!dir.exists(file.path(taxonomyDir))){stop("Specified taxonomy folder does not exist.")}
  
  ## NEXT SUBSAMPLE THE DATA FOR QC AND DEFINE OFF-TARGET TYPES

  ## Subsample and filter metadata and data
  kpSamp2  = subsampleCells(AIT.anndata$obs[,subclass.column], subclass.subsample)
  goodSamp = !is.na(AIT.anndata$obs$class_label)    # For back-compatibility; usually not used
  kpSamp2  = kpSamp2 & goodSamp                     # For back-compatibility; usually not used
  annoQC   = AIT.anndata$obs[kpSamp2,]
  annoQC$subclass_label = annoQC[,subclass.column]  # For compatibility with existing code.
  annoQC$class_label = annoQC[,class.column]        # For compatibility with existing code.
  annoQC$subclass_label = make.names(annoQC$subclass_label)
  datQC    = as.matrix(Matrix::t(AIT.anndata$raw$X))[,kpSamp2]
  colnames(datQC) = AIT.anndata$obs_names[kpSamp2]; rownames(datQC) = AIT.anndata$var_names

  ## Determine off.target.types if not provided
  if((mode.name %in% names(AIT.anndata$uns$filter))&(is.null(off.target.types))){
    omit <- AIT.anndata$uns$filter[kpSamp2]
    off.target.types <- c(setdiff(annoQC$class_label,annoQC$class_label[omit]),
                          setdiff(annoQC$subclass_label,annoQC$subclass_label[omit]))
  } 
  ## Confirm if off target types are valid and define them
  if(is.null(off.target.types)){
    stop("off.target.types must be provided for new modes and for modes containing all subclasses.")
  }
  offTarget = is.element(annoQC$class_label, off.target.types) | is.element(annoQC$subclass_label, off.target.types)
  if(sum(offTarget)==0){
    stop("No valid off-target classes or subclasses are provided. Please update off.target.types accordingly.")
  }
  
  ## DEFINE QC PARAMETERS AND ASSOCIATED MARKER GENES ANDSAVE TO uns
  
  ## Define class and subclass vectors
  classBr   = annoQC$subclass_label
  classBr[!offTarget] = annoQC$class_label[!offTarget]
  classBr   = factor(classBr)
  subclassF = factor(annoQC$subclass_label)
  
  print("Define and output marker genes for each broad class and off-target subclass.") 
  ## -- These are selected using some reasonable approach that could probably be improved, if needed.    
  markers    = defineClassMarkers(datQC, subclassF, classBr, numMarkers = 50)
  allMarkers = unique(unlist(markers))
  rownames(datQC) = make.names(rownames(datQC))
  countsQC   = datQC[allMarkers,]
  cpmQC      = cpm(datQC)[allMarkers,]  ## Only use of scrattch.hicat in this function

  ## Save patchseqQC information to uns
  AIT.anndata$uns$QC_markers[[mode.name]] = list("allMarkers" = allMarkers,
                                                 "markers"    = markers,
                                                 "countsQC"   = countsQC,
                                                 "cpmQC"      = cpmQC,
                                                 "classBr"    = classBr,
                                                 "subclassF"  = subclassF,
                                                 "qc_samples" = colnames(countsQC),
                                                 "qc_genes"   = rownames(countsQC))
  
  ## CREATE THE TAXONOMY MODE, IF NEEDED
  
  if(!(mode.name %in% names(AIT.anndata$uns$filter))){
    ## Filter out off target cells along with additional cells beyond those subsampled
    ## ---  If the mode.name is new, this filter is used for the "retain" input for buildTaxonomyMode
    filter = is.element(AIT.anndata$obs$class_label, off.target.types) | is.element(AIT.anndata$obs$subclass_label, off.target.types)
    filter = !((!filter) & ((subsampleCells(AIT.anndata$obs$cluster_label,subsample)))) # NEW, for subsampling
    
    ## Define the taxonomy mode here. 
    ## --- There are a LOT of parameters hidden in the '...', but the key piece here is we use the off.target.types for filtering
    AIT.anndata <- buildTaxonomyMode(AIT.anndata = AIT.anndata, mode.name = mode.name, retain.cells = !filter, 
                                     retain.clusters = NULL, add.dendrogram.markers = FALSE, ...)
  }  
  
  ## Add dendrogram markers and membership tables, if requested
  if(add.dendrogram.markers){
    ## Determine the cluster column and warn if not "cluster_id"
    hierarchy = AIT.anndata$uns$hierarchy
    hierarchy = hierarchy[order(unlist(hierarchy))]
    if(is.null(hierarchy)) stop("Hierarchy must be included in the standard AIT mode in proper format to create a mode.  Please run checkTaxonomy().")
    celltypeColumn = names(hierarchy)[length(hierarchy)][[1]]
    if(celltypeColumn!="cluster_id") warning("AIT schema requires clusters to be in 'cluster_id' slot. We recommend calling the finest level of the hierarch as 'cluster_id'.")
    
    print("===== Adding dendrogram markers and membership tables for tree mapping =====")
    tryCatch({
      AIT.anndata = addDendrogramMarkers(AIT.anndata, 
                                         mode=mode.name,  # "standard",     # NOTE: WE MIGHT NEED THIS FOR "standard" TOO, but skip for now. 
                                         celltypeColumn = celltypeColumn,
                                         ...)
    }, error = function(e) {
      print("===== Error adding dendrogram markers. Skipping this step. =====")
      print(e)
    })
  }
  
  
  ## OUTPUT AND RETURN THE RESULTS

  ## Save patch-seq mode into taxonomy anndata
  AIT.anndata$write_h5ad(file.path(taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  
  return(AIT.anndata)
}



#' Save marker genes and QC metrics for patchSeqQC
#'
#' This is the same as addPatchseqQCMetrics, but the name buildPatchseqTaxonomy is kept for back-compatibility
#'
#' @param ... Variables to be passed to `buildPatchseqTaxonomy`
#' 
#' @return AIT.anndata An updated AIT.anndata variable.
#'
#' @export
buildPatchseqTaxonomy <- function(...) { addPatchseqQCMetrics(...) }
