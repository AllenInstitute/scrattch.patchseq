#' Save marker genes for patchSeqQC
#'
#' This function saves all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm) to AIT.anndata$uns. This is only used for patch-seq analysis.  Requirements for input include:
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param variable.genes Set of variable genes used for correlation and Seurat mapping. We recommend providing this gene set, but if not provided the gene set will default to the highly.variable.genes from mode="standard" (not recommended).  This can also be set after the fact using `updateHighlyVariableGenes`.
#' @param subsample The number of cells to retain per cluster (default = 100). See details.
#' @param subclass.column Column name corresponding to the moderate-resolution cell types used for the cell types of interest (default = "subclass_label"). See details.
#' @param class.column Column name corresponding to the low-resolution cell types used for the off-target cell types (default = "class_label"). See details.
#' @param off.target.types A character vector of off-target (also known as 'contamination') cell types.  This must include at least one of the cell types found in "class.column" and/or "subclass.column" (both columns are checked)
#' @param mode.name A name to identify the new taxonomy version.
#' @param subclass.subsample The number of cells to retain for PatchseqQC contamination calculation (default = 100, probably no need to change).
#' @param num.markers The maximum number of markers to calculate per node per direction (default = 50)
#' @param taxonomyDir The location to save shiny output (default = current working directory).
#' @param addMapMyCells If TRUE (default), will also prep this mode of the taxonomy for hierarchical mapping
#' @param hierarchy Column names of annotations to map against with hierarchical mapping.  Note that this only works for metadata that represent clusters or groups of clusters (e.g., subclass, supertype, neighborhood, class). Will default to whatever is in AIT.anndata$uns$hierarchy and is required if addMapMyCells=TRUE
#' @param ... Additional variables to be passed to `addDendrogramMarkers`
#' 
#' **How filtering works**
#' 
#' Currently taxonomies are subsetted (or 'filtered') using two methods.  First, at most a random set of `subsample` cells are included per cluster (and therefore the rest of the cells per cluster are filtered).  If you do not want to subsample, then set subsample=Inf.  Second, off target clusters are removed. More specifically, any cells with values of `off.target.types` in either `subclass.column` or `class.column` are removed. Note that cells labeled as TRUE in the filter are the ones REMOVED.
#' 
#' **Notes on PatchseqQC**
#' 
#' PatchseqQC requires cell types defined at two resolutions to work, currently defined in the class `subclass.column` and `class.column` for the higher and lower-resolution cell types respectively. PatchseqQC has the concept of on target vs. off target gene expression.  For on target types, it uses the lower resolution types and for off-target types it uses the higher-resolution types. QC metrics can then estimate the quality of the cell along with how much contamination there is from other types of cells. PatchseqQC isn't called in this function, but several of the additions listed below are input files required to call it later.
#' 
#' **Additions to AIT.anndata**
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
#' $memb[[mode.name]]  
#' ...$memb.ref,  
#' ...$map.df.ref
#' 
#' The following variable is added to AIT.anndata$var
#' $highly_variable_genes_{AIT.anndata$uns$mode}
#' 
#' @import patchseqtools
#' @import scrattch.hicat
#' @import reticulate
#'
#' @return AIT.anndata An updated AIT.anndata variable with the above content added to AIT.anndata$uns for the relevant mode.name.
#'
#' @export
buildPatchseqTaxonomy = function(AIT.anndata,
                                 variable.genes = NULL,
                                 mode.name = "patchseq", 
                                 subsample = 100,
                                 subclass.column = "subclass_label",
                                 class.column = "class_label",
                                 off.target.types = c("Glia","glia","non-neuronal","Non-neuronal"), ## "Glia", "NN"
                                 subclass.subsample = 100,
                                 num.markers = 50,
                                 taxonomyDir = file.path(AIT.anndata$uns$taxonomyDir),
                                 addMapMyCells = TRUE,
                                 hierarchy = AIT.anndata$uns$hierarchy,
                                 ...
){

  ## Capture current mode so after output the mode is unchanged
  currentMode = AIT.anndata$uns$mode
  
  ## Ensure filtering mode doesn't already exist
  if(mode.name %in% names(AIT.anndata$uns$filter)){ print(paste0("Mode ", mode.name, " already in Taxonomy, you will be overwriting the previous mode files.")) }

  ## Create the required files for patchSeqQC and determine offtarget cells
  if(!is.element(subclass.column, colnames(AIT.anndata$obs))){stop(paste(subclass.column,"is not a column in the metadata data frame."))}
  if(!is.element(class.column, colnames(AIT.anndata$obs))){stop(paste(class.column,"is not a column in the metadata data frame."))}
  if(!dir.exists(file.path(taxonomyDir))){"Specified taxonomy folder does not exist."}

  ## Determine taxonomy mode directory (Move to utility function)
  if(mode.name == "standard"){ taxonomyModeDir = file.path(taxonomyDir) } else { taxonomyModeDir = file.path(file.path(taxonomyDir), mode.name) }
  if(!dir.exists(taxonomyModeDir)){ dir.create(taxonomyModeDir, showWarnings = TRUE) }

  ## Copy metadata
  metadata = AIT.anndata$obs

  ## Ensure variable naming scheme matches assumptions  
  # --- NOTE: we should consider renaming these variables.
  metadata$subclass_label = AIT.anndata$obs[,subclass.column]  # For compatibility with existing code.
  metadata$class_label = AIT.anndata$obs[,class.column]  # For compatibility with existing code.
  
  ## Subsample and filter metadata and data
  kpSamp2  = subsampleCells(metadata$subclass_label, subclass.subsample)
  goodSamp = !is.na(metadata$class_label)  # For back-compatibility; usually not used
  kpSamp2  = kpSamp2 & goodSamp            # For back-compatibility; usually not used
  annoQC   = metadata[kpSamp2,]
  annoQC$subclass_label = make.names(annoQC$subclass_label)
  datQC    = as.matrix(Matrix::t(AIT.anndata$raw$X))[,kpSamp2]
  colnames(datQC) = AIT.anndata$obs_names[kpSamp2]; rownames(datQC) = AIT.anndata$var_names

  ## Define class and subclass labels
  ## --- We wrap on-target types by class but retain off-target types by subclass
  offTarget = is.element(annoQC$class_label, off.target.types) | is.element(annoQC$subclass_label, off.target.types)
  if(sum(offTarget)==0){stop("No valid off-target classes or subclasses are provided. Please update off.target.types accordingly.")}
  
  ## 
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

  ## Filter out off target cells along with additional cells beyond those subsampled
  AIT.anndata$uns$filter[[mode.name]] = is.element(metadata$class_label, off.target.types) | is.element(metadata$subclass_label, off.target.types)
  AIT.anndata$uns$filter[[mode.name]] = !((!AIT.anndata$uns$filter[[mode.name]]) & ((subsampleCells(metadata$cluster_label,subsample)))) # NEW, for subsampling

  ## Filter stats files if they exist
  cluster.keep <- is.element(AIT.anndata$uns$stats[["standard"]][["clusters"]],AIT.anndata$obs$cluster_label[!AIT.anndata$uns$filter[[mode.name]]])
  AIT.anndata$uns$stats[[mode.name]][["medianExpr"]] = AIT.anndata$uns$stats[["standard"]][["medianExpr"]][,cluster.keep]
  AIT.anndata$uns$stats[[mode.name]][["features"]] = AIT.anndata$uns$stats[["standard"]]$features
  AIT.anndata$uns$stats[[mode.name]][["clusters"]] =  AIT.anndata$uns$stats[["standard"]][["clusters"]][cluster.keep]

  ## Save patchseqQC information to uns
  AIT.anndata$uns$QC_markers[[mode.name]] = list("allMarkers" = allMarkers,
                                                  "markers" = markers,
                                                  "countsQC" = countsQC,
                                                  "cpmQC" = cpmQC,
                                                  "classBr" = classBr,
                                                  "subclassF" = subclassF,
                                                  "qc_samples" = colnames(countsQC),
                                                  "qc_genes" = rownames(countsQC))
  

  ##################
  ## ------- Modify the dendrogram and save
  ##

  ## Load the complete dendrogram, always from standard mode
  dend = json_to_dend(AIT.anndata$uns$dend[["standard"]])

  ## Prune dendrogram to remove off.target types
  dend = prune(dend, setdiff(labels(dend), unique(AIT.anndata$obs$cluster_label[!AIT.anndata$uns$filter[[mode.name]]])))

  ## Save dendrogram
  saveRDS(dend, file.path(taxonomyModeDir, "dend.RData"))

  ## Store the pruned dendrogram, in dend list under "patchseq" mode.name
  AIT.anndata$uns$dend[[mode.name]] = toJSON(dend_to_json(dend))

  ## Save patch-seq mode into taxonomy anndata
  AIT.anndata$write_h5ad(file.path(taxonomyDir, paste0(AIT.anndata$uns$title, ".h5ad")))
  
  ## Update the log file and check the taxonomy for proper quality
  if(!checkTaxonomy(AIT.anndata,taxonomyDir)){
    stop(paste("Taxonomy has some breaking issues.  Please check checkTaxonomy_log.txt in",taxonomyDir,"for details"))
  }

  ## Update markers after pruning
  AIT.anndata = addDendrogramMarkers(AIT.anndata, mode=mode.name, taxonomyDir=taxonomyDir, ...)
  # The reference probability matrix for the subsetted taxonomy is defined and outputted in this function as well
  # $memb[[mode.name]]
  # ...$memb.ref,
  # ...$map.df.ref
  
  ## Add MapMyCells (hierarchical mapping) functionality, if requested
  if(addMapMyCells) {
    if(is.character(hierarchy)) hierarchy <- as.list(hierarchy)
    if((length(hierarchy)==0)|(sum(class(hierarchy)=="list")<1)){
      warning("hierarchy must be a list of term_set_labels in the reference taxonomy ordered from most gross to most fine included in AIT_anndata or provided separately. Since this is NOT the case, addMapMyCells is being skipped")
    } else{
      AIT.anndata = mappingMode(AIT.anndata, mode=mode.name)
      AIT.anndata = addMapMyCells(AIT.anndata, hierarchy, force=TRUE)
      AIT.anndata = mappingMode(AIT.anndata, mode=currentMode)
    }
  }
  
  ## Add highly variable genes. If not provided, default to global highly variable genes.
  AIT.anndata = updateHighlyVariableGenes(AIT.anndata, variable.genes, mode.name)
  AIT.anndata = mappingMode(AIT.anndata, mode=currentMode)
  
  ##
  return(AIT.anndata)
}
