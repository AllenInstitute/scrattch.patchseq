# Tutorial: Building a patchseq Shiny taxonomy 

In this tutorial we demonstrate how to setup a Patchseq Shiny taxonomy using scrattch.patchseq for viewing on MolGen Shiny (AIBS internal) and running mapping algorithms against. This tutorial assumes you've already created a reference taxonomy from the [Tasic et al. 2016 study](https://www.nature.com/articles/nn.4216) by following the [scrattch.taxonomy example](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_taxonomy.md). Query data from the [Gouwens, Sorensen, et al 2020 study](https://doi.org/10.1016/j.cell.2020.09.057) on GABAergic interneurons in primary visual cortex is used as an example. 

*We strongly encourage running this code within the scrattch docker environment.  This example was created using docker://jeremyinseattle/scrattch:1.2 and will likely fail if run using any earlier scrattch versions.*


## Load libraries and reference taxonomy
```R
## Load libraries
suppressPackageStartupMessages({
  library(scrattch.taxonomy)
  library(scrattch.mapping)
  library(scrattch.patchseq)
  library(reticulate)
})
cell_type_mapper <- import("cell_type_mapper") # For hierarchical mapping

## Specify which reference taxonomy to map against.
## -- Replace folder and file name with correct locations
taxonomyDir = getwd() 
AIT.anndata = loadTaxonomy(taxonomyDir = taxonomyDir, anndata_file="Tasic2016.h5ad")
```

### Download, read in, and format query data
```R
## Download some example patch-seq data to map
## -- These data are from Gouwens et al 2020 and would be replaced by your query data
patchFolder  <- "https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/"
counts_url   <- "20200513_Mouse_PatchSeq_Release_count.csv.tar"
download.file(paste0(patchFolder,counts_url),counts_url)
untar(counts_url)
query.counts <- as.matrix(data.table::fread("20200513_Mouse_PatchSeq_Release_count/20200513_Mouse_PatchSeq_Release_count.csv"),rownames=1)
query.logCPM <- logCPM(query.counts)

## Download metadata from Gouwens et al 2020 (to match data from scrattch.mapping tutorial) 
## -- These data and metadata would be replaced by your query data
metadata_url <- "20200625_patchseq_metadata_mouse.csv.tar"
download.file(paste0(patchFolder,metadata_url),metadata_url)
untar(metadata_url)
query.anno   <- as.data.frame(data.table::fread("20200625_patchseq_metadata_mouse/20200625_patchseq_metadata_mouse.csv"))
rownames(query.anno) <- query.anno$transcriptomics_sample_id
query.anno   <- query.anno[colnames(query.counts),]
```

### Define off target cell types

Let's start by defining the cell classes which are off target for tasic2016 patchseq mapping and patchseqQC. For neocortex this is typically defined at the "class" level and are the "Non-neuronal" cell classes/types.
```R
## Identify the offtarget cell types manually.
print(unique(AIT.anndata$obs$broad_type))

## Add in the off.target annotation.
AIT.anndata$obs$off_target = AIT.anndata$obs$broad_type

## First we need to add our off.target annotation to the factor levels
levels(AIT.anndata$obs$off_target) = c(levels(AIT.anndata$obs$off_target), "Non-neuronal")

## Now lets set all Non-neuronal cells to the "Non-neuronal" off target annotation.
AIT.anndata$obs$off_target[!is.element(AIT.anndata$obs$off_target, c("GABA-ergic Neuron","Glutamatergic Neuron"))] = "Non-neuronal"
AIT.anndata$obs$off_target <- droplevels(AIT.anndata$obs$off_target)
```

### Build the patchseq taxonomy

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.taxonomy::mappingMode()` as dicusssed next.**  Note that we also need to explicitly set the mapping mode prior to calculating MapMyCells statistics and prior to mapping to this new taxonomy mode. 

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "primary_type_label", ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                    class.column = "off_target", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = taxonomyDir)  ## Replace with location to store taxonomy
```
The `buildPatchseqTaxonomy` function does the following:

* Create a new mode of the taxonomy corresponding to a subset of clusters in the taxonomy, with baggage needed for all mapping algorithms
* Returns updated AIT.anndata object to allow for patchseq mapping as well as for the QC steps below.

### Map against the patchseq taxonomy
```R
## Set scrattch.mapping mode - DO NOT SKIP THIS STEP
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")

# This function is part of the 'scrattch.mapping' library
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                 query.data = query.logCPM,
                                 genes.to.use="marker_genes_binary",
                                 corr.map   = TRUE, # Flags for which mapping algorithms to run
                                 tree.map   = TRUE, 
                                 seurat.map = TRUE, 
                                 mapmycells.hierarchical.map = TRUE,
                                 mapmycells.flat.map = TRUE)

## If you want the mapping data.frame and associated scores from the S4 mappingClass
mapping.results = getMappingResults(query.mapping, scores=TRUE)
```

### QC patch-seq data and create shiny directory

This section performs additional query data QC and then outputs the files necessary for visualization of Patch-seq data with molgen-shiny tools (AIBS internal) into a folder. If `return.metrics = TRUE`, this function also returns an updated annotation table that includes original metadata, mapping results, and all of the additional metrics and statistics described below.

More specifically, the function buildMappingDirectory:
* Creates a new folder to deposit all relevant query files (this is what you copy into the the molgen-shiny directory)
* Performs patchSeqQC on the query data to define Normalized Marker Sum (NMS) and other QC scores (see https://github.com/PavlidisLab/patchSeqQC for details)
* Calculates KL Divergence and associated Core/I1/I2/I3/PoorQ calls used to help assess Patch-seq quality
* Calculates and outputs tree mapping probabilities of every cell mapping to every cluster and tree node
* Calculates and outputs UMAP coordinates for query cells integrated into reference UMAP space
* Optionally returns the updated metadata file

```R
updated.query.anno <- buildMappingDirectory(
    AIT.anndata    = AIT.anndata, 
    mappingFolder  = file.path(taxonomyDir,"patchseq"),  ## Put the correct file path for output here
    query.data     = query.counts, ## Counts are required here (NOT cpm or logCPM)
    query.metadata = query.anno,
    query.mapping  = query.mapping, ## This has to be an S4 mappingClass from scrattch.mapping.
    genes.to.use   = "marker_genes_binary",
    doPatchseqQC   = TRUE,  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
    return.metrics = TRUE  ## Set to TRUE to return the updated metrics table
)
```
