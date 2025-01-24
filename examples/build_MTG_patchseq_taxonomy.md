
# Tutorial: Building a Patch-seq Shiny taxonomy (start to finish) - Human MTG

In this example we demonstrate how to setup a Patch-seq Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. This tutorial parallels the other tutorial for building a Patch-seq Shiny taxonomy, but using the SEA-AD taxonomy based only on adult neurotypical donors [Gabitto, Travaglini et al taxonomy](https://www.nature.com/articles/s41593-024-01774-5) as reference and query Patch-seq data from [Berg et al 2020](https://www.nature.com/articles/s41586-021-03813-8). This example also includes a bit more explanatory text and does not assume you've created the taxonomy yet. If you are bringing your own data to the tutorial you can replace all sections that say "**FOR EXAMPLE ONLY**" with your own data munging steps. 

For creating a standard taxonomy and mapping against it, the following input variables are required: 

To create a taxonomy in AIT format:
* Reference count matrix (gene x cell), with genes as rownames and sample identifiers as colnames
* Reference annotation data.frame (cell x field), with sample identifiers as rownames (this should include cluster calls)
* Variable and/or marker genes (vector), can be calculated if not provided
* UMAP coordinates (cell x 2), can be calculated if not provided
* Dendrogram of clusters (dendrogram), optional (except for tree mapping), can be calculated if not provided

To map patch-seq (or other) data to the taxonomy:
* Query count matrix and metadata (of the format above)
* Query annotation data.frame (cell x field), optional, will be passed through to MolGen Shiny directory

### Before starting, load the necessary libraries

```R
## Load libraries
library(scrattch.taxonomy)
library(scrattch.mapping)
library(scrattch.patchseq)
library(reticulate) # For hierarchical mapping
cell_type_mapper <- import("cell_type_mapper") # For hierarchical mapping

## Specify which reference taxonomy to map against.
## -- Replace folder and file name with correct location
taxonomyDir = getwd() 
```

## Part 1: Building the taxonomy

This section describes how to build a reference taxonomy for mapping. **If you have already built a taxonomy in AIT format, you can skip to Part 2: Mapping to the taxonomy.**

### 1.1: Read in the reference data

The first step is to read in your cell (columns) by gene (rows) data matrix of counts alongside the metadata associated with each cell.  These are saved as taxonomy.counts and taxonomy.metadata, respectively in the section below.  We use data from Gabitto, Travaglini et al 2024 **FOR EXAMPLE ONLY**.  Note that colnames(taxonomy.counts) and rownames(taxonomy.metadata) should be identical and correspond to cell IDs, and rownames(taxonomy.counts) should be the gene names (in the is case, gene symbols).  It is also okay to swap the rows and columns the in the cellxgene matrix, so long as the correponding row and column names still match.

```R
## Download the reference data to the working directory and read it in
seaad_url  <- "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
dend_url   <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0f/37/0f3755cb-3acb-4b93-8a62-5d6adc74c673/dend.rds"
#download.file(seaad_url,"Reference_MTG_data.h5ad")  # NOTE: we recommend downloading via the web browser, as this command may fail
#download.file(dend_url,"Reference_MTG_dend.rds")
seaad_data <- read_h5ad("Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
seaad_dend <- readRDS("Reference_MTG_dend.rds")

## Subsample data (this can be done either here, or within buildTaxonomy)
keepCells <- subsampleCells(seaad_data$obs$cluster_label,100,seed=42)

## Select only a few clusters to speed up calculations (NOTE: YOU WOULD NORMALLY NEVER DO THIS)
#keepCl    <- c(paste0("Lamp5_",1:3),paste0("Sst_",4:5),paste0("Pvalb_",6:7),
#               paste0("L2/3 IT_",1:2),paste0("L6b_",1:2),paste0("L5 ET_",1:2),
#			          paste0("Astro_",1:2),paste0("Oligo_",1:2),paste0("Micro-PVM_",1:2))
#seaad_dend<- prune(seaad_dend,setdiff(labels(seaad_dend),keepCl))
#keepCells <- keepCells&is.element(seaad_data$obs$cluster_label,keepCl)

## Get (subsampled) subset data and annotations
taxonomy.counts = (seaad_data$X)[keepCells,]
cn <- c("sample_name","cluster_label","cluster_confidence","subclass_label","class_label",
        "external_donor_name_label","age_label","donor_sex_label")
taxonomy.metadata = seaad_data$obs[keepCells,cn]

## Ensure count matrix and annotations are in the same order (this shouldn't be needed)
taxonomy.metadata = taxonomy.metadata[match(rownames(taxonomy.counts), taxonomy.metadata$sample_name),]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))

## Transpose the counts matrix (... for now; code in process to avoid transposing large matrices)
taxonomy.counts <- t(taxonomy.counts)
taxonomy.counts <- as(taxonomy.counts, "dgCMatrix")
```

### 1.2: Create the (parent) AIT Taxonomy 

This section will create the parent taxonomy for the reference data.  In this case, we include up 1000 cells for **every** cell type defined in Hodge et al 2019, along with their associated metadata, and will subsample the clusters and cells further at a later step.

This code block loads the scrattch.taxonomy library, and then calculates variables genes and defines a UMAP **using a very basic approach**.  If variable genes and/or 2-dimensional coordinates already exist, they can be provided to buildTaxonomy below rather than calculated in this way. 

```R
## Compute top 1000 binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.metadata$cluster, 1000)

## Compute UMAP coordinates
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames (Required!)
rownames(umap.coords) = colnames(taxonomy.counts)
```

The next step builds the parent taxonomy using a single call to the function buildTaxonomy.  After running this script, your taxonomy will contain all of the data and metadata in standard formats, and will be ready for correlation and tree mapping.  However, **you still need to run buildPatchseqMapping in the next section** for tree mapping, subsetting, and other QC metrics.  

```R
## Set up the levels of hierarchy for all mapping functions later
## -- This MUST be from broadest to most specific types, and NOT vice versa
hierarchy = list("class_label", "subclass_label", "cluster_label")

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(
      counts        = taxonomy.counts,
      meta.data     = taxonomy.metadata,
      dend          = seaad_dend,  # If this is omitted buildTaxonomy will generate a dendrogram
      feature.set   = binary.genes,
      umap.coords   = umap.coords,
      taxonomyTitle = "AI_taxonomy",  # Determines the file name
      taxonomyDir   = taxonomyDir,
      subsample     = 100, # A lot of subsampling to speed up calculations
      hierarchy     = hierarchy
)

## If you also want to create a shiny directory for the REFERENCE data set, uncomment this section
## Create Shiny directory (AIBS-internal)
#createShiny(AIT.anndata,
#            shinyDir = getwd(),  # Replace location with desired location for shiny directory output
#            metadata_names = NULL)

## Alternatively, if you have already created the taxonomy, you can load it using "loadTaxonomy"
## Load the taxonomy (from h5ad file name)
#AIT.anndata = loadTaxonomy(taxonomyDir, "AI_taxonomy.h5ad")
```

### 1.3: Create a child taxonomy for patch-seq mapping

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.taxonomy::mappingMode()` as discussed next.**

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(
                 AIT.anndata,
                 mode.name = "patchseq", ## Give a name to off.target filtered taxonomy
                 subsample = 100, ## Subsampling for the new taxonomy.
                 subclass.column = "subclass_label", 
                 class.column = "class_label", ## The column by which off-target types are determined.
                 off.target.types = c("Non-neuronal and Non-neural"), ## The off-target class.column labels for patchseqQC.
                 subclass.subsample = 100, ## Subsampling is for PatchseqQC contamination calculation.
                 num.markers = 50, ## Number of markers for each annotation in `class_label`
                 taxonomyDir = taxonomyDir ## Replace with location to store taxonomy
) 
```

The `buildPatchseqTaxonomy` function does the following, updating the anndata variable and file accordingly:
* Creates a "child" taxonomy that only includes the cells and clusters requested (e.g., neurons for patch-seq mapping)
* Defines marker genes for each node of the dendrogram for use with tree mapping
* Creates a subsetted version of the parent dendrogram that also includes node marker genes
* Creates the required marker and expression variables for 'QC_markers' for use with patchseqQC
* Creates a table of cell to cluster probability (e.g., 'membership') values for calculation of KL divergence and creation of constellation diagrams
Currently a subfolder for the child taxonomy is also created, but everything in that folder is also stored in the anndata, and so it can be safely ignored.

**At this point the reference taxonomy is created and ready for Patch-seq mapping.**


## Part 2: Mapping to the taxonomy and patch-seq quality control

The rest of this example demonstrates how to map patch-seq data to the reference defined above and apply quality control metrics to patch-seq data.  These are the primary use cases of the `scrattch.mapping` and `scrattch.patchseq` libraries, respectively. This part is entirely distinct from creation of a reference taxonomy in Part 1 above. 


### 2.1: Read in the QUERY (e.g., Patch-seq) data

The first step is to read in your cell (columns) by gene (rows) data matrix of **log-normalized** query counts alongside (optional) metadata associated with each cell.  These are saved as query.logCPM and query.anno, respectively in the section below.  We use data from Berg et al 2020 **FOR EXAMPLE ONLY**.  Note that colnames(query.logCPM) and rownames(query.anno) should be identical and correspond to cell IDs, and rownames(query.logCPM) should be the gene names (in the is case, gene symbols).  **Only common gene names in query.logCPM and taxonomy.counts will be used for mapping.**

```R
## Download data and metadata from GitHub repo for Berg et al 2022
download.file("https://github.com/AllenInstitute/patchseq_human_L23/raw/master/data/input_patchseq_data_sets.RData", "patchseq.RData", mode="wb")

## Load the data
load("patchseq.RData")

## Rename the query data and metadata for convenience
query.anno = annoPatch  # Some cell annotations for all cells from the paper 
query.logCPM = datPatch # logCPM values for all cells from the paper
```


### 2.2: Set scrattch.mapping mode

**Do not skip this step!** Here we set taxonomy mode to the relevant "child" taxonomy defined in 1.3 above.  This will let the mapping functions know to  only consider the unfiltered cells and cell types (e.g., use subsetted cells from neuronal clusters for Patch-seq mapping).

```R
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
```

### 2.3: Map query cells to Patch-seq reference

This is the step that performs the mapping.  Since we have very different query and reference data sets, Seurat mapping in this example may not be reliable. Note that each mapping algorithm can map to multiple levels of the taxonomy. 

```R
# This function is part of the 'scrattch.mapping' library
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                 query.data = query.logCPM,
                                 label.cols = hierarchy,  # Will default to AIT.anndata$uns$hierarchy if not provided
                                 corr.map   = TRUE, # Flags for which mapping algorithms to run
                                 tree.map   = TRUE, 
                                 seurat.map = TRUE,
                                 hierarchical.map = TRUE)

## If you want the mapping data.frame and associated scores from the S4 mappingClass
mapping.results = getMappingResults(query.mapping, scores=TRUE)
```

### 2.4: QC patch-seq data and create shiny directory

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
    query.data     = query.logCPM,  ## Counts are required here (NOT cpm or logCPM), but it will convert to linear space automatically
    query.metadata = query.anno,
    query.mapping  = query.mapping, ## This has to be an S4 mappingClass from scrattch.mapping.
    doPatchseqQC   = TRUE,  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
    return.metrics = TRUE  ## Set to TRUE to return the updated metrics table
)
```
