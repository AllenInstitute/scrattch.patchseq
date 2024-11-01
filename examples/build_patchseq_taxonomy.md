# Tutorial: Building a patchseq Shiny taxonomy 

In this tutorial we demonstrate how to setup a patchseq Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. 

### Required inputs:

* Standard Shiny taxonomy setup following the "build_taxonomy" tutorial.
* Query patchseq count matrix and metadata.

### Additional prerequisites:

* Installation of the `tasic2016data` data package [from here](https://github.com/AllenInstitute/tasic2016data/).
* Installation of `scrattch.mapping` [from here](https://github.com/AllenInstitute/scrattch.mapping) for data mapping. 

### Load in test data from Tasic 2016:
```R
## Load libraries
library(scrattch.taxonomy)
library(scrattch.mapping)
library(scrattch.patchseq)
library(tasic2016data)

## Load in the tasic2016 data and wrangle as a query data set.
query.anno = tasic_2016_anno
query.counts = tasic_2016_counts 
query.anno = query.anno[match(colnames(query.counts),query.anno$sample_name),]
rownames(query.anno) = query.anno$sample_name  
keep = query.anno$broad_type!="Unclassified"
query.counts = query.counts[,keep]
query.logCPM = logCPM(query.counts)
query.anno = query.anno[keep,]
```

### Load and setup the base taxonomy:
```R
## Standard shiny taxonomy
# NOTE: replace 'taxonomy' location below with output folder from the "build_taxonomy" tutorial
taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016"

## Load in the taxonomy
AIT.anndata = loadTaxonomy(taxonomy, anndata_file="Tasic2016.h5ad")
```

### Define off target cell types

Let's start by defining the cell classes which are off target for tasic2016 patchseq mapping and patchseqQC. For neocortex this is typically defined at the "class" level and are the "Non-neuronal" cell classes/types.
```R
## Identify the offtarget cell types manually.
print(unique(AIT.anndata$obs$broad_type_label))

## Add in the off.target annotation.
AIT.anndata$obs$off_target = AIT.anndata$obs$broad_type_label

## First we need to add our off.target annotation to the factor levels
levels(AIT.anndata$obs$off_target) = c(levels(AIT.anndata$obs$off_target), "Non-neuronal")

## Now lets set all Non-neuronal cells to the "Non-neuronal" off target annotation.
AIT.anndata$obs$off_target[!is.element(AIT.anndata$obs$off_target, c("GABA-ergic Neuron","Glutamatergic Neuron", "Astrocyte"))] = "Non-neuronal"
```

### Build the patchseq taxonomy:

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.taxonomy::mappingMode()` as dicusssed next.**

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "cluster_label", ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                    class.column = "off_target", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016")
```
The `buildPatchseqTaxonomy` function return/created the following:

* An updated AIT.anndata object for patchseq mapping and QC steps.
* Created the required marker and expression files for patchseqQC and save under 'mode.name' sub-directory

### Set scrattch.mapping mode

Now we will set scrattch.mapping to use only cells not in the off.target.types, this will filter the taxonomy and adjust the dendrogram to remove any `cluster` in `off.target.types`.

```R
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
```

### Map against the patchseq taxonomy:
```R
# This function is part of the 'scrattch.mapping' library
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.logCPM,
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = FALSE, 
                                  hierarchical.map = FALSE, ## Not that this will not work correctly with the patchseq mode. TODO.
                                  label.cols = c("cluster_label", "broad_type_label")) # Columns to map against from AIT.anndata$obs
## If you want the mapping data.frame from the S4 mappingClass
mapping.results = getMappingResults(query.mapping, scores=FALSE)
```

### Determine patchseq contamination with PatchseqQC:
```R
patchseq.qc = applyPatchseqQC(AIT.anndata, ## A patchseq taxonomy object.
                                query.counts, ## Counts are required here.
                                query.anno, ## Query annotations to add PatchSeqQC onto.
                                verbose=FALSE)
```

### Setup the patchseq Shiny taxonomy files for MolGen Shiny:
```R
buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = '/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016/patchseq_mapping',
                      query.data     = query.counts, ## Counts are required here.
                      query.metadata = query.anno,
                      query.mapping  = query.mapping, ## This has to be an S4 mappingClass from scrattch.mapping.
                      doPatchseqQC   = FALSE)  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
```
