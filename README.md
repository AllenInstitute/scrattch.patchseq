# scrattch.patchseq

Generalized scripts for patchseq analysis at the Allen Institute.

## Documentation

You can find a detail description of all scrattch.patchseq functions here: ![Documentation](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/scrattch.patchseq_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/VERSIONS.md)

## Installation

### Using docker (recommended)
We have setup a docker environment for scrattch.taxonomy, scrattch.mapping, and scrattch.patchseq that contains all the required dependencies and the current version of all scrattch packages. **See [the readme](https://github.com/AllenInstitute/scrattch/blob/master/README.md#using-docker) for [the parent scrattch package](https://github.com/AllenInstitute/scrattch) for the most up-to-date docker information.**

### Directly from GitHub (strongly discouraged)

While we advise using the provided docker, you can also install scrattch.patchseq directly from GitHub as follows:

```
devtools::install_github("AllenInstitute/scrattch.patchseq")
```

This strategy **might not work** due to complicated dependencies. Also note that `doMC` may need to be installed manually from [HERE](https://r-forge.r-project.org/R/?group_id=947) if you use Windows. Vignettes are provided below.

## Usage examples

1. [**Map against a small mouse PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/examples/build_patchseq_taxonomy.md) This example provides the basics for updating a taxonomy to be compatible with patchseq style mapping and visualization on (internal Allen Institute) MolGen Shiny tools, and for collecting quality control metrics of potential use to everyone. This is a continuation of examples in scrattch.taxonomy and scrattch.mapping building data from [Tasic et al 2016](https://www.nature.com/articles/nn.4216) as reference and data from [Gouwens, Sorensen, et al 2020](https://www.sciencedirect.com/science/article/pii/S009286742031254X) as query.

2. [**Build and map against a human MTG PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/examples/build_MTG_patchseq_taxonomy.md) This example provides shows how to create a standard taxonomy, update it to be compatible with patchseq style mapping and visualization on (internal Allen Institute) MolGen Shiny tools, and for collecting quality control metrics of potential use to everyone. This example essentially combines examples 1 and 2 and applies them to human neocortical data sets.  Reference data is from the [Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)](https://portal.brain-map.org/explore/seattle-alzheimers-disease) project [(Gabitto, Travaglini et al., 2024)](https://www.nature.com/articles/s41593-024-01774-5) for human middle temporal gyrus (MTG) and patch-seq query data include layer 2-3, excitatory neurons from [Berg et al (2021)](https://www.nature.com/articles/s41586-021-03813-8). 

   
## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-patchseq/issues).

## TODO

- [ ] Update documentation.

## Done
