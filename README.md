# scrattch.patchseq

Generalized scripts for patchseq analysis at the Allen Institute.

## Documentation

You can find a detail description of all scrattch.patchseq functions here: ![Documentation](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/scrattch.patchseq_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/VERSIONS.md)

## Docker

We have setup a docker environemnt for scrattch.taxonomy and scrattch.mapping that contains all the required dependencies and the current version of all scrattch packages. This docker is accessible through docker hub via: `njjai/scrattch_mapping:0.6.6`.

#### HPC usage:

##### Non-interactive
`singularity exec --cleanenv docker://njjai/scrattch_mapping:0.6.6 Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.6`


## Usage examples

1. [**Build and map against a small mouse PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/examples/build_patchseq_taxonomy.md) This example provides the basics for updating a taxonomy to be compatible with patchseq style mapping and visualization on (internal Allen Institute) MolGen Shiny tools.

2. [**Build and map against a human MTG PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch.patchseq/blob/main/examples/build_MTG_patchseq_taxonomy.md) This example provides shows how to create a standard taxonomy and update it to be compatible with patchseq style mapping and visualization on (internal Allen Institute) MolGen Shiny tools. This example essentially combines examples 1 and 2 and applies them to human neocortical data sets.  Data is from Hodge et al. (2019) for human MTG and patch-seq examples are from Berg et al (2021) (layer 2-3, excitatory neurons). 

   
## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-patchseq/issues).

## TODO

- [ ] Update documentation.

## Done
