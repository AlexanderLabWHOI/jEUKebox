# jEUKebox: Simulating metatranscriptomes for environmental eukaryotic microbes
_Pronounced jukebox_

## Protistan microbial communities: the forgotten validation in metatranscriptome assessment

Protistan communities are diverse. They are also not generally the focus of assessment and validation procedures for 'omics techniques. Bacteria lead the way as part of the "pack" of microbes that need to be assembled, identified, and otherwise validated. This is because of the relative simplicity (speaking loosely) of bacterial genomes and transcriptomes, but also because of the existing research that already supports prokaryoticgenetic assessment. This is particularly poignant due to the role of bacteria in human disease and health, which makes them a priority for study. 

However, metatranscriptomes (among other tools) are increasingly being applied to eukaryote-dominated environmental communities. Though validation has been performed using other genetic tools, e.g. PCR approaches, or through microscopic counts, there is an urgent need for general validation of omics tools for assembling and annotating eukaryote-dominated metatranscriptomes. 

For this reason, we have developed this companion tool to the `eukrhythmic` microbial eukaryote metatranscriptome analysis pipeline!

![the jEUKebox pipeline](files/jEUKebox-pipeline.png)

## How to run the pipeline

In order to run the pipeline, you need to have installed:

- Snakemake > 7
- mamba

The remainder of the dependent packages required will be downloaded for you during execution.

`jEUKebox` is written in [`Snakemake`](https://snakemake.readthedocs.io/en/stable/). The workflow can be executed just like an ordinary `Snakemake` workflow. The `Snakefile` is named `jEUKebox`, and lists all of the outputs expected from running the pipeline, calls upon methods written in the `rules` folder, which handles the execution of individual software tools needed to create the output files. 

To run `jEUKebox`, one would execute:

```
snakemake --cores <number-of-cores> -s jEUKebox
```

This would run the pipeline *locally*. In order to execute the pipeline on a high-performance computing cluster, users can modify the file in `submit/snake-submit.sh`. This is designed for a `SLURM` HPC manager, but can be adapted to meet the needs of other systems. 

## The configuration file

The configuration file allows the user to customize the pipeline to suit the particular collection of organisms of interest and their abundance in the communities.

- transcriptome_selections: files can be specified individually in list form; if this field is included in the configuration file, "protein_files" should have the same length.
- protein_selections: analogous protein files for the transcriptomes in the first configuration file field. Alternatively, both of these fields can be circumvented by providing a directory.
- community_spec: a CSV file containing important specification information for each of the generated communities. See the section below for more detailed information on how to set up this file.
- outputdir: the base directory where the output of the processing pipeline should be stored. 


## The `community_spec` file in the configuration

This file is the way to control the desired community composition of the generated _assemblies_. Note that this does not control the raw reads that might come from each of the contigs (its representation in the sample), which is done based on an exponential distribution for the various _transcripts_ which were simulated as having been sequenced, weighted by the abundance of each organism. 

There is a default csv file included in the `files` directory of the repository, but this can be created manually using the following column fields:

- Community - the number of each community to be included (should be distinct).
- NumberOrganisms - the number of _distinct_ organisms to be included in each community, from the nucleotide and protein files as specified in the configuration.
- ProportionCommunity - a list specifying the proportion of contigs within its community that each organism should represent in the simulated assembly.
- NumberHighSimilarity - a number (less than total NumberOrganisms for each community row) of the organisms within the community that should have high ANI similarity. The most similar organisms per ANI similarity will be chosen first. 
- GroupsHighSimilarity - the number of _groups_ with high similarity that are expected to be found. This allows for the case in which you wish to have multiple high-similarity groups within your data, but want those two groups to be relatively unrelated to one another. The NumberHighSimilarity value should be divisible by the GroupsHighSimilarity value, else the total number of organism in the community will be reduced to accommodate even divisibility into the group size.

In the future, we intend to create a mode where users can just explicitly specify the output community file, such that the intermediate step is bypassed. If this feature is of interest to you, please submit an issue on the **Issues** [tab](https://github.com/AlexanderLabWHOI/jEUKebox/issues) describing your use case, and we'll get back to you as soon as we can!

**Note**: if there is no preference for either or both of _NumberHighSimilarity_ or _GroupsHighSimilarity_, they can be set to -1 and will be ignored.