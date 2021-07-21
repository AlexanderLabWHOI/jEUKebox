# jEUKebox: Simulating metatranscriptomes for environmental eukaryotic microbes
_Pronounced jukebox_

## Protistan microbial communities: the forgotten validation in metatranscriptome assessment

Protistan communities are diverse. They are also not generally the focus of assessment and validation procedures for 'omics techniques. Bacteria lead the way as part of the "pack" of microbes that need to be assembled, identified, and otherwise validated. This is because of the relative simplicity (speaking loosely) of bacterial genomes and transcriptomes, but also because of the existing research that already supports prokaryoticgenetic assessment. This is particularly poignant due to the role of bacteria in human disease and health, which makes them a priority for study. 

However, metatranscriptomes (among other tools) are increasingly being applied to eukaryote-dominated environmental communities. Though validation has been performed using other genetic tools, e.g. PCR approaches, or through microscopic counts, there is an urgent need for general validation of omics tools for assembling and annotating eukaryote-dominated metatranscriptomes. 

For this reason, we have developed this companion tool to the `eukrhythmic` microbial eukaryote metatranscriptome analysis pipeline!

## The `community_spec` file in the configuration

This file is the way to control the desired community composition of the generated _assemblies_. Note that this does not control the raw reads that might come from each of the contigs (its representation in the sample), which is done based on an exponential distribution for the various _transcripts_ which were simulated as having been sequenced, weighted by the abundance of each organism. 

There is a default csv file included in the `files` directory of the repository, but this can be created manually using the following column fields:

- Community - the number of each community to be included (should be distinct).
- NumberOrganisms - the number of _distinct_ organisms to be included in each community, from the nucleotide and protein files as specified in the configuration.
- ProportionCommunity - a list specifying the proportion of contigs within its community that each organism should represent in the simulated assembly.
- NumberHighSimilarity - a number (less than total NumberOrganisms for each community row) of the organisms within the community that should have high ANI similarity. The most similar organisms per ANI similarity will be chosen first. 
- GroupsHighSimilarity - the number of _groups_ with high similarity that are expected to be found. This allows for the case in which you wish to have multiple high-similarity groups within your data, but want those two groups to be relatively unrelated to one another.

**Note**: if there is no preference for either or both of _NumberHighSimilarity_ or _GroupsHighSimilarity_, they can be set to -1 and will be ignored.