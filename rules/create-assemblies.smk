configfile: "config.yaml"

import io
import os
from os import listdir
from os.path import isfile, join
import sys
import pandas as pd
import numpy as np
import pathlib
import random
from Bio import SeqIO
import datetime
from snakemake.exceptions import print_exception, WorkflowError  

mydate = datetime.datetime.now()
datestring = mydate.strftime("%b%d")

rule create_assemblies:
    input:
        community_files_prot = os.path.join(config["outputdir"], "community_{comm}_prot", "{comm}_complete.csv"),
        community_files = os.path.join(config["outputdir"], "community_{comm}", "{comm}_complete.csv"),
        communities_file = os.path.join(config["outputdir"], "communities.csv"),
        orthfinder_genct = os.path.join(config["outputdir"],"orthofinder","orthofinder_{comm}",
                                        "Results_Jul17","Orthogroups","Orthogroups.GeneCount.tsv"),
        orthfinder_OGs = os.path.join(config["outputdir"],"orthofinder","orthofinder_{comm}",
                                        "Results_Jul17","Orthogroups","Orthogroups.tsv")
    output:
        mock_assembly = os.path.join(config["outputdir"], "mock_assemblies",
                                                 "mock_assembly_{comm}",
                                                 "mock_assembly_{comm}.fasta")
    params:
        community_dir = os.path.join(config["outputdir"], "community_{comm}")
    run:
        community_spec = pd.read_csv(input.communities_file)
        files_inassembly = os.listdir(params.community_dir)
        community_spec = community_spec[community_spec.MMETSP_inds.isin(files_inassembly)] 
                                            #community_spec.MMETSP_inds in files_inassembly,:]
                                            
        # we want to use the OrthoFinder output to ensure that each assembly contains contigs from
        # OG spanned by ALL species in assembly.
        
        gene_ct = pd.read_csv(input.orthfinder_genct,sep="\t")
        orthogroups = pd.read_csv(input.orthfinder_OGs,sep="\t")
        gene_ct = gene_ct.replace(0, np.nan)
        gene_ct = gene_ct.dropna(how='any', axis=0)
        to_write=[]
        for row_ind in range(len(community_spec.index)):
            #record_dict = SeqIO.to_dict(SeqIO.parse(list(community_spec.Organism)[row_ind], "fasta"))
            record_list = list(SeqIO.parse(list(community_spec.Organism)[row_ind], "fasta"))
            percentage = list(community_spec.Proportion)[row_ind]
            total_items = len(record_list) #len(record_dict.items())
            #random_reads = random.sample(record_dict.items(), int(total_items * percentage))
            random_num = list(np.random.randint(0, high=total_items, size=int(total_items * percentage)))
            random_num = random.choice(list(range(0,total_items)), size=int(total_items * percentage), replace=False)
            to_write.extend([record_list[random_num_curr] for random_num_curr in random_num])
        with open(output.mock_assembly, 'w') as handle:
            SeqIO.write(to_write, handle, 'fasta')