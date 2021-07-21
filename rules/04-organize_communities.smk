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
from snakemake.exceptions import print_exception, WorkflowError  

# put the FASTA files of interest in their appropriate locations
rule organize_folders:
    input:
        communities_file = os.path.join(config["outputdir"], "03-community_spec", "communities.csv")
    output:
        community_files_prot = os.path.join(config["outputdir"], "04-communities", "prot", "{comm}_complete.csv"),
        community_files_nucl = os.path.join(config["outputdir"], "04-communities", "nucl", "{comm}_complete.csv")
    params:
        community = "{comm}",
        directory_prot = os.path.join(config["outputdir"], "04-communities", "prot"),
        directory_nucl = os.path.join(config["outputdir"], "04-communities", "nucl")
    run:
        community_spec = pd.read_csv(input.communities_file)
        community_spec = community_spec[community_spec.Community == int(params.community)]
        os.system("mkdir -p " + str(params.directory_prot)) 
        os.system("mkdir -p " + str(params.directory_nucl))
        for com_ind in range(len(community_spec.index)):
            protein_file = protein_files[np.where([(curr == str(list(community_spec.Organism)[com_ind])) \
                                                   for curr in transcriptomes])[0][0]]
            nucl_file = list(community_spec.Organism)[com_ind]
            protein_ind = nucl_file.split("/")[-1].split(".")[0]
            os.system("cp " + str(protein_file) + " " + os.path.join(params.directory_prot,
                                                                     str(protein_ind) + ".pep.fa"))
            os.system("cp " + str(nucl_file) + " " + os.path.join(params.directory_nucl,
                                                                  str(protein_ind) + ".fasta"))
            
        if len(community_spec.index) > 0:
            os.system("touch " + str(output.community_files_prot))
            os.system("touch " + str(output.community_files_nucl))