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

rule orthofinder:
    input:
        community_files_prot = os.path.join(config["outputdir"], "04-communities",
                                            "prot", "{comm}_complete.csv"),
        community_files_nucl = os.path.join(config["outputdir"], "04-communities",
                                            "nucl", "{comm}_complete.csv")
    output:
        orthfinder_genct = os.path.join(config["outputdir"], "05-orthofinder", 
                                        "orthofinder_{comm}",
                                        "Results_" + str(datestring),"Orthogroups",
                                        "Orthogroups.GeneCount.tsv"),
        orthfinder_OGs = os.path.join(config["outputdir"], "05-orthofinder",
                                      "orthofinder_{comm}",
                                      "Results_" + str(datestring),
                                      "Orthogroups","Orthogroups.tsv"),
        orthfinder_SCGs = os.path.join(config["outputdir"], "05-orthofinder",
                                      "orthofinder_{comm}",
                                      "Results_" + str(datestring),
                                      "Orthogroups", "Orthogroups_SingleCopyOrthologues.txt")
      
    params:
        outfold = os.path.join(config["outputdir"], "orthofinder", "orthofinder_{comm}"),
        outfolder = config["outputdir"],
        directory_prot = os.path.join(config["outputdir"], "04-communities", "prot")
    conda:
        os.path.join("envs","workflow-env.yaml")
    shell:
        '''
        rm -rf {params.outfold}
        orthofinder -t 108 -o {params.outfold} -f {params.directory_prot}
        '''