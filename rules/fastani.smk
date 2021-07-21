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

rule gen_transcriptome_names:
    input:
        transcriptomes = transcriptomes
    output:
        queries_file = os.path.join(config["outputdir"], "fastani_setup", "fastani_queries.txt"),
        references_file = os.path.join(config["outputdir"], "fastani_setup", "fastani_references.txt")
    params:
        transcriptomes = "\n".join(transcriptomes)
    shell:
        '''
        for transcriptome in {input.transcriptomes}; do
            echo -e $transcriptome >> {output.queries_file}
            echo -e $transcriptome >> {output.references_file}
        done
        '''
        
rule fastani:
    input:
        queries_file = os.path.join(config["outputdir"], "fastani_setup", "fastani_queries.txt"),
        references_file = os.path.join(config["outputdir"], "fastani_setup", "fastani_references.txt")
    output:
        tsv_file = os.path.join(config["outputdir"], "fastani", "fastani_comparison.tsv")
    params:
        transcriptomes = " ".join(transcriptomes)
    conda:
        os.path.join("..", "envs", "workflow-env.yaml")
    shell:
        '''
        fastANI --ql {input.queries_file} --rl {input.references_file} -o {output.tsv_file}
        '''