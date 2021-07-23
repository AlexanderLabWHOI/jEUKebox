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
import pickle
from Bio import SeqIO
import datetime
from snakemake.exceptions import print_exception, WorkflowError  

rule rSubRead:
    input:
        mock_assembly = os.path.join(config["outputdir"], "06-designer_assemblies",
                                     "designer_assembly_{comm}.fasta"),
        assembly_spec = os.path.join(config["outputdir"], "06-designer_assemblies", "specs",
                                     "designer_assembly_{comm}_spec.csv"),
    output:
        raw_reads_1 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R2.fastq.gz"),
        raw_spec = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "spec", "raw_read_spec_{comm}.csv")
    params:
        Rscript = os.path.join("..", "scripts", "RsubReadsRun.R")
    conda:
        os.path.join("..", "envs", "R-env.yaml")
    shell:
        '''
        Rscript --vanilla {params.Rscript} {input.mock_assembly} {input.assembly_spec} {output.raw_reads_1} {output.raw_spec}
        '''