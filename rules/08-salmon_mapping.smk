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


rule salmon_index:
    input:
        raw_reads_1 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R2.fastq.gz"),
        mock_assembly = os.path.join(config["outputdir"], "06-designer_assemblies",
                                     "designer_assembly_{comm}.fasta")
    output:
        os.path.join(config["outputdir"], "08-salmon_mapping",\
                     "{comm}_index", "versionInfo.json")   
    params:
        libtype = "A",
        indexname = os.path.join(config["outputdir"], "08-salmon_mapping",\
                     "{comm}_index"),
        kval = 31
    conda:
        os.path.join("..", "envs", "post-process-env.yaml")
    log:
        err = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_index.err"),
        out = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_index.out")
    shell:
        '''
        salmon index -t {input.mock_assembly} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        '''

rule salmon_quant:
    input:
        raw_reads_1 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R2.fastq.gz"),
        mock_assembly = os.path.join(config["outputdir"], "06-designer_assemblies",
                                     "designer_assembly_{comm}.fasta"),
        index_complete = os.path.join(config["outputdir"], "08-salmon_mapping",\
                                      "{comm}_index", "versionInfo.json")  
    output: 
        os.path.join(config["outputdir"], "08-salmon_mapping",\
                     "{comm}_quant", "quant.sf")  
    params:
        libtype = "A",
        indexname = os.path.join(config["outputdir"], "08-salmon_mapping",\
                     "{comm}_index"),
        outdir = os.path.join(config["outputdir"], "08-salmon_mapping",\
                     "{comm}_quant")
    conda:
        os.path.join("..", "envs", "post-process-env.yaml")
    log:
        err = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_quant.err"),
        out = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_quant.out")
    shell:
        '''
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.raw_reads_1} -2 {input.raw_reads_2} --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        '''
        
rule salmon_index_protein:
    input:
        raw_reads_1 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R2.fastq.gz"),
        mock_assembly = os.path.join(config["outputdir"], "06-designer_assemblies",
                                     "protein",
                                     "designer_assembly_{comm}.pep.fasta")
    output:
        os.path.join(config["outputdir"], "10-salmon_mapping-prot",\
                     "{comm}_index", "versionInfo.json")   
    params:
        libtype = "A",
        indexname = os.path.join(config["outputdir"], "10-salmon_mapping-prot",\
                     "{comm}_index"),
        kval = 21
    conda:
        os.path.join("..", "envs", "post-process-env.yaml")
    log:
        err = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_index_prot.err"),
        out = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_index_prot.out")
    shell:
        '''
        salmon index -t {input.mock_assembly} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        '''

rule salmon_quant_protein:
    input:
        raw_reads_1 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R1.fastq.gz"),
        raw_reads_2 = os.path.join(config["outputdir"], "07-simulated_raw_reads",
                                   "sim_raw_reads_{comm}_R2.fastq.gz"),
        mock_assembly = os.path.join(config["outputdir"], "06-designer_assemblies",
                                     "protein",
                                     "designer_assembly_{comm}.pep.fasta"),
        index_complete = os.path.join(config["outputdir"], "10-salmon_mapping-prot",\
                                      "{comm}_index", "versionInfo.json")  
    output: 
        os.path.join(config["outputdir"], "10-salmon_mapping-prot",\
                     "{comm}_quant", "quant.sf")  
    params:
        libtype = "A",
        indexname = os.path.join(config["outputdir"], "10-salmon_mapping-prot",\
                     "{comm}_index"),
        outdir = os.path.join(config["outputdir"], "10-salmon_mapping-prot",\
                     "{comm}_quant")
    conda:
        os.path.join("..", "envs", "post-process-env.yaml")
    log:
        err = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_quant_prot.err"),
        out = os.path.join(config["outputdir"], "logs", "salmon", "{comm}_quant_prot.out")
    shell:
        '''
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.raw_reads_1} -2 {input.raw_reads_2} --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        '''
