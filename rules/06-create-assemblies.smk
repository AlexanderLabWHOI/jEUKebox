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
        community_files_prot = os.path.join(config["outputdir"], "04-communities",
                                            "prot", "{comm}_complete.csv"),
        community_files_nucl = os.path.join(config["outputdir"], "04-communities",
                                            "nucl", "{comm}_complete.csv")
        communities_file = os.path.join(config["outputdir"], "03-community_spec",
                                        "communities.csv"),
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
    output:
        mock_assembly = os.path.join(config["outputdir"], "06-designer_assemblies",
                                     "designer_assembly_{comm}.fasta"),
        pickled_dict = os.path.join(config["outputdir"], "04-communities", "community_id_dicts",
                                     "{comm}.pickle")
    params:
        directory_prot = os.path.join(config["outputdir"], "04-communities", "prot"),
        community_dir = os.path.join(config["outputdir"], "04-communities", "nucl")
    run:
        community_spec = pd.read_csv(input.communities_file)
        files_inassembly = os.listdir(params.community_dir)
        community_spec = community_spec[community_spec.MMETSP_inds.isin(files_inassembly)] 
        
        # We need a way to translate between the peptide IDs in the OrthoFinder output 
        # and the nucleotide IDs we wish to use.
        peptide_files = os.listdir(params.directory_prot)
        nucleotide_files = os.listdir(params.community_dir)
        peptide_ids = []
        nucle_ids_fromprot = []
        peptide_dict = dict()
        from_pep_dict = dict()
        
        for peptide_file in peptide_files:
            for record in SeqIO.parse(os.path.join(params.directory_prot, peptide_file), "fasta"):
                curr_string = record.description
                peptide_ids.append(record.id)
                nucle_ids_fromprot.append(curr_string.split("/NCGR_PEP_ID=")[-1].split("_")[0])
                peptide_dict[curr_string.split("/NCGR_PEP_ID=")[-1].split("_")[0]] = record.id
                from_pep_dict[record.id] = curr_string.split("/NCGR_PEP_ID=")[-1].split("_")[0]
        
        with open(output.pickled_dict,'wb') as outfile:
            pickle.dump(peptide_dict,outfile)
                                            
        # we want to use the OrthoFinder output to ensure that each assembly contains contigs from
        # OG spanned by ALL species in assembly.
        
        gene_ct = pd.read_csv(input.orthfinder_genct, sep="\t")
        orthogroups = pd.read_csv(input.orthfinder_OGs, sep="\t")
        single_copy_OGs = pd.read_csv(input.orthfinder_SCGs, sep="\t",header=None,names=["SCs"])
        single_copy_OGs = list(single_copy_OGs.SCs)
        gene_ct = gene_ct.replace(0, np.nan)
        gene_ct = gene_ct.dropna(how='any', axis=0)
        to_write=[]
        
        # we want to include 10% of the single copy OGs at random.
        indices_OGs = random.choice(list(range(0,len(single_copy_OGs))), size=int(len(single_copy_OGs)*0.10), replace=False)
        selected_OGs = [single_copy_OGs[curr] for curr in indices_OGs]
        selected_peps = []
        [selected_peps.extend(list(orthogroups.loc[[curr in selected_OGs for curr in orthogroups.Orthogroup],:].\
                                   drop("Orthogroup",axis=1).loc[ind])) for ind in len(selected_OGs)]
        selected_peps = [tester for tester in selected_peps if str(tester) != 'nan']
        
        # for each member of the community we need to do something different
        for row_ind in range(len(community_spec.index)):
            org_id = list(community_spec.File_Inds)[row_ind]

            record_list = list(SeqIO.parse(list(community_spec.Organism)[row_ind], "fasta"))
            
            # the percentage of the assembly that we wish to be taken up by this organism 
            percentage = list(community_spec.Proportion)[row_ind]

            # is this organism among the related group?
            related_check = list(community_spec["Related?"])[row_ind]
            
            # orthogroups that contain this organism as well as others (regardless of how many)
            shared_ogs = orthogroups.loc[(~pd.isna(orthogroups[org_id])) & \
                                         (~pd.isna(orthogroups.drop(["Orthogroup",org_id]))).any(axis=1),:]
            
            # orthogroups that contain this organism as well as others
            not_shared_ogs = orthogroups.loc[(~pd.isna(orthogroups[org_id])) & \
                                          pd.isna(orthogroups.drop(["Orthogroup",org_id])).all(axis=1),:]
                
            # the total number of contigs in the transcriptome assembly for this organism
            total_items = len(record_list) 
            
            if related_check == 1:
                # the number of contigs we wish to include from shared OGs
                shared_num = int(total_items * 0.75)
                # the number of contigs we wish to include from species-specific OGs
                not_shared_num = int(total_items * 0.25)
            else:
                # the number of contigs we wish to include from shared OGs
                shared_num = int(total_items * 0.75)
                # the number of contigs we wish to include from species-specific OGs
                not_shared_num = int(total_items * 0.25)
            
            random_num = list(np.random.randint(0, high=total_items, size=int(total_items * percentage)))
            random_num = random.choice(list(range(0,total_items)), size=int(total_items * percentage), replace=False)
            to_write.extend([record_list[random_num_curr] for random_num_curr in random_num])
            
        with open(output.mock_assembly, 'w') as handle:
            SeqIO.write(to_write, handle, 'fasta')