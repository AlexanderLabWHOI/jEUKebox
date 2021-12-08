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

mydate = datetime.datetime.now()
datestring = mydate.strftime("%b%d")
total_size_assembly = config["total_size"]

rule create_assemblies:
    input:
        community_files_prot = os.path.join(config["outputdir"], "04-communities",
                                            "prot", "{comm}", "{comm}_complete.csv"),
        community_files_nucl = os.path.join(config["outputdir"], "04-communities",
                                            "nucl", "{comm}", "{comm}_complete.csv"),
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
        mock_assembly_prot = os.path.join(config["outputdir"], "06-designer_assemblies", "protein",
                                     "designer_assembly_{comm}.pep.fasta"),
        assembly_spec = os.path.join(config["outputdir"], "06-designer_assemblies", "specs",
                                     "designer_assembly_{comm}_spec.csv"),
        pickled_dict = os.path.join(config["outputdir"], "04-communities", "community_id_dicts",
                                     "{comm}.pickle")
    params:
        directory_prot = os.path.join(config["outputdir"], "04-communities", "prot", "{comm}"),
        community_dir = os.path.join(config["outputdir"], "04-communities", "nucl", "{comm}"),
        outdir = config["outputdir"],
        total_size_metatranscriptome = total_size_assembly
    run:
        community_spec = pd.read_csv(input.communities_file)
        files_inassembly = [curr.split(".")[0] for curr in os.listdir(params.community_dir)]
        community_spec = community_spec[community_spec.File_Inds.isin(files_inassembly)] 
        
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
                if curr_string.split("/NCGR_PEP_ID=")[-1].split("_")[0] not in peptide_dict:
                    peptide_dict[curr_string.split("/NCGR_PEP_ID=")[-1].split("_")[0]] = record.id
                if record.id not in from_pep_dict:
                    from_pep_dict[record.id] = curr_string.split("/NCGR_PEP_ID=")[-1].split("_")[0]
        
        
        os.makedirs(os.path.dirname(output.pickled_dict), exist_ok=True)
        with open(output.pickled_dict,'wb') as outfile:
            pickle.dump(peptide_dict,outfile)
        logfile = open(os.path.join(params.outdir, "04-communities", "community_id_dicts", "logfile.txt"),'w')
        logfile.write("starting...\n")
                                            
        # we want to use the OrthoFinder output to ensure that each assembly contains contigs from
        # OG spanned by ALL species in assembly.
        
        gene_ct = pd.read_csv(input.orthfinder_genct, sep="\t")
        orthogroups = pd.read_csv(input.orthfinder_OGs, sep="\t")
        single_copy_OGs = pd.read_csv(input.orthfinder_SCGs, sep="\t",header=None,names=["SCs"])
        single_copy_OGs = list(single_copy_OGs.SCs)
        gene_ct = gene_ct.replace(0, np.nan)
        gene_ct = gene_ct.dropna(how='any', axis=0)
        to_write=[]
        to_write_prot=[]
        percent_dict = dict()
        
        # we want to include 10% of the single copy OGs at random.
        indices_OGs = np.random.choice(list(range(0,len(single_copy_OGs))), size=int(len(single_copy_OGs)*0.10), replace=False)
        selected_OGs = [single_copy_OGs[curr] for curr in indices_OGs]
        selected_peps = []
        [selected_peps.extend(list(orthogroups.loc[[curr in selected_OGs for curr in orthogroups.Orthogroup],:].\
                                   drop("Orthogroup",axis=1).reset_index(drop=True).loc[ind])) for ind \
                              in range(len(selected_OGs))]
        selected_peps = [tester for tester in selected_peps if str(tester) != 'nan']
        selected_nucls = [from_pep_dict[curr_pep] for curr_pep in selected_peps]
        
        # save the info about what organism each contig came from
        concordance = pd.DataFrame(columns=["Contig","Organism"])
        org_ids = []
        percentages = []
        
        leftover_size = int(params.total_size_metatranscriptome)# - len(selected_peps)
        all_record_list = list()
        all_record_list_prot = list()
        # for each member of the community we need to do something different
        for row_ind in range(len(community_spec.index)):
            org_id = list(community_spec.File_Inds)[row_ind]

            record_list = list(SeqIO.parse(list(community_spec.Organism)[row_ind], "fasta"))
            all_record_list.extend(record_list)
            record_list_prot = list(SeqIO.parse(list(community_spec.Protein_Files)[row_ind], "fasta"))
            all_record_list_prot.extend(record_list_prot)
            
            # the percentage of the assembly that we wish to be taken up by this organism 
            percentage = list(community_spec.Proportion)[row_ind]

            # is this organism among the related group?
            related_check = list(community_spec["Related?"])[row_ind]
            
            # orthogroups that contain this organism as well as others (regardless of how many)
            shared_ogs = orthogroups.loc[(~pd.isna(orthogroups[org_id + ".pep"])) & \
                                         (~pd.isna(orthogroups.drop(["Orthogroup",org_id + ".pep"],\
                                                                     axis=1))).any(axis=1),:]
            shared_og_campeps = []
            [shared_og_campeps.extend([curr2 for curr2 in curr.split(",") if curr2 not in selected_peps]) for curr in \
                     list(shared_ogs[org_id + ".pep"]) if str(curr) != "nan"]
            shared_og_campeps = list(set(shared_og_campeps))
            
            # orthogroups that contain this organism as well as others
            not_shared_ogs = orthogroups.loc[(~pd.isna(orthogroups[org_id + ".pep"])) & \
                                          pd.isna(orthogroups.drop(["Orthogroup",org_id + ".pep"],\
                                                                   axis=1)).all(axis=1),:]
            not_shared_og_campeps = []
            [not_shared_og_campeps.extend([curr2 for curr2 in curr.split(",") if curr2 not in shared_og_campeps]) for curr in \
                     list(not_shared_ogs[org_id + ".pep"]) if (str(curr) != "nan") & \
                     (str(curr) not in shared_og_campeps)]
            not_shared_og_campeps = list(set(not_shared_og_campeps))
                
            # the total number of contigs in the transcriptome assembly for this organism
            total_items_transcriptome = len(record_list) 
            
            # the desired size of the total metatranscriptome assembly times the percentage we want
            # this one to represent
            total_items = int(leftover_size * percentage)
            
            # this number can't exceed contigs available
            if total_items > total_items_transcriptome:
                total_items = total_items_transcriptome
            
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
            
            total_shared = len(shared_og_campeps)
            total_not_shared = len(not_shared_og_campeps)
            if total_shared <= shared_num:
                shared_num = total_shared - 1
            if total_not_shared <= not_shared_num:
                not_shared_num = total_not_shared - 1
            
            random_num_shared = np.random.choice(list(range(0,total_shared)), size=shared_num, replace=False)
            random_num_not_shared = np.random.choice(list(range(0,total_not_shared)), size=not_shared_num, replace=False)
            #random_num_not_shared = np.random.choice(list(set(list(range(0,total_not_shared))).\
            #                                              difference(set(random_num_shared))),\
            #                                      size=not_shared_num, replace=False)
            selected_nucls.extend([from_pep_dict[shared_og_campeps[random_num_curr].strip()] for \
                                  random_num_curr in random_num_shared])
            selected_nucls.extend([from_pep_dict[not_shared_og_campeps[random_num_curr].strip()] for \
                                  random_num_curr in random_num_not_shared])
            
            if org_id not in percent_dict:
                percent_dict[org_id] = percentage
            
        chosen_indices = list(np.random.choice(list(range(0,len(selected_nucls))), 
                                               params.total_size_metatranscriptome, replace=False))
        selected_nucls = list(set([selected_nucls[curr] for curr in \
                                   range(len(selected_nucls)) if curr in chosen_indices]))
        selected_peps = set([peptide_dict[curr] for curr in selected_nucls])
        org_ids = [curr.split("|")[-2] if len(curr.split("|")) > 1 else curr.split("|")[0] for curr in selected_nucls]
        percentages = [percent_dict[curr] if curr in percent_dict else np.random.random(1)[0] for curr in org_ids]
        #to_write.extend([curr for curr in record_list if \
        #                 any([nucl_sel in str(curr.id) for nucl_sel in selected_nucls])])
        check_1 = [curr for curr in all_record_list if (len((curr.id).split("|")) > 1) & 
                   ("|".join((curr.id).split("|")[2:4]) in selected_nucls)]
        to_write.extend(check_1)
        check_1_prot = [curr for curr in all_record_list_prot if curr.id in selected_peps]
        to_write_prot.extend(check_1_prot)
        
        #org_ids.extend([org_id]*len(to_write))
        #percentages.extend([percentage]*len(to_write))
            #to_write_prot.extend([curr for curr in record_list_prot if \
            #                      any([nucl_sel in str(curr.id) for nucl_sel in selected_peps])])\
        
        logfile.close()
        os.makedirs(os.path.dirname(output.mock_assembly), exist_ok=True)   
        os.makedirs(os.path.dirname(output.mock_assembly_prot), exist_ok=True) 
        to_write_distinct = []
        to_write_distinct_prot = []
        unique_ids = []
        unique_prots = []
        org_ids_distinct = []
        percentages_distinct = []
        for curr,curr_prot,percentage,org_id in zip(to_write,to_write_prot,percentages,org_ids):
            if (curr.id not in unique_ids) & (curr_prot.id not in unique_prots):
                unique_prots.append(curr_prot.id)
                unique_ids.append(curr.id)
                to_write_distinct.append(curr)
                to_write_distinct_prot.append(curr_prot)
                org_ids_distinct.append(org_id)
                percentages_distinct.append(percentage)
                     
        with open(output.mock_assembly, 'w') as handle:
            SeqIO.write(to_write_distinct, handle, 'fasta')   
        with open(output.mock_assembly_prot, 'w') as handle:
            SeqIO.write(to_write_distinct_prot, handle, 'fasta')
            
        
        concordance = concordance.append(pd.DataFrame({"Contig":unique_ids,
                                                       "Organism":org_ids_distinct,
                                                       "Proportion":percentages_distinct}),ignore_index=True)
            
        concordance.to_csv(output.assembly_spec)
