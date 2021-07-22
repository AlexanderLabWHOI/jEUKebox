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
import math
from Bio import SeqIO
import datetime
import ast
from snakemake.exceptions import print_exception, WorkflowError  

mydate = datetime.datetime.now()
datestring = mydate.strftime("%b%d")

# Generate community information using fastANI similarities and provided csv file.
rule create_communities:
    input:
        tsv_file = os.path.join(config["outputdir"], "02-fastani", "fastani_comparison.tsv")
    output:
        communities_file = os.path.join(config["outputdir"], "03-community_spec", "communities.csv")
    params:
        transcriptomes = transcriptomes,
        comm_df_loc = config["community_spec"]
    run:
        ani_file = pd.read_csv(input.tsv_file, sep = "\t",
                               header=None,names = ["file1","file2","ani","XX","XXX"])
        ani_file = ani_file[ani_file.file1 != ani_file.file2]
        unrelated = list(set(params.transcriptomes) - set(ani_file.file1) - set(ani_file.file2))
        ani_file = ani_file.sort_values(by = "ani")
        communities = pd.read_csv(params.comm_df_loc)
        output_file = pd.DataFrame(columns = ["Community", "Organism", "Proportion", "File_Inds", "Related?"])
        for community_ind in range(len(communities.index)):
            number_orgs = int(list(communities.NumberOrganisms)[community_ind])
            number_highsim = list(communities.NumberHighSimilarity)[community_ind]
            groups_highsim = list(communities.GroupsHighSimilarity)[community_ind]
            community_curr = list(communities.Communities)[community_ind]
            organisms = []
            proportions = []
            high_relation = []
           
            if number_highsim > 0:

                # the number of high-similarity groups cannot be <1 in this case
                if groups_highsim < 1:
                    groups_highsim = 1
                num_pergroup = math.floor(number_orgs / groups_highsim)
                while groups_highsim > 0:
                    top_related = []
                    ani_rows = 0
                    while (ani_rows < len(ani_file.index)) & (len(top_related) < num_pergroup):
                        candidate_1 = list(ani_file.file1)[ani_rows]
                        candidate_2 = list(ani_file.file2)[ani_rows]
                        if (len(top_related) == 0) | (candidate_1 in top_related) | (candidate_2 in top_related):
                            top_related.extend([candidate_1,candidate_2])
                            top_related = list(set(top_related))
                        ani_rows+=1
                   
                    groups_highsim-=1
                    organisms.extend(top_related)
                    
            high_relation.extend([1]*len(organisms))
                
            number_orgs-=len(organisms)
            high_relation.extend([0]*number_orgs)
            if number_orgs <= len(unrelated):
                organisms.extend(np.random.choice(unrelated, number_orgs, replace=False))
            else:
                organisms.extend(np.random.choice(unrelated, number_orgs, replace=True))
            
            prop_list = ast.literal_eval(communities.ProportionCommunity[community_ind])
            #prop_list = [n.strip() for n in prop_list]
            proportions.extend(prop_list)
            mmetsp_inds = [curr.split("/")[-1].split(".")[0] for curr in organisms]
            output_file = output_file.append(pd.DataFrame({"Community": community_curr,
                                             "Organism": organisms,
                                             "File_Inds": mmetsp_inds,
                                             "Proportion": proportions,
                                             "Related?": high_relation}))
            
        output_file.to_csv(output.communities_file)