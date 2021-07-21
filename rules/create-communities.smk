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

# Generate community information using fastANI similarities and provided csv file.
rule create_communities:
    input:
        tsv_file = os.path.join(config["outputdir"], "fastani", "fastani_comparison.tsv")
    output:
        communities_file = os.path.join(config["outputdir"], "communities.csv")
    params:
        transcriptomes = transcriptomes,
        comm_df_loc = config["community_spec"]
    run:
        ani_file = pd.read_csv(input.tsv_file, sep = "\t",
                               header=None,names = ["file1","file2","ani","XX","XXX"])
        ani_file = ani_file[ani_file.file1 != ani_file.file2]
        unrelated = list(set(params.transcriptomes) - set(ani_file.file1) - set(ani_file.file2))
        sorted_by_ani = ani_file.sort_values(by = "ani")
        communities = pd.read_csv(params.comm_df_loc)
        output_file = pd.DataFrame(columns = ["Community", "Organism", "Proportion"])
        for community_ind in range(len(communities.index)):
            number_orgs = int(list(communities.NumberOrganisms)[community_ind])
            number_highsim = list(communities.NumberHighSimilarity)[community_ind]
            community_curr = list(communities.Communities)[community_ind]
            organisms = []
            proportions = []
            if community_curr != 6:
                number_orgs = number_orgs - number_highsim
            else:
                unrelated = list(set(params.transcriptomes))
            if number_highsim == 2:
                organisms.extend(list(ani_file.loc[:,["file1","file2"]].iloc[0]))
            elif (number_highsim == 4) & (community_curr == 3):
                # this is the community where we want all 4 to be related to _each other_
                top_related = []
                ani_rows = 0
                while (ani_rows < len(ani_file.index)) & (len(top_related) < 4):
                    candidate_1 = list(ani_file.file1)[ani_rows]
                    candidate_2 = list(ani_file.file2)[ani_rows]
                    if (len(top_related) == 0) | (candidate_1 in top_related) | (candidate_2 in top_related):
                        top_related.extend([candidate_1,candidate_2])
                        top_related = list(set(top_related))
                    ani_rows+=1
                        
                 
                organisms.extend(top_related)
            elif (number_highsim == 4) & (community_curr == 5):
                # this is the community where we want all 4 to be related, but not to each other
                top_related = []
                ani_rows = 0
                while (ani_rows < len(ani_file.index)) & (len(top_related) < 4):
                    candidate_1 = list(ani_file.file1)[ani_rows]
                    candidate_2 = list(ani_file.file2)[ani_rows]
                    if (len(top_related) == 0) | ((candidate_1 not in top_related) & (candidate_2 not in top_related)):
                        top_related.extend([candidate_1,candidate_2])
                        top_related = list(set(top_related))
                        
                    ani_rows+=1
                        
                organisms.extend(top_related)
                
            if number_orgs <= len(unrelated):
                organisms.extend(np.random.choice(unrelated, number_orgs, replace=False))
            else:
                organisms.extend(np.random.choice(unrelated, number_orgs, replace=True))
            
            proportions.extend(list(communities.ProportionCommunity[community_ind]))
            mmetsp_inds = [curr.split("/")[-1] for curr in organisms]
            output_file = output_file.append(pd.DataFrame({"Community": community_curr,
                                             "Organism": organisms,
                                             "MMETSP_inds": mmetsp_inds,
                                             "Proportion": proportions}))
            
        output_file.to_csv(output.communities_file)