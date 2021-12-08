import os
import pandas as pd
import numpy as np

mmetsp_dir = "/vortexfs1/omics/alexander/data/mmetsp/"
mmetsp_files = os.listdir(mmetsp_dir)

mmetsp_prot_files = sorted([curr for curr in mmetsp_files if "clean.pep" in curr])
mmetsp_nucl_files = [curr.split("_")[0] + "_clean.fasta" for curr in mmetsp_prot_files] 
#sorted([curr for curr in mmetsp_files if "clean.fasta" in curr])

rand_mmetsp_inds = np.random.choice(list(range(len(mmetsp_prot_files))),12,replace=False)
rand_mmetsp_inds = list(range(len(mmetsp_prot_files)))
np.random.shuffle(rand_mmetsp_inds)
rand_mmetsp_inds = list(rand_mmetsp_inds)[0:12] 

rand_mmetsps_prot = [os.path.join(mmetsp_dir,mmetsp_prot_files[curr]) for curr in list(range(len(mmetsp_prot_files))) if curr in rand_mmetsp_inds]
rand_mmetsps_nucl = [os.path.join(mmetsp_dir,mmetsp_nucl_files[curr]) for curr in list(range(len(mmetsp_nucl_files))) if curr in rand_mmetsp_inds]

with open("mmetsp_ids.txt","w") as f:
    f.writelines("proteinfiles:")
    f.writelines(["- "+curr+"\n" for curr in rand_mmetsps_prot])
    f.writelines("nuclfiles:")
    for curr in rand_mmetsps_nucl:
        f.writelines("- "+curr+"\n")
print(rand_mmetsps_nucl)
print(rand_mmetsps_prot)
