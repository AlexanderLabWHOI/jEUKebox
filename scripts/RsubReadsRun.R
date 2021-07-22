install.packages("pacman", repos='http://cran.us.r-project.org')
pacman::p_load(BiocManager,Biostrings,Rsubread)

args = commandArgs(trailingOnly=TRUE)

assembly_file = args[1]
communities_file = args[2]
assembly_spec = args[3]
output_raw_reads = args[4]