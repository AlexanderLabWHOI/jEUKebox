install.packages("pacman", repos='http://cran.us.r-project.org')
pacman::p_load(BiocManager,Biostrings,Rsubread)

args = commandArgs(trailingOnly=TRUE)

assembly_file = args[1]
#communities_file = args[2]
assembly_spec = args[2]
output_raw_reads = args[3]
output_raw_spec = args[4]

raw_read_stub = strsplit(output_raw_reads,"\\.")[1]
output_name = file.path(rawreads_dir,assembly_stub)
assembly_spec = read.csv(assembly_spec)

dir.create(output_name, showWarnings = FALSE, recursive=TRUE)
scanned_info = Rsubread::scanFasta(

        # the file containing the transcript database
        fasta_file,

        # manipulating transcript names
        simplify.transcript.names = FALSE,

        # miscellaneous options
        quiet = FALSE)

nsequences <- nrow(scanned_info) - sum(scanned_info$Duplicate)

# Assign a random TPM value to each non-duplicated transcript sequence
expressionlevels <- rep(0, nrow(scanned_info))
scaling_for_expr <- assembly_spec$Proportion
expressionlevels[!scanned_info$Duplicate] <- rexp(nsequences) # exponential distribution


scaling_for_expr <- c(scaling_for_expr, rep(0.5,length(expressionlevels)-length(scaling_for_expr)))
expressionlevels <- expressionlevels * scaling_for_expr

assembly_spec$ExpressionLevel <- expressionlevels
write.csv(assembly_spec, output_raw_spec)

transcript_info = Rsubread::simReads(

    # the transcript database and the wanted abundances
    assembly_file,
    #expression.levels,
    expressionlevels,

    # the name of the output
    raw_read_stub,

    # options on the output
    library.size = 1000000,
    read.length = 75,
    # logical indicating if the true mapping location of reads or
    # read-pairs should be encoded into the read names
    truth.in.read.names = TRUE,

    # simulating sequencing errors
    simulate.sequencing.error = TRUE,
    quality.reference = NULL,

    # parameters for generating paired-end reads.
    paired.end = TRUE,
    fragment.length.min = 100,
    fragment.length.max = 500,
    fragment.length.mean = 180,
    fragment.length.sd = 40,

    # manipulating transcript names
    simplify.transcript.names = FALSE)