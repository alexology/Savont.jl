# Command-line savont:
#
#     julia -t 4 --project=path/to/Savont.jl path/to/Savont.jl/bin/savont.jl asv reads.fastq.gz -o savont-out
#
# Subcommands and options are those of savont 0.8.0 (`--help` lists them):
# asv, classify, sintax, export, download.

using Savont

exit(Savont.main(ARGS))
