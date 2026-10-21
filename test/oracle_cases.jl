# Cases compared with the upstream oracle (oracle/rust_upstream). Shared by
# test/rust_oracle.jl and tools/make_reference.jl.

const DATA_DIRECTORY = joinpath(@__DIR__, "data")
const REFERENCE_DIRECTORY = joinpath(@__DIR__, "reference")

const PRIMITIVE_INPUTS = ["primitives.fastq", "primitives.fasta"]
const PRIMITIVE_PARAMETERS = [
    (k=5, c=1, w=3, l=2, minimum_bq=25),
    (k=17, c=11, w=10, l=3, minimum_bq=25),
    (k=27, c=2, w=28, l=4, minimum_bq=30),
]

const AMPLICON_INPUTS = ["amplicons_sample1.fastq.gz", "amplicons_sample2.fastq.gz"]
# (name, upstream `savont asv` arguments, equivalent SavontConfig keywords)
const MARKER_RUNS = [
    ("default_blockmers", ["--use-blockmers"], (; use_blockmers=true)),
    ("single_strand", ["--single-strand", "--use-blockmers"], (; single_strand=true, use_blockmers=true)),
    ("k13_bq20", ["--kmer-size", "13", "--minimum-base-quality", "20"], (; kmer_size=13, minimum_base_quality=20)),
]

primitive_reference_name(input, p) = "primitives_$(splitext(input)[2][2:end])_k$(p.k)_c$(p.c)_w$(p.w)_l$(p.l)_q$(p.minimum_bq).tsv.gz"
marker_reference_name(name) = "markers_$(name).tsv.gz"
const STATS_REFERENCE = "stats.tsv.gz"

primitive_command(oracle, input, p) =
    `$oracle primitives $(joinpath(DATA_DIRECTORY, input)) $(p.k) $(p.c) $(p.w) $(p.l) $(p.minimum_bq)`
marker_command(oracle, arguments) =
    `$oracle markers $(joinpath.(DATA_DIRECTORY, AMPLICON_INPUTS)) --threads 2 $arguments`
stats_command(oracle) = `$oracle stats`

# Twin reads, primary and polymorphic-marker clustering: (name, inputs, upstream
# arguments, SavontConfig keywords). Blockmers need k + l in upstream's list of
# supported sizes.
const NOISY_1 = ["noisy_1pct.fastq.gz"]
const NOISY_2 = ["noisy_2pct.fastq.gz"]
const CLUSTER_RUNS = [
    ("default", AMPLICON_INPUTS, String[], (;)),
    ("blockmers", AMPLICON_INPUTS, ["--use-blockmers"], (; use_blockmers=true)),
    ("k13", AMPLICON_INPUTS, ["--kmer-size", "13"], (; kmer_size=13)),
    ("strict_threshold", AMPLICON_INPUTS, ["--primary-clustering-threshold", "0.985", "--min-cluster-size", "3"],
     (; primary_clustering_threshold=0.985, min_cluster_size=3)),
    ("single_strand_c7", AMPLICON_INPUTS, ["--single-strand", "-c", "7"], (; single_strand=true, sampling_rate=7)),
    ("low_polymorphism", AMPLICON_INPUTS, ["--low-polymorphism"], (; low_polymorphism=true)),
    ("noisy1", NOISY_1, ["--min-read-length", "1000"], (; min_read_length=1000)),
    ("noisy1_blockmers", NOISY_1, ["--min-read-length", "1000", "--use-blockmers"], (; min_read_length=1000, use_blockmers=true)),
    ("noisy2", NOISY_2, ["--min-read-length", "1000"], (; min_read_length=1000)),
    ("noisy2_blockmers", NOISY_2, ["--min-read-length", "1000", "--use-blockmers"], (; min_read_length=1000, use_blockmers=true)),
    ("noisy2_one_iteration", NOISY_2, ["--min-read-length", "1000", "--max-iterations-recluster", "1"],
     (; min_read_length=1000, max_iterations_recluster=1)),
]
cluster_reference_name(name) = "cluster_$(name).tsv.gz"
cluster_command(oracle, inputs, arguments) =
    `$oracle cluster $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`
