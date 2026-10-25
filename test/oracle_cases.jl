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
const MIXED = ["mixed.fastq.gz"]
const CHIMERAS = ["chimeras.fastq.gz"]
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
    ("mixed", MIXED, ["--low-polymorphism", "--min-read-length", "80"], (; low_polymorphism=true, min_read_length=80)),
    ("mixed_hpc_masked", MIXED, ["--low-polymorphism", "--min-read-length", "80", "--use-hpc", "--mask-low-quality"],
     (; low_polymorphism=true, min_read_length=80, use_hpc=true, mask_low_quality=true)),
    ("chimeras", CHIMERAS, String[], (;)),
]
cluster_reference_name(name) = "cluster_$(name).tsv.gz"
cluster_command(oracle, inputs, arguments) =
    `$oracle cluster $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`
# Initial POA consensus of each final cluster (stage 4a), on the same runs.
consensus_reference_name(name) = "consensus_$(name).tsv.gz"
consensus_command(oracle, inputs, arguments) =
    `$oracle consensus $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`
# Pileups and polishing (stage 4b), on the same runs.
polish_reference_name(name) = "polish_$(name).tsv.gz"
polish_command(oracle, inputs, arguments) =
    `$oracle polish $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`
# Merging of similar consensuses (stage 5), on the same runs.
merge_reference_name(name) = "merge_$(name).tsv.gz"
merge_command(oracle, inputs, arguments) =
    `$oracle merge $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`
# Chimera detection (stage 6), on the same runs.
chimera_reference_name(name) = "chimera_$(name).tsv.gz"
chimera_command(oracle, inputs, arguments) =
    `$oracle chimera $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`

# The whole pipeline with its output files (stages 5-7 and outputs as main.rs),
# on the clustering runs and on runs of the options handled only there.
const ASV_RUNS = vcat([(name, inputs, arguments, keywords) for (name, inputs, arguments, keywords) in CLUSTER_RUNS], [
    ("pooled", AMPLICON_INPUTS, ["--pooled-samples"], (; pooled_samples=true)),
    ("pooled_low_polymorphism", AMPLICON_INPUTS, ["--pooled-samples", "--low-polymorphism"],
     (; pooled_samples=true, low_polymorphism=true)),
    ("skip_chimera_detection", AMPLICON_INPUTS, ["--skip-chimera-detection"], (; skip_chimera_detection=true)),
])
asv_reference_name(name) = "asv_$(name).tsv.gz"
asv_command(oracle, inputs, arguments) =
    `$oracle asv $(joinpath.(DATA_DIRECTORY, inputs)) --threads 2 $arguments`

# Partial order alignment with savont's settings: (input, "band" or "noband", isa).
# `isa = :sse41` sets SPOARS_FORCE_ISA so that spoars uses 128-bit vectors.
const POA_RUNS = [(input, mode, isa) for input in ["poa_plain.fastq.gz", "poa_shifted.fastq.gz", "poa_shifted_short.fastq.gz"]
                  for (mode, isa) in (("noband", :avx2), ("band", :avx2), ("band", :sse41))]
poa_reference_name(input, mode, isa) = "$(first(split(input, '.')))_$(mode)_$(isa).tsv.gz"
poa_command(oracle, input, mode, isa) =
    addenv(`$oracle poa $(joinpath(DATA_DIRECTORY, input)) $mode`, "SPOARS_FORCE_ISA" => (isa === :sse41 ? "sse41" : ""))

# Taxonomy subcommands on test/data/taxonomy (`taxonomy` command of the oracle,
# which prints every output file). The database is copied to a temporary
# directory because upstream caches its minimap2 index next to the FASTA; one
# thread makes upstream's ties between references deterministic.
const TAXONOMY_DIRECTORY = joinpath(DATA_DIRECTORY, "taxonomy")
# (name, ASV directory, database, upstream arguments, keywords)
const CLASSIFY_RUNS = [
    ("emu_single", "asv_single", "emu-1", String[], (;)),
    ("emu_pooled", "asv_pooled", "emu-1", String[], (;)),
    ("emu_headers", "asv_headers", "emu-1", String[], (;)),
    ("emu_thresholds", "asv_single", "emu-1", ["--species-threshold", "99.8", "--genus-threshold", "96", "--detailed-unclassified"],
     (; species_threshold=99.8, genus_threshold=96.0, detailed_unclassified=true)),
    ("silva_single", "asv_single", "silva-138.2", String[], (;)),
    ("gg2_single", "asv_single", "greengenes2-2024.09", String[], (;)),
    ("unite_pooled", "asv_pooled", "unite-10.0", String[], (;)),
]
const SINTAX_RUNS = [
    ("emu_single", "asv_single", "emu-1", String[], (;)),
    ("emu_options", "asv_single", "emu-1", ["--n-iter", "50", "--min-bootstrap", "0.6", "--detailed-unclassified"],
     (; n_iter=50, min_bootstrap=0.6, detailed_unclassified=true)),
    ("silva_single", "asv_single", "silva-138.2", String[], (;)),
    ("gg2_pooled", "asv_pooled", "greengenes2-2024.09", String[], (;)),
    ("unite_headers", "asv_headers", "unite-10.0", String[], (;)),
]
# (name, run directories, upstream arguments, keywords)
const EXPORT_RUNS = [
    ("default", ["export_a", "export_b", "export_c"], String[], (;)),
    ("no_fuzzy", ["export_a", "export_b", "export_c"], ["--no-fuzzy"], (; fuzzy=false)),
    ("relabel", ["export_a", "export_b"], ["--relabel", "x", "y", "z"], (; relabel=["x", "y", "z"])),
    ("single", ["asv_single"], String[], (;)),
]
classify_reference_name(name) = "classify_$(name).tsv.gz"
sintax_reference_name(name) = "sintax_$(name).tsv.gz"
export_reference_name(name) = "export_$(name).tsv.gz"
function taxonomy_db_command(oracle, subcommand, input, db, arguments)
    directory = mktempdir()
    cp(joinpath(TAXONOMY_DIRECTORY, db), joinpath(directory, db))
    return `$oracle taxonomy $subcommand -i $(joinpath(TAXONOMY_DIRECTORY, input)) -o $(joinpath(directory, "out"))
            -d $(joinpath(directory, db)) -t 1 $arguments`
end
classify_command(oracle, input, db, arguments) = taxonomy_db_command(oracle, "classify", input, db, arguments)
sintax_command(oracle, input, db, arguments) = taxonomy_db_command(oracle, "sintax", input, db, arguments)
export_command(oracle, inputs, arguments) =
    `$oracle taxonomy export -i $(joinpath.(TAXONOMY_DIRECTORY, inputs)) -o $(joinpath(mktempdir(), "out")) $arguments`
