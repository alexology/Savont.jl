# `savont download` (upstream `databases.rs` and `download.rs`): reference
# databases fetched into `<location>/<keyword>`.

const DATABASE_SOURCES = Dict(
    "emu-1" => (description="EMU default 16S rRNA database",
                urls=["https://osf.io/8qcwd/download"]),
    "silva-138.2" => (description="SILVA SSU Ref NR99 v138.2",
                      urls=["https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz",
                            "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz"]),
    "greengenes2-2024.09" => (description="GreenGenes2 2024.09 species-level trainset from DADA2",
                              urls=["https://zenodo.org/records/14169078/files/gg2_2024_09_toSpecies_trainset.fa.gz"]),
    "unite-10.0" => (description="UNITE ITS 10.0, dynamic divergence thresholding, all eukaryotes",
                     urls=["https://s3.hpc.ut.ee/plutof-public/original/e861a3d6-54f4-42dc-882a-5f129beac39a.tgz"]),
)

# Extract a gzip-compressed tar archive into `destination`, moving the archive's
# top-level entries there (as `tar -xzf archive -C destination`).
function _extract_tgz(archive, destination)
    staging = mktempdir(destination)
    try
        open(archive) do io
            Tar.extract(GzipDecompressorStream(io), staging)
        end
        for name in readdir(staging)
            mv(joinpath(staging, name), joinpath(destination, name); force=true)
        end
    finally
        rm(staging; force=true, recursive=true)
    end
end

function _gunzip(path)
    target = chopsuffix(path, ".gz")
    open(path) do input
        open(target, "w") do output
            write(output, GzipDecompressorStream(input))
        end
    end
    rm(path)
    return target
end

"""
    preprocess_unite(source, destination)

Rewrite a UNITE FASTA (`Species|Accession|SH|type|k__…;s__…` headers) as
`unite_sequences.fasta`, named by the SH accession only (minimap2 rejects names
longer than 255 characters), and `unite_taxonomy.tsv` (`sh_id`, lineage). Other
headers are kept, cut to 255 characters.
"""
function preprocess_unite(source::AbstractString, destination::AbstractString)
    open(joinpath(destination, "unite_sequences.fasta"), "w") do fasta
        open(joinpath(destination, "unite_taxonomy.tsv"), "w") do taxonomy
            println(taxonomy, "sh_id\ttaxonomy")
            for line in eachline(source)
                if startswith(line, '>')
                    parts = split(line[2:end], '|'; limit=5)
                    if length(parts) == 5
                        println(fasta, '>', parts[3])
                        println(taxonomy, parts[3], '\t', parts[5])
                    else
                        println(fasta, '>', first(line[2:end], 255))
                    end
                else
                    println(fasta, line)
                end
            end
        end
    end
    return destination
end

function _download_one(keyword, destination, urls)
    if keyword == "emu-1"
        archive = joinpath(destination, "emu_default.tar.gz")
        Downloads.download(urls[1], archive)
        _extract_tgz(archive, destination); rm(archive)
        subdirectory = joinpath(destination, "emu_default")
        for name in readdir(subdirectory)
            mv(joinpath(subdirectory, name), joinpath(destination, name); force=true)
        end
        rm(subdirectory; recursive=true)
    elseif keyword == "silva-138.2"
        for url in urls
            Downloads.download(url, joinpath(destination, last(split(url, '/'))))
        end
        _gunzip(joinpath(destination, last(split(urls[2], '/'))))
    elseif keyword == "greengenes2-2024.09"
        Downloads.download(urls[1], joinpath(destination, last(split(urls[1], '/'))))
    elseif keyword == "unite-10.0"
        archive = joinpath(destination, "unite-10.0.tgz")
        Downloads.download(urls[1], archive)
        _extract_tgz(archive, destination); rm(archive)
        names = filter(n -> endswith(n, ".fasta") && n != "unite_sequences.fasta" && isfile(joinpath(destination, n)), readdir(destination))
        isempty(names) && error("No .fasta file found after UNITE extraction")
        original = joinpath(destination, first(names))
        preprocess_unite(original, destination); rm(original)
    else
        throw(_unknown_keyword(keyword))
    end
end

"""
    download_database(location, keywords...; urls=Dict()) -> Vector{String}

Download the databases `keywords` (`emu-1`, `silva-138.2`,
`greengenes2-2024.09`, `unite-10.0`) as upstream `savont download`: each into
`<location>/<keyword>`, unpacked and prepared for [`load_database`](@ref), with
a `.savont_db` marker. `urls` overrides the download addresses of a keyword
(a vector in upstream's order). Returns the database directories.
"""
function download_database(location::AbstractString, keywords::AbstractString...; urls::AbstractDict=Dict{String,Vector{String}}())
    for keyword in keywords
        haskey(DATABASE_SOURCES, keyword) || throw(_unknown_keyword(keyword))
    end
    directories = String[]
    for keyword in keywords
        destination = mkpath(joinpath(location, keyword))
        @info "Downloading '$keyword' ($(DATABASE_SOURCES[keyword].description)) to $destination"
        _download_one(keyword, destination, get(urls, keyword, DATABASE_SOURCES[keyword].urls))
        write(joinpath(destination, DATABASE_MARKER), keyword)
        push!(directories, destination)
    end
    return directories
end
