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

Both files are written to temporary files and published together only when the
whole source has been read, so a missing or unreadable
source, a source that is one of the outputs, an output path that is a directory,
or a failure while publishing never leaves the existing files truncated, deleted
or mismatched (deviation D19).
"""
function preprocess_unite(source::AbstractString, destination::AbstractString)
    isfile(source) || throw(ArgumentError("UNITE source FASTA $(repr(source)) does not exist"))
    outputs = joinpath.(destination, ["unite_sequences.fasta", "unite_taxonomy.tsv"])
    _check_publishable(outputs)
    fasta_path, fasta = mktemp(destination; cleanup=false)
    taxonomy_path, taxonomy = mktemp(destination; cleanup=false)
    try
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
        close(fasta); close(taxonomy)
        _publish_files([fasta_path, taxonomy_path], outputs)
    finally
        close(fasta); close(taxonomy)
        rm(fasta_path; force=true); rm(taxonomy_path; force=true)
    end
    return destination
end

# Output paths must be absent or regular files: anything else (a directory, for
# example) is refused rather than replaced.
function _check_publishable(targets)
    for target in targets
        ispath(target) && !isfile(target) &&
            throw(ArgumentError("cannot write $(repr(target)): the path exists and is not a regular file"))
    end
end

"""
    _publish_files(sources, targets)

Move each complete file `sources[i]` to `targets[i]` as one step for the set: the
existing targets are first moved aside, and if any move fails the new files are
removed and the previous ones restored before the error is rethrown (if even the
restore fails, the error names the backups that were kept). Targets must be absent
or regular files. This keeps the outputs consistent with each other; it does not
make the set atomic for readers running at the same time.
"""
function _publish_files(sources, targets)
    # Plain renames only: `Base.mv` falls back to copying and deleting when a
    # rename fails, which can empty a target before failing (review SV16).
    _check_publishable(targets)
    suffix = ".previous-" * string(rand(UInt32); base=16)
    backups = Pair{String,String}[]          # target => backup
    placed = String[]
    try
        for target in targets
            isfile(target) || continue
            _rename(target, target * suffix)
            push!(backups, target => target * suffix)
        end
        for (source, target) in zip(sources, targets)
            _rename(source, target)
            push!(placed, target)
        end
    catch
        kept = String[]
        for target in placed
            try rm(target; force=true) catch end
        end
        for (target, backup) in backups
            try
                ispath(target) && rm(target; force=true)
                _rename(backup, target)
            catch
                push!(kept, backup)
            end
        end
        isempty(kept) || @error "Could not restore previous files; backups kept: $(join(kept, ", "))"
        rethrow()
    end
    for (_, backup) in backups
        try rm(backup; force=true) catch end
    end
    return targets
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

Each database is downloaded and prepared in a staging directory next to its
destination, which replaces the previous directory only once complete and
loadable (a minimal check of its taxonomy and FASTA): a failed download, or one that delivers
something other than the expected files (an HTML error page, for example), leaves
the previous database as it was (deviation D20).
"""
function download_database(location::AbstractString, keywords::AbstractString...; urls::AbstractDict=Dict{String,Vector{String}}())
    for keyword in keywords
        haskey(DATABASE_SOURCES, keyword) || throw(_unknown_keyword(keyword))
    end
    mkpath(location)
    directories = String[]
    for keyword in keywords
        destination = joinpath(location, keyword)
        staging = mktempdir(location; prefix=".$keyword.download-", cleanup=false)
        try
            @info "Downloading '$keyword' ($(DATABASE_SOURCES[keyword].description)) to $destination"
            _download_one(keyword, staging, get(urls, keyword, DATABASE_SOURCES[keyword].urls))
            write(joinpath(staging, DATABASE_MARKER), keyword)
            _validate_database(staging)
            _replace_directory(staging, destination)
        finally
            rm(staging; force=true, recursive=true)
        end
        push!(directories, destination)
    end
    return directories
end

"""
    _validate_database(directory)

Minimal content check of a prepared database: it loads with
[`load_database`](@ref), its taxonomy is not empty, the first record of its FASTA
parses with a non-empty sequence, and that record has a taxonomy entry. Throws an
`ArgumentError` otherwise. The marker only types the directory; it is not a
certificate of integrity.
"""
function _validate_database(directory)
    db = try
        load_database(directory)
    catch err
        throw(ArgumentError("downloaded database in $directory cannot be loaded: " * sprint(showerror, err)))
    end
    isempty(db.taxonomy) && throw(ArgumentError("downloaded database in $directory has no taxonomy entries"))
    header = try
        _first_fasta_record(db.fasta_path)
    catch err
        throw(ArgumentError("downloaded FASTA $(db.fasta_path) is not valid: " * sprint(showerror, err)))
    end
    header === nothing && throw(ArgumentError("downloaded FASTA $(db.fasta_path) has no sequence"))
    haskey(db.taxonomy, database_key(db, header)) ||
        throw(ArgumentError("the first reference of $(db.fasta_path) ($(repr(header))) has no taxonomy entry"))
    return db
end

# Description of the first FASTA record (plain or gzip), `nothing` without one;
# a record without sequence is an error.
function _first_fasta_record(path)
    open(path) do raw
        stream = !eof(raw) && peek(raw) == 0x1f ? GzipDecompressorStream(raw) : raw
        reader = FASTX.FASTA.Reader(stream)
        try
            for record in reader
                isempty(FASTX.sequence(String, record)) && throw(ArgumentError("first record has an empty sequence"))
                return String(FASTX.description(record))
            end
            return nothing
        finally
            close(reader)
        end
    end
end

# Put the complete directory `staging` at `destination`, keeping the previous
# `destination` until the new one is in place (restored if the move fails). Plain
# renames, in the same parent directory: a rename that fails changes nothing.
function _replace_directory(staging, destination)
    ispath(destination) || return _rename(staging, destination)
    backup = destination * ".previous-" * basename(staging)
    _rename(destination, backup)
    try
        _rename(staging, destination)
    catch
        try
            _rename(backup, destination)
        catch
            @error "Could not restore the previous database; it is kept in $backup"
        end
        rethrow()
    end
    try rm(backup; force=true, recursive=true) catch end
    return destination
end

# Rename without `Base.mv`'s copy-and-delete fallback (libuv `uv_fs_rename`, as
# `Base.Filesystem.rename` does from Julia 1.12 on); throws an `IOError` and leaves
# both paths unchanged on failure.
function _rename(source::AbstractString, target::AbstractString)
    status = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), source, target)
    status < 0 && Base.uv_error("rename($(repr(source)), $(repr(target)))", status)
    return target
end
