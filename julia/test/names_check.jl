# 🔤🔗 Checks that keep translation/dictionary.csv and oPUMAS.jl in sync, in both directions:
#    1. the header is the agreed one, ids are unique, and no Fortran, Julia or meaning cell is empty;
#    2. every Julia path exists (a workspace field, or a `Settings` field for `settings.*`);
#    3. every entry of the name tables in src/names.jl and of NAMELIST_FIELDS has its row;
#    4. every Fortran name still occurs in the Fortran sources (catches upstream renames);
#    5. every exported function of TWINS.csv names its Fortran routine in a `Fortran twin:`
#       (or `Twin of`) docstring line, and every other export is a listed package utility.
#    Used by test/test_names.jl (Pkg.test) and tools/check_names.jl.

const REPO = dirname(dirname(@__DIR__))
const DICTIONARY = joinpath(REPO, "translation", "dictionary.csv")
const HEADER = ["id", "kind", "scope", "units", "meaning", "Fortran", "Julia", "standard_name", "notes"]

"Fields of one CSV line, honouring double quotes around fields with commas."
function csv_fields(line)
    fields, field, quoted = String[], IOBuffer(), false
    for c in line
        if c == '"'
            quoted = !quoted
        elseif c == ',' && !quoted
            push!(fields, String(take!(field)))
        else
            print(field, c)
        end
    end
    push!(fields, String(take!(field)))
    return fields
end

"Rows of a CSV file as Dicts keyed by the header."
function read_rows(path)
    lines = readlines(path)
    header = csv_fields(first(lines))
    return [Dict(zip(header, csv_fields(line))) for line in lines[2:end] if !isempty(line)]
end

"Exports that are package plumbing, not twins of a Fortran routine."
const UTILITIES = (:ColumnInputs, :read_inputs, :allocate, :output_paths, :field_at,
    :read_settings, :ventilated_diffusion, :SaturationConstants, :HydrometeorProperties,
    :HostConstants, :Settings, :Parameters)

"Every (Fortran name, Julia path) pair that oPUMAS.jl itself defines."
function julia_names()
    pairs = Set{Tuple{String, String}}()
    for (table, container) in ((oPUMAS.LOCAL_RATE_PATHS, "rates"),
        (oPUMAS.PROCESS_RATE_PATHS, "process_rates"), (oPUMAS.TENDENCY_PATHS, "tendencies"),
        (oPUMAS.DIAGNOSTIC_PATHS, "diagnostics"), (oPUMAS.FLUX_PATHS, "diagnostics"),
        (oPUMAS.SURFACE_PATHS, "diagnostics"))
        for (name, path) in Base.pairs(table)
            push!(pairs, (String(name), join((container, String.(path)...), '.')))
        end
    end
    for (name, field) in oPUMAS.NAMELIST_FIELDS
        push!(pairs, (name, "settings.$field"))
    end
    return pairs
end

"Lower-case identifiers of every Fortran file of the repo (sources, interfaces, column driver)."
function fortran_identifiers()
    words = Set{String}()
    for (dir, _, files) in walkdir(REPO)
        occursin(r"/(\.git|julia|build|data)(/|$)", dir) && continue
        for file in files
            endswith(file, ".F90") || continue
            path = joinpath(dir, file)
            islink(path) && continue
            for m in eachmatch(r"[A-Za-z_]\w*", read(path, String))
                push!(words, lowercase(m.match))
            end
        end
    end
    return words
end

"All errors of the five checks, empty when dictionary and twin agree."
function name_errors()
    errors = String[]
    header = csv_fields(first(readlines(DICTIONARY)))
    header == HEADER || push!(errors, "dictionary header $header, expected $HEADER")
    dictionary = read_rows(DICTIONARY)

    # 🧾 1. format
    ids = [row["id"] for row in dictionary]
    for id in unique(filter(id -> count(==(id), ids) > 1, ids))
        push!(errors, "dictionary: duplicate id $id")
    end
    for row in dictionary, column in ("Fortran", "Julia", "meaning")
        isempty(row[column]) && push!(errors, "dictionary $(row["id"]): empty $column")
    end

    # 🗂️ 2. Julia paths
    work = allocate((; T = zeros(1, 1)))
    for row in dictionary
        path = row["Julia"]
        parts = Symbol.(split(path, '.'))
        if first(parts) == :settings
            length(parts) == 2 && hasfield(Settings, parts[2]) ||
                push!(errors, "dictionary $(row["id"]): no Settings field $path")
            continue
        end
        node = work
        for part in parts
            if node isa NamedTuple && hasproperty(node, part)
                node = getproperty(node, part)
            else
                push!(errors, "dictionary $(row["id"]): no field $path")
                break
            end
        end
    end

    # 🔁 3. names.jl and namelist.jl entries have their rows
    rows = Set((row["Fortran"], row["Julia"]) for row in dictionary)
    for (name, path) in sort!(collect(julia_names()))
        (name, path) in rows || push!(errors, "dictionary: no row for $name => $path")
    end

    # 🧱 4. Fortran names
    words = fortran_identifiers()
    for row in dictionary
        lowercase(row["Fortran"]) in words ||
            push!(errors, "dictionary $(row["id"]): Fortran name $(row["Fortran"]) not in the sources")
    end

    # 📜 5. twin docstrings
    twins = read_rows(joinpath(REPO, "julia", "TWINS.csv"))
    routines = Dict{Symbol, Vector{String}}()
    for row in twins
        push!(get!(routines, Symbol(row["julia_function"]), String[]), row["fortran_routine"])
    end
    for name in filter(name -> name != :oPUMAS, names(oPUMAS))
        if haskey(routines, name)
            # 🧾 raw docstring text: rendering Markdown would eat underscores (`ρ_air`)
            docs = Base.Docs.meta(oPUMAS)[Base.Docs.Binding(oPUMAS, name)].docs
            text = replace(join((join(d.text) for d in values(docs)), ' '), r"\s+" => " ")
            occursin("Fortran twin", text) || occursin("Twin of", text) ||
                push!(errors, "$name: no `Fortran twin` docstring line")
            for routine in routines[name]
                routine = replace(first(split(routine, ' ')), r"_(vect|vec)$" => "")
                occursin(routine, text) ||
                    push!(errors, "$name: docstring does not name its twin $routine")
            end
        elseif !(name in UTILITIES)
            push!(errors, "$name: exported but neither in TWINS.csv nor a listed utility")
        end
    end
    for name in keys(routines)
        isdefined(oPUMAS, name) || push!(errors, "TWINS.csv: $name is not defined in oPUMAS")
    end
    return errors
end
