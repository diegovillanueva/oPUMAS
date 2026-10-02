# 🔤🔗 Check the names that tie oPUMAS.jl to the Fortran (plan M7):
#    1. every `Julia` path of translation/dictionary.csv exists in the workspace;
#    2. every exported function of TWINS.csv names its Fortran routine in a
#       `Fortran twin:` (or `Twin of`) docstring line, and every other export is a listed
#       package utility.
#    Usage: julia/tools/julia.sh --project=julia julia/tools/check_names.jl
using oPUMAS

const REPO = dirname(dirname(@__DIR__))

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

errors = String[]

# 🗂️ 1. dictionary paths
work = allocate((; T = zeros(1, 1)))
dictionary = read_rows(joinpath(REPO, "translation", "dictionary.csv"))
haskey(first(dictionary), "Julia") ||
    error("translation/dictionary.csv has no Julia column")
paths = 0
for row in dictionary
    path = row["Julia"]
    isempty(path) && continue
    global paths += 1
    node = work
    for part in Symbol.(split(path, '.'))
        if node isa NamedTuple && hasproperty(node, part)
            node = getproperty(node, part)
        else
            push!(errors, "dictionary $(row["id"]): no field $path")
            break
        end
    end
end

# 🔁 2. twin docstrings
twins = read_rows(joinpath(REPO, "julia", "TWINS.csv"))
routines = Dict{Symbol, Vector{String}}()
for row in twins
    push!(get!(routines, Symbol(row["julia_function"]), String[]), row["fortran_routine"])
end
exported = filter(name -> name != :oPUMAS, names(oPUMAS))
for name in exported
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

foreach(e -> println("❌ ", e), errors)
println(isempty(errors) ? "✅" : "❌", " $paths dictionary paths, ",
    count(in(keys(routines)), exported), " exported twins, $(length(errors)) error(s)")
exit(isempty(errors) ? 0 : 1)
