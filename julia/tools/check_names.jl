# 🔤🔗 Check that translation/dictionary.csv, oPUMAS.jl and the Fortran agree, outside Pkg.test
#    (the checks themselves are in test/names_check.jl, plan M7).
#    Usage: julia/tools/julia.sh --project=julia julia/tools/check_names.jl
using oPUMAS

include(joinpath(dirname(@__DIR__), "test", "names_check.jl"))

errors = name_errors()
foreach(e -> println("❌ ", e), errors)
println(isempty(errors) ? "✅" : "❌", " $(length(read_rows(DICTIONARY))) dictionary rows, ",
    "$(length(errors)) error(s)")
exit(isempty(errors) ? 0 : 1)
