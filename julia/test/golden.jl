# 📦🧪 Access to the golden files written by tests/column/run.sh.

using NCDatasets

"Directory of the golden files (`OPUMAS_GOLDEN`, default `tests/column/data`)."
const GOLDEN = get(
    ENV, "OPUMAS_GOLDEN", joinpath(@__DIR__, "..", "..", "tests", "column", "data")
)

function golden_path(name)
    path = joinpath(GOLDEN, name)
    isfile(path) || error("missing $path: run tests/column/run.sh first")
    return path
end

"Read one captured variable: a scalar, or an array in (col, lev) order."
function captured(ds, name)
    v = ds[name]
    return ndims(v) == 0 ? v[] : Array(v)
end

"""
Relative error of `a` against reference `b`, zero where both are zero. Exact zeros in the
reference must stay exact: a non-zero `a` there counts as an infinite error.
"""
function relative_error(a, b)
    b == 0 && return a == 0 ? 0.0 : Inf
    return abs(a - b) / abs(b)
end
