# ⏱️📈 Run the offline single-column time loop and record profiles for plotting.
#   tools/julia.sh tools/column_loop.jl settings.nml inputs.nc output.nc [steps]
# 🧾📂 output.nc holds, per step, the state the step started from (PUMAS input names),
#    the tendencies and the surface precipitation, with dimensions (col, lev, step).

using NCDatasets
using oPUMAS

const STATE = (t = (:T,), q = (:q_vapor,), qcn = (:cloud_liquid, :q),
    ncn = (:cloud_liquid, :N), qin = (:cloud_ice, :q), nin = (:cloud_ice, :N),
    qrn = (:rain, :q), nrn = (:rain, :N), qsn = (:snow, :q), nsn = (:snow, :N),
    qgr = (:graupel, :q), ngr = (:graupel, :N))
const RECORDED = (:tlat, :qvlat, :qctend, :qitend, :qrtend, :qstend, :qgtend, :effc,
    :effi, :prect, :preci)

function main(settings_file, inputs_file, output_file, steps)
    parameters = setup(read_settings(settings_file))
    work = allocate(read_inputs(inputs_file))
    paths = Dict(output_paths())
    columns, levels = size(work.inputs.T)
    state = Dict(n => zeros(columns, levels, steps) for n in keys(STATE))
    outputs = Dict(n => zeros(size(field_at(work[paths[n][1]], paths[n][2:end]))..., steps)
                   for n in RECORDED)
    elapsed = @elapsed time_loop!(work, parameters, steps;
        after_step = function (work, step)
            for (name, path) in pairs(STATE)
                state[name][:, :, step] .= field_at(work.inputs, path)
            end
            for name in RECORDED
                (container, path...) = paths[name]
                selectdim(outputs[name], ndims(outputs[name]), step) .=
                    field_at(work[container], path)
            end
        end)
    NCDataset(output_file, "c") do ds
        ds.attrib["deltat"] = work.inputs.time_step
        ds.attrib["settings"] = abspath(settings_file)
        ds.attrib["inputs"] = abspath(inputs_file)
        ds.attrib["note"] = "state at the start of each step; host forcing held fixed"
        p = NCDataset(inputs_file) do source
            Array{Float64}(source["p"].var[:, :])
        end
        defVar(ds, "p", p, ("col", "lev"))
        for (name, values) in state
            defVar(ds, String(name), values, ("col", "lev", "step"))
        end
        for (name, values) in outputs
            dims = ndims(values) == 2 ? ("col", "step") : ("col", "lev", "step")
            defVar(ds, String(name), values, dims)
        end
    end
    println("⏱️ $steps steps of $columns columns in $(round(elapsed; digits = 1)) s")
    println(abspath(output_file))
end

main(ARGS[1], ARGS[2], ARGS[3], length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 144)
