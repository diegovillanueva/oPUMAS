# 📐🧮 Real powers, exp and log from the system libm, as gfortran calls them, so that the
#    twin rounds like the Fortran (FINDINGS.md). Integer and Float32 powers stay Julia's.

const LIBM = "libm.so.6"

^(x::Float64, y::Float64) = ccall((:pow, LIBM), Float64, (Float64, Float64), x, y)
^(x, y) = Base.:^(x, y)
exp(x::Float64) = ccall((:exp, LIBM), Float64, (Float64,), x)
exp(x) = Base.exp(x)
log(x::Float64) = ccall((:log, LIBM), Float64, (Float64,), x)
log(x) = Base.log(x)
log10(x::Float64) = ccall((:log10, LIBM), Float64, (Float64,), x)
log10(x) = Base.log10(x)
