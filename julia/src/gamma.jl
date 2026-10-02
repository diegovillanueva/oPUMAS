# 🧮🔢 The gamma function, taken from the same C library the Fortran uses.

"""
    gamma_function(x)

Euler's gamma function Γ(x) for real `x`.

🧊🔗 With gfortran, `pumas_gamma` is the compiler intrinsic, which calls the C library's
`tgamma`; this function calls the same `tgamma`, so both languages get identical bits.

Fortran twin: `pumas_gamma` in `pumas_gamma_function.F90`.
"""
gamma_function(x::Float64) = ccall((:tgamma, "libm.so.6"), Float64, (Float64,), x)
