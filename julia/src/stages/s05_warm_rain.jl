# 🌧️🔁 Stage 5, autoconversion of cloud liquid to rain (micro_pumas_v1.F90 lines
#    2028–2165). KK2000 is always computed as a diagnostic.

"""
    autoconversion!(work, parameters)

🌧️ Stage 5. Autoconversion by `kk2000` or `sb2001`, chosen by `warm_rain`.
"""
function autoconversion!(work, parameters)
    (; settings) = parameters
    F = work.fortran
    (; qcic, ncic, qric, rho, relvar, pgam, prc, nprc, nprc1) = F
    (; qctend_KK2000, nrtend_KK2000, nctend_KK2000, qrtend_KK2000) = F
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        kk = kk2000_autoconversion(parameters, qcic[i, k], ncic[i, k], rho[i, k],
            relvar[i, k])
        qctend_KK2000[i, k] = kk.mass
        nrtend_KK2000[i, k] = kk.rain_number
        nctend_KK2000[i, k] = kk.droplet_number
    end
    if settings.warm_rain == :kk2000
        for k in 1:levels, i in 1:columns
            prc[i, k] = qctend_KK2000[i, k]
            nprc1[i, k] = nctend_KK2000[i, k]
            nprc[i, k] = nrtend_KK2000[i, k]
            qrtend_KK2000[i, k] = -qctend_KK2000[i, k]
        end
    elseif settings.warm_rain == :sb2001
        (; qctend_SB2001, nrtend_SB2001, nctend_SB2001, qrtend_SB2001) = F
        for k in 1:levels, i in 1:columns
            sb = sb2001_autoconversion(pgam[i, k], qcic[i, k], ncic[i, k], qric[i, k],
                rho[i, k])
            qctend_SB2001[i, k] = sb.mass
            nrtend_SB2001[i, k] = sb.rain_number
            nctend_SB2001[i, k] = sb.droplet_number
        end
        for k in 1:levels, i in 1:columns
            prc[i, k] = qctend_SB2001[i, k]
            nprc[i, k] = nrtend_SB2001[i, k]
            nprc1[i, k] = nctend_SB2001[i, k]
            qrtend_SB2001[i, k] = -qctend_SB2001[i, k]
        end
    end
    return nothing
end
