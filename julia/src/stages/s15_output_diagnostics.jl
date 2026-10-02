# 📡🌧️ Stage 15, output diagnostics (micro_pumas_v1.F90 lines 4307–4547): mean
#    diameters and effective radii of precipitation, radar reflectivities and the ice
#    mass fraction.

"""
    output_diagnostics!(work, parameters)

📡 Stage 15. Fills the `*out2`, `d*out*`, `freq*`, `reff_*`, radar and `nfice`
diagnostics from the final state.
"""
function output_diagnostics!(work, parameters)
    settings = parameters.settings
    rhogtmp = dense_precipitation(parameters).density
    F = work.fortran
    (; qc, nc, qi, rho, lcldm, icldm, precip_frac, nctend, rflx) = F
    (; qrout, nrout, qsout, nsout, qgout, ngout, qrout2, nrout2, qsout2, nsout2) = F
    (; qgout2, ngout2, drout2, dsout2, dgout2, dsout, dgout, freqr, freqs, freqg) = F
    (; reff_rain, reff_snow, reff_grau, refl, arefl, areflz, frefl, csrfl, acsrfl) = F
    (; fcsrfl, refl10cm, reflz10cm, nfice, lamr, n0r, lami, lams, lamg, n0g) = F
    (; dumni0A2D, dumns0A2D, dum_2D, dumi) = F
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)

    # 🌧️ rain
    for k in 1:levels, i in 1:columns
        drout2[i, k] = mean_diameter(qrout[i, k], nrout[i, k], rho[i, k], rhow)
    end
    for k in 1:levels, i in 1:columns
        if qrout[i, k] > 1.0e-7 && nrout[i, k] > 0.0
            qrout2[i, k] = qrout[i, k] * precip_frac[i, k]
            nrout2[i, k] = nrout[i, k] * precip_frac[i, k]
            freqr[i, k] = precip_frac[i, k]
            reff_rain[i, k] = 1.5 * drout2[i, k] * 1.0e6
        else
            qrout2[i, k] = 0.0
            nrout2[i, k] = 0.0
            drout2[i, k] = 0.0
            freqr[i, k] = 0.0
            reff_rain[i, k] = 0.0
        end
        if nsout[i, k] == 0.0
            nsout[i, k] = 1.0e-34
        end
    end

    # ❄️ snow
    for k in 1:levels, i in 1:columns
        dsout2[i, k] = mean_diameter(qsout[i, k], nsout[i, k], rho[i, k], rhosn)
    end
    for k in 1:levels, i in 1:columns
        if qsout[i, k] > 1.0e-7 && nsout[i, k] > 0.0
            qsout2[i, k] = qsout[i, k] * precip_frac[i, k]
            nsout2[i, k] = nsout[i, k] * precip_frac[i, k]
            freqs[i, k] = precip_frac[i, k]
            dsout[i, k] = 3.0 * rhosn / rhows * dsout2[i, k]
            reff_snow[i, k] = 1.5 * dsout2[i, k] * 1.0e6
        else
            dsout[i, k] = 0.0
            qsout2[i, k] = 0.0
            nsout2[i, k] = 0.0
            dsout2[i, k] = 0.0
            freqs[i, k] = 0.0
            reff_snow[i, k] = 0.0
        end
        if ngout[i, k] == 0.0
            ngout[i, k] = 1.0e-34
        end
    end

    # 🟤 graupel or hail
    dense_on = settings.do_hail || settings.do_graupel
    for k in 1:levels, i in 1:columns
        dgout2[i, k] = dense_on ?
                       mean_diameter(qgout[i, k], ngout[i, k], rho[i, k], rhogtmp) : 0.0
    end
    for k in 1:levels, i in 1:columns
        if qgout[i, k] > 1.0e-7 && ngout[i, k] > 0.0
            qgout2[i, k] = qgout[i, k] * precip_frac[i, k]
            ngout2[i, k] = ngout[i, k] * precip_frac[i, k]
            freqg[i, k] = precip_frac[i, k]
            dgout[i, k] = 3.0 * rhogtmp / rhows * dgout2[i, k]
            reff_grau[i, k] = 1.5 * dgout2[i, k] * 1.0e6
        else
            dgout[i, k] = 0.0
            qgout2[i, k] = 0.0
            ngout2[i, k] = 0.0
            dgout2[i, k] = 0.0
            freqg[i, k] = 0.0
            reff_grau[i, k] = 0.0
        end
    end

    # 📡 analytic radar reflectivity (Shupe), in-precipitation, mm⁶ m⁻³
    rthrsh = 0.0001 / 3600.0
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= qsmall && (nc[i, k] + nctend[i, k] * deltat) > 10.0
            x = qc[i, k] / lcldm[i, k] * rho[i, k] * 1000.0
            dum = (x * x) /
                  (0.109 * (nc[i, k] + nctend[i, k] * deltat) / lcldm[i, k] * rho[i, k] /
                   1.0e6) * lcldm[i, k] / precip_frac[i, k]
        else
            dum = 0.0
        end
        if qi[i, k] >= qsmall
            dum1 = (qi[i, k] * rho[i, k] / icldm[i, k] * 1000.0 / 0.1)^(1.0 / 0.63) *
                   icldm[i, k] / precip_frac[i, k]
        else
            dum1 = 0.0
        end
        if qsout[i, k] >= qsmall
            dum1 = dum1 + (qsout[i, k] * rho[i, k] * 1000.0 / 0.1)^(1.0 / 0.63)
        end
        refl[i, k] = dum + dum1
        # 🌧️ Z = a R^b, R in mm h⁻¹ (Comstock et al. 2004)
        if rflx[i, k+1] >= rthrsh
            dum = 32.0 * (rflx[i, k+1] * 3600.0)^1.4
        else
            dum = 0.0
        end
        refl[i, k] = refl[i, k] + dum
        areflz[i, k] = refl[i, k] * precip_frac[i, k]
        if refl[i, k] > radar.minimum_reflectivity
            refl[i, k] = 10.0 * log10(refl[i, k])
        else
            refl[i, k] = -9999.0
        end
        if refl[i, k] > radar.minimum_dbz
            arefl[i, k] = refl[i, k] * precip_frac[i, k]
            frefl[i, k] = precip_frac[i, k]
        else
            arefl[i, k] = 0.0
            areflz[i, k] = 0.0
            frefl[i, k] = 0.0
        end
        csrfl[i, k] = min(radar.cloudsat_max_dbz, refl[i, k])
        if csrfl[i, k] > radar.cloudsat_min_dbz
            acsrfl[i, k] = refl[i, k] * precip_frac[i, k]
            fcsrfl[i, k] = precip_frac[i, k]
        else
            acsrfl[i, k] = 0.0
            fcsrfl[i, k] = 0.0
        end
    end

    # 📡 10 cm reflectivity (Morrison; ice dielectric factor from Smith 1984)
    minrefl10 = radar.minimum_reflectivity_10cm
    for k in 1:levels, i in 1:columns
        dum1 = minrefl10
        dum2 = minrefl10
        dum3 = minrefl10
        dum4 = minrefl10
        if lamr[i, k] > 0.0
            x = lamr[i, k]
            dum1 = rho[i, k] * n0r[i, k] * 720.0 / (x * x * x) / (x * x * x) / x
            dum1 = max(dum1, minrefl10)
        end
        if lami[i, k] > 0.0
            x = lami[i, k]
            x2 = x * x
            dum2 = rho[i, k] * (0.176 / 0.93) * 720.0 * dumni0A2D[i, k] *
                   ((rhoi / 900.0) * (rhoi / 900.0)) / ((x * x2) * (x2 * x2))
            dum2 = max(dum2, minrefl10)
        end
        if lams[i, k] > 0.0
            dum3 = rho[i, k] * (0.176 / 0.93) * 720.0 * dumns0A2D[i, k] *
                   ((rhosn / 900.0) * (rhosn / 900.0)) / lams[i, k]^7.0
            dum3 = max(dum3, minrefl10)
        end
        # 🟤 `do_hail .or. do_graupel .and. lamg > 0` binds as hail or (graupel and λ > 0)
        if settings.do_hail || (settings.do_graupel && lamg[i, k] > 0.0)
            dum4 = rho[i, k] * (0.176 / 0.93) * 720.0 * n0g[i, k] *
                   ((rhogtmp / 900.0) * (rhogtmp / 900.0)) / lamg[i, k]^7.0
            dum4 = max(dum4, minrefl10)
        end
        reflz10cm[i, k] = (dum1 + dum2 + dum3 + dum4) * precip_frac[i, k]
        dum = reflz10cm[i, k] * 1.0e18
        refl10cm[i, k] = 10.0 * log10(dum)

        # 🧊 ice mass fraction
        dum_2D[i, k] = qsout[i, k] + qrout[i, k] + qc[i, k] + qi[i, k]
        dumi[i, k] = qsout[i, k] + qi[i, k]
        if dumi[i, k] > qsmall && dum_2D[i, k] > qsmall
            nfice[i, k] = min(dumi[i, k] / dum_2D[i, k], 1.0)
        else
            nfice[i, k] = 0.0
        end
    end
    return nothing
end
