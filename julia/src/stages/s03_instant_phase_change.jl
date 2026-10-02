# 🌡️🔁 Stage 3, instantaneous melting and freezing, first pass (micro_pumas_v1.F90
#    lines 1784–1898). The state is changed directly; stage 10 moves it into tendencies.

"""
    instant_phase_change!(work, parameters)

🌡️ Stage 3. Snow and graupel melt above +2 °C and rain freezes below the rain freezing
temperature, each limited so that the latent heat does not move the temperature across
the threshold.
"""
function instant_phase_change!(work, parameters)
    (; settings) = parameters
    (; cpp, xlf, snowmelt, rainfrze) = host_short_names(parameters)
    F = work.fortran
    (; t, qr, nr, qs, ns, qg, ng, tlat) = F
    (; minstsm, ninstsm, minstgm, ninstgm, minstrf, ninstrf) = F
    (; meltsdttot, meltstot, meltgtot, frzrdttot) = F
    rdeltat = 1.0 / work.inputs.time_step
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        # ❄️➡️💧 melting of snow at +2 °C
        if t[i, k] > snowmelt && qs[i, k] > 0.0
            dum = -xlf / cpp * qs[i, k]
            if t[i, k] + dum < snowmelt
                dum = (t[i, k] - snowmelt) * cpp / xlf
                dum = dum / qs[i, k]
                dum = max(0.0, dum)
                dum = min(1.0, dum)
            else
                dum = 1.0
            end
            minstsm[i, k] = dum * qs[i, k]
            ninstsm[i, k] = dum * ns[i, k]
            dum1 = -xlf * minstsm[i, k] * rdeltat
            tlat[i, k] = tlat[i, k] + dum1
            meltsdttot[i, k] = meltsdttot[i, k] + dum1
            meltstot[i, k] = minstsm[i, k] * rdeltat
            qs[i, k] = max(qs[i, k] - minstsm[i, k], 0.0)
            ns[i, k] = max(ns[i, k] - ninstsm[i, k], 0.0)
            qr[i, k] = max(qr[i, k] + minstsm[i, k], 0.0)
            nr[i, k] = max(nr[i, k] + ninstsm[i, k], 0.0)
        end

        # 🟤➡️💧 melting of graupel at +2 °C
        if t[i, k] > snowmelt && qg[i, k] > 0.0
            dum = -xlf / cpp * qg[i, k]
            if t[i, k] + dum < snowmelt
                dum = (t[i, k] - snowmelt) * cpp / xlf
                dum = dum / qg[i, k]
                dum = max(0.0, dum)
                dum = min(1.0, dum)
            else
                dum = 1.0
            end
            minstgm[i, k] = dum * qg[i, k]
            ninstgm[i, k] = dum * ng[i, k]
            dum1 = -xlf * minstgm[i, k] * rdeltat
            tlat[i, k] = tlat[i, k] + dum1
            meltsdttot[i, k] = meltsdttot[i, k] + dum1
            meltgtot[i, k] = minstgm[i, k] * rdeltat
            qg[i, k] = max(qg[i, k] - minstgm[i, k], 0.0)
            ng[i, k] = max(ng[i, k] - ninstgm[i, k], 0.0)
            qr[i, k] = max(qr[i, k] + minstgm[i, k], 0.0)
            nr[i, k] = max(nr[i, k] + ninstgm[i, k], 0.0)
        end

        # 💧➡️🧊 freezing of rain, to graupel (or hail) if on, else to snow
        if t[i, k] < rainfrze && qr[i, k] > 0.0
            dum = xlf / cpp * qr[i, k]
            if t[i, k] + dum > rainfrze
                dum = -(t[i, k] - rainfrze) * cpp / xlf
                dum = dum / qr[i, k]
                dum = max(0.0, dum)
                dum = min(1.0, dum)
            else
                dum = 1.0
            end
            minstrf[i, k] = dum * qr[i, k]
            ninstrf[i, k] = dum * nr[i, k]
            dum1 = xlf * minstrf[i, k] * rdeltat
            tlat[i, k] = tlat[i, k] + dum1
            frzrdttot[i, k] = frzrdttot[i, k] + dum1
            qr[i, k] = max(qr[i, k] - minstrf[i, k], 0.0)
            nr[i, k] = max(nr[i, k] - ninstrf[i, k], 0.0)
            if settings.do_hail || settings.do_graupel
                qg[i, k] = max(qg[i, k] + minstrf[i, k], 0.0)
                ng[i, k] = max(ng[i, k] + ninstrf[i, k], 0.0)
            else
                qs[i, k] = max(qs[i, k] + minstrf[i, k], 0.0)
                ns[i, k] = max(ns[i, k] + ninstrf[i, k], 0.0)
            end
        end
    end
    return nothing
end
