# ♻️🔁 Stage 10, re-apply the instantaneous adjustments (micro_pumas_v1.F90 lines
#    3172–3220). Stage 2 and 3 changed the state before the process rates; here the state
#    returns to its start-of-step values and the changes move into the tendencies.

"""
    reapply_instant_changes!(work, parameters)

♻️ Stage 10. Droplet activation, snow and graupel melting and rain freezing become
tendencies; snow sublimation and production join `nevapr` and `prain`.
"""
function reapply_instant_changes!(work, parameters)
    F = work.fortran
    (; nc, qs, ns, qr, nr, qg, ng, ncn, qsn, nsn, qrn, nrn, qgr, ngr, npccn, dum_2D) = F
    (; nctend, qstend, nstend, qrtend, nrtend, qgtend, ngtend, nevapr, prain) = F
    rdeltat = 1.0 / work.inputs.time_step
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        nc[i, k] = ncn[i, k]
        nctend[i, k] = nctend[i, k] + npccn[i, k]
        dum_2D[i, k] = qs[i, k]
        qs[i, k] = qsn[i, k]
        qstend[i, k] = qstend[i, k] + (dum_2D[i, k] - qs[i, k]) * rdeltat
        dum_2D[i, k] = ns[i, k]
        ns[i, k] = nsn[i, k]
        nstend[i, k] = nstend[i, k] + (dum_2D[i, k] - ns[i, k]) * rdeltat
        dum_2D[i, k] = qr[i, k]
        qr[i, k] = qrn[i, k]
        qrtend[i, k] = qrtend[i, k] + (dum_2D[i, k] - qr[i, k]) * rdeltat
        dum_2D[i, k] = nr[i, k]
        nr[i, k] = nrn[i, k]
        nrtend[i, k] = nrtend[i, k] + (dum_2D[i, k] - nr[i, k]) * rdeltat
        dum_2D[i, k] = qg[i, k]
        qg[i, k] = qgr[i, k]
        qgtend[i, k] = qgtend[i, k] + (dum_2D[i, k] - qg[i, k]) * rdeltat
        dum_2D[i, k] = ng[i, k]
        ng[i, k] = ngr[i, k]
        ngtend[i, k] = ngtend[i, k] + (dum_2D[i, k] - ng[i, k]) * rdeltat
        nevapr[i, k] = nevapr[i, k] + F.evapsnow[i, k]
        prain[i, k] = prain[i, k] + F.prodsnow[i, k]
    end
    return nothing
end
