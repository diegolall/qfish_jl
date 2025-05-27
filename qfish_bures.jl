using LinearAlgebra

function qfish_bures(ρ::AbstractArray,ρ_end::AbstractArray,dtheta::BigFloat)
    fidelity=real(tr(√(ρ*ρ_end))^2)
    bures_dist=2*(1-sqrt(fidelity))
    return 4*bures_dist/(dtheta^2)
end #function