function getFlux_TOA(
    T :: Float64;
    T_ref :: Union{Float64, Nothing} = nothing,
)
    if T_ref == nothing
        return σ_boltz * T^4
    else
        return 4 * σ_boltz * T_ref^3 * (T - T_ref)
    end

end


function getFlux_atmocn(
    
)

end
