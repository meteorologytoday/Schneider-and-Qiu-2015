function step!(
    m :: Model,
    dt :: Float64,
)

    st = m.st
    

 
    A_Γ, F_Γ = getΓMatrix(m)
    A_B, F_B = getBMatrix(m)
    
    A = A_Γ + A_B
    F = F_Γ + F_B    
   
#    print(Matrix(A))
     
    
    _X = st._X
    _X[:] = ( I - dt * A ) \ ( _X + dt * F )
    
    diagΨ!(m)

end



function getΓMatrix(
    m :: Model,
    return_jacobian :: Bool = false,
)
    pp = m.ev.pp
    amo = m.co.amo
    zero_vec = zeros(Float64, m.ev.gd.Ny)
    zero_mtx = spdiagm( 0 => zero_vec )
    
    Jk2 = pp.J * pp.k^2

    A_YDiff = blockdiag(m.co.ops[:ydiff], zero_mtx)
    F_YDiff = vcat(zero_vec, zero_vec)

    A_linear = blockdiag(
        - Jk2 * (1 .+ amo.T_f_T * amo.T_f_T / Jk2^2),
        zero_mtx,
    )
 
    F_linear = vcat(zero_vec, zero_vec)
    
    A_baroclinic = sparse(
        [ 
            zero_mtx   amo.T_interp_V * (- amo.V_f_V / Jk2) * amo.V_∂y_T ;
            zero_mtx   zero_mtx  ;
        ]
    )
        
    #print(Matrix(amo.T_interp_V * (amo.V_f_V / Jk2) * amo.V_∂y_T))

    F_baroclinic = vcat(zero_vec, zero_vec)
    
    A = A_YDiff + A_linear + A_baroclinic
    F = F_YDiff + F_linear + F_baroclinic
   
    return A, F 
    
end



function getBMatrix(
    m :: Model,
    return_jacobian :: Bool = false,
)
    
    Ny = m.ev.gd.Ny
    zero_vec = zeros(Float64, Ny)
    zero_mtx = spdiagm( 0 => zero_vec )

    A_DIVΨ = - m.ev.pp.N * sparse( vcat( spzeros(Ny, 2*Ny), m.co.amo.T_DIVy_V * m.co.ops[:M_Ψ]) )
    F_DIVΨ = vcat(zero_vec, zero_vec)
   
    A_F_TOA, F_F_TOA = B_getF_TOAMatrix(m)
    A_F_AO,  F_F_AO  = B_getF_AOMatrix(m)
    
    A_YDiff = blockdiag(zero_mtx, m.co.ops[:ydiff])
    F_YDiff = vcat(zero_vec, zero_vec)
  
    A_F_TOA = blockdiag(zero_mtx, A_F_TOA) 
    A_F_AO  = blockdiag(zero_mtx, A_F_AO)

    F_F_TOA = vcat(zero_vec, F_F_TOA)
    F_F_AO  = vcat(zero_vec, F_F_AO)
 
    
    A = A_F_TOA + A_F_AO + A_YDiff + A_DIVΨ
    F = F_F_TOA + F_F_AO + F_YDiff + F_DIVΨ
    
    return A, F
    
end




function diagΨ!(
    m :: Model,
)

    m.st.Ψ[:] = m.co.ops[:M_Ψ] * m.st._X
    m.st.Ψ_wgt[:] = m.co.amo.V_Δx_V * m.st.Ψ

end


# 

function B_getYDiffMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    A = m.co.ops[:ydiff]
    F = m.st.B * 0

    jacobian = m.co.ops[:ydiff]


    if return_jacobian

        return A, F, jacobian
    else

        return A, F

    end

end

function computeTS(
    m :: Model,
    B :: AbstractArray,
)

    pp = m.ev.pp
    gd = m.ev.gd


    return θ0 * (1 .+ (B .- (pp.N * gd.H) / 2) / g0)

end



function computeT_TOA(
    m :: Model,
    B :: AbstractArray,
)

    pp = m.ev.pp
    gd = m.ev.gd


    return θ0 * (1 .+ (B .+ (pp.N * gd.H) / 2) / g0)

end


function B_getF_TOAMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    ev = m.ev
    pp = ev.pp
    gd = ev.gd

    Npts = gd.Ny

    factor = g0 / ( θ0 * c_p * ρ0 * gd.H ) * (gd.H * pp.k / 2) 

    T_TOA = computeT_TOA(m, m.st.B)
    tmp = σ_boltz * T_TOA.^3
    A   = spdiagm(0 => (- factor * tmp * θ0 / g0))
    F   = - factor * tmp * ( θ0 * ( 1 + pp.N * gd.H / (2 * g0)) )  

    jacobian_A = A
    jacobian_F = factor * spdiagm( 0 => - 3 * σ_boltz * T_TOA.^2 * ( θ0 * ( 1 + pp.N * gd.H / (2 * g0) ) ) * (θ0/g0) )


    if return_jacobian
        
        return A, F, jacobian_A + jacobian_F
        
    else
        return A, F
    end
end

function B_getF_AOMatrix(
    m :: Model,
    return_jacobian :: Bool = false
)

    ev = m.ev
    pp = ev.pp
    gd = ev.gd
    st = m.st

    Npts = gd.Ny
    
    factor = g0 / ( θ0 * c_p * ρ0 * gd.H ) * (gd.H * pp.k / 2) 

    # F_ao
    A  = - factor * pp.C * θ0 / g0 * sparse(I, Npts, Npts)
    F = factor * ( st.SST .- θ0 * ( 1 - pp.N * gd.H / (2 * g0)) ) * pp.C 
    #println(Matrix(A))
    jacobian_A = A
    jacobian_F = 0 * A


    if return_jacobian
        
        return A, F, jacobian_A + jacobian_F
        
    else
        return A, F
    end
    
end
