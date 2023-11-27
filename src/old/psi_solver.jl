mutable struct StreamfunctionSolver

    amo          :: AdvancedMatrixOperators

    eW_send_W    :: AbstractArray{Float64, 2}
    W_send_eW    :: AbstractArray{Float64, 2}

    V_solveΨ_fromB_T    :: AbstractArray{Float64, 2}
    V_solveΨ_from∂B∂y_V :: AbstractArray{Float64, 2}
 
    T_solveΓ_fromB_T    :: AbstractArray{Float64, 2}
    T_solveΓ_from∂B∂y_V :: AbstractArray{Float64, 2}
    
    function StreamfunctionSolver(;
        pp :: PhyParams,
        amo :: AdvancedMatrixOperators,
    )
    
        bmo = amo.bmo
        Ny = bmo.Ny

        # Create eV
        # need a mask excluding bnd points
        mask_eff_V = reshape(amo.V_mask_V * ones(Float64, bmo.V_pts), bmo.V_dim...)

        # Create coversion matrix and its inverse
        V_num = reshape(collect(1:length(mask_eff_V)), size(mask_eff_V)...)
        active_num_eff_V = V_num[ mask_eff_V .==1 ]

        eV_send_V = bmo.V_I_V[active_num_eff_V, :]; dropzeros!(eV_send_V)
        V_send_eV = sparse(eV_send_V')
 
        # construct operators on the left-hand-side
        T_WgtLAPy_T = sparse(amo.T_DIVy_V * spdiagm(0=>pp.L) * amo.V_∂y_T )
        T_L_T = sparse(
            pp.J * pp.k^2 * T_WgtLAPy_T - ( amo.T_f_T * amo.T_f_T .+ pp.J^2 * pp.k^4 )
        )
        
        T_invL_T = inv(Matrix(T_L_T))
        
        T_solveΓ_fromB_T    = sparse( T_invL_T * amo.T_interp_V * amo.V_f_V * amo.V_∂y_T )
        T_solveΓ_from∂B∂y_V = sparse( T_invL_T * amo.T_interp_V * amo.V_f_V )

        V_solveΨ_fromB_T = - amo.V_mask_V * (pp.J * pp.k^4)^(-1) * (
              amo.V_∂y_T
            + amo.V_interp_T * amo.T_f_T * T_solveΓ_fromB_T
        )


        V_solveΨ_from∂B∂y_V = - amo.V_mask_V * (pp.J * pp.k^4)^(-1) * (
              amo.bmo.V_I_V
            + amo.V_interp_T * amo.T_f_T * T_solveΓ_from∂B∂y_V
        )
        
        return new(
            amo,

            eV_send_V,
            V_send_eV,

            V_solveΨ_fromB_T,
            V_solveΨ_from∂B∂y_V,

            T_solveΓ_fromB_T,
            T_solveΓ_from∂B∂y_V,
        ) 

    end
end

