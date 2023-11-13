mutable struct AdvancedMatrixOperators

    bmo       :: BasicMatrixOperators

    T_DIVx_U    :: AbstractArray{Float64, 2}
    T_DIVy_V    :: AbstractArray{Float64, 2}
    T_DIVz_W    :: AbstractArray{Float64, 2}
 
    U_∂x_T      :: AbstractArray{Float64, 2}
    V_∂y_T      :: AbstractArray{Float64, 2}
    W_∂z_T      :: AbstractArray{Float64, 2}

    T_∂x_U      :: AbstractArray{Float64, 2}
    T_∂y_V      :: AbstractArray{Float64, 2}
    T_∂z_W      :: AbstractArray{Float64, 2}
    
    U_∂x_U      :: AbstractArray{Float64, 2}
    V_∂y_V      :: AbstractArray{Float64, 2}
    W_∂z_W      :: AbstractArray{Float64, 2}

    T_LAPx_T      :: AbstractArray{Float64, 2}
    T_LAPy_T      :: AbstractArray{Float64, 2}

    U_interp_T  :: AbstractArray{Float64, 2} 
    V_interp_T  :: AbstractArray{Float64, 2} 
    W_interp_T  :: AbstractArray{Float64, 2} 

    T_interp_U  :: AbstractArray{Float64, 2} 
    T_interp_V  :: AbstractArray{Float64, 2} 
    T_interp_W  :: AbstractArray{Float64, 2} 
 
   
    V_interp_U  :: AbstractArray{Float64, 2} 
    U_interp_V  :: AbstractArray{Float64, 2} 

    W_interp_U  :: AbstractArray{Float64, 2} 
    W_interp_V  :: AbstractArray{Float64, 2} 

    U_interp_W  :: AbstractArray{Float64, 2} 
    V_interp_W  :: AbstractArray{Float64, 2} 
    
    T_NEWSavg_T :: AbstractArray{Float64, 2} 



    T_mask_T       :: AbstractArray{Float64, 2}
    U_mask_U       :: AbstractArray{Float64, 2}
    V_mask_V       :: AbstractArray{Float64, 2}
    W_mask_W       :: AbstractArray{Float64, 2}

    T_Δx_T :: AbstractArray{Float64, 2}
    U_Δx_U :: AbstractArray{Float64, 2}
    V_Δx_V :: AbstractArray{Float64, 2}
    W_Δx_W :: AbstractArray{Float64, 2}

    T_Δy_T :: AbstractArray{Float64, 2}
    U_Δy_U :: AbstractArray{Float64, 2}
    V_Δy_V :: AbstractArray{Float64, 2}
    W_Δy_W :: AbstractArray{Float64, 2}

    T_Δz_T :: AbstractArray{Float64, 2}
    U_Δz_U :: AbstractArray{Float64, 2}
    V_Δz_V :: AbstractArray{Float64, 2}
    W_Δz_W :: AbstractArray{Float64, 2}

    T_Δv_T :: AbstractArray{Float64, 2}


    function AdvancedMatrixOperators(;
        gd :: Grid,
    )


        Nx = gd.Nx
        Ny = gd.Ny
        Nz = gd.Nz
        
        cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))
        
        @time bmo = BasicMatrixOperators(Ny=Ny, Nz=Nz, Nx=Nx)
        
        
        mask_flat = view(gd.mask, :)

        onU_if_unblocked_west_onT = bmo.U_E_T  * mask_flat
        onU_if_unblocked_east_onT = bmo.U_W_T  * mask_flat

        onV_if_unblocked_north_onT = bmo.V_S_T  * mask_flat
        onV_if_unblocked_south_onT = bmo.V_N_T  * mask_flat
 
        onW_if_unblocked_up_onT    = bmo.W_DN_T * mask_flat
        onW_if_unblocked_dn_onT    = bmo.W_UP_T * mask_flat

        U_mask = onU_if_unblocked_east_onT .* onU_if_unblocked_west_onT
        V_mask = onV_if_unblocked_north_onT .* onV_if_unblocked_south_onT
        W_mask = onW_if_unblocked_up_onT    .* onW_if_unblocked_dn_onT

        # For debug
        #U_mask .= 1
        #V_mask .= 1
        #W_mask .= 1

        T_mask_T   = mask_flat |> cvtDiagOp
        U_mask_U   = U_mask    |> cvtDiagOp
        V_mask_V   = V_mask    |> cvtDiagOp
        W_mask_W   = W_mask    |> cvtDiagOp

        # ===== [ BEGIN grid length, area and volume ] =====
        
        Δx_T = gd.Δx_T
        Δx_U = gd.Δx_U
        Δx_V = gd.Δx_V
        Δx_W = gd.Δx_W
 
        Δy_T = gd.Δy_T
        Δy_U = gd.Δy_U
        Δy_V = gd.Δy_V
        Δy_W = gd.Δy_W
 
        Δz_T = gd.Δz_T
        Δz_U = gd.Δz_U
        Δz_V = gd.Δz_V
        Δz_W = gd.Δz_W
        
        Δax_T  = Δy_T  .* Δz_T
        Δax_U  = Δy_U  .* Δz_U
        Δax_V  = Δy_V  .* Δz_V
        Δax_W  = Δy_W  .* Δz_W
        
        Δay_T  = Δx_T  .* Δz_T
        Δay_U  = Δx_U  .* Δz_U
        Δay_V  = Δx_V  .* Δz_V
        Δay_W  = Δx_W  .* Δz_W

        Δaz_T  = Δx_T  .* Δy_T
        Δaz_U  = Δx_U  .* Δy_U
        Δaz_V  = Δx_V  .* Δy_V
        Δaz_W  = Δx_W  .* Δy_W

        Δv_T   = Δx_T  .* Δy_T  .* Δz_T
        Δv_U   = Δx_U  .* Δy_U  .* Δz_U
        Δv_V   = Δx_V  .* Δy_V  .* Δz_V
        Δv_W   = Δx_W  .* Δy_W  .* Δz_W
        
        T_Δax_T   = Δax_T  |> cvtDiagOp
        U_Δax_U   = Δax_U  |> cvtDiagOp
        V_Δax_V   = Δax_V  |> cvtDiagOp
        W_Δax_W   = Δax_W  |> cvtDiagOp

        T_Δay_T   = Δay_T  |> cvtDiagOp
        U_Δay_U   = Δay_U  |> cvtDiagOp
        V_Δay_V   = Δay_V  |> cvtDiagOp
        W_Δay_W   = Δay_W  |> cvtDiagOp

        T_Δaz_T   = Δaz_T  |> cvtDiagOp
        U_Δaz_U   = Δaz_U  |> cvtDiagOp
        V_Δaz_V   = Δaz_V  |> cvtDiagOp
        W_Δaz_W   = Δaz_W  |> cvtDiagOp

        T_Δv_T     = Δv_T  |> cvtDiagOp
        U_Δv_U     = Δv_U  |> cvtDiagOp
        V_Δv_V     = Δv_V  |> cvtDiagOp
        W_Δv_W     = Δv_W  |> cvtDiagOp

        T_invΔv_T  = (Δv_T.^(-1))  |> cvtDiagOp
        U_invΔv_U  = (Δv_U.^(-1))  |> cvtDiagOp
        V_invΔv_V  = (Δv_V.^(-1))  |> cvtDiagOp
        W_invΔv_W  = (Δv_W.^(-1))  |> cvtDiagOp
       
        T_Δx_T = ( Δx_T ) |> cvtDiagOp
        U_Δx_U = ( Δx_U ) |> cvtDiagOp
        V_Δx_V = ( Δx_V ) |> cvtDiagOp
        W_Δx_W = ( Δx_W ) |> cvtDiagOp

        T_Δy_T = ( Δy_T ) |> cvtDiagOp
        U_Δy_U = ( Δy_U ) |> cvtDiagOp
        V_Δy_V = ( Δy_V ) |> cvtDiagOp
        W_Δy_W = ( Δy_W ) |> cvtDiagOp

        T_Δz_T = ( Δz_T ) |> cvtDiagOp
        U_Δz_U = ( Δz_U ) |> cvtDiagOp
        V_Δz_V = ( Δz_V ) |> cvtDiagOp
        W_Δz_W = ( Δz_W ) |> cvtDiagOp


        T_invΔx_T = ( Δx_T.^(-1) ) |> cvtDiagOp
        T_invΔy_T = ( Δy_T.^(-1) ) |> cvtDiagOp
        T_invΔz_T = ( Δz_T.^(-1) ) |> cvtDiagOp
 
        U_invΔx_U = ( Δx_U.^(-1) ) |> cvtDiagOp
        U_invΔy_U = ( Δy_U.^(-1) ) |> cvtDiagOp
        U_invΔz_U = ( Δz_U.^(-1) ) |> cvtDiagOp
 
        V_invΔx_V = ( Δx_V.^(-1) ) |> cvtDiagOp
        V_invΔy_V = ( Δy_V.^(-1) ) |> cvtDiagOp
        V_invΔz_V = ( Δz_V.^(-1) ) |> cvtDiagOp

        W_invΔx_W = ( Δx_W.^(-1) ) |> cvtDiagOp
        W_invΔy_W = ( Δy_W.^(-1) ) |> cvtDiagOp
        W_invΔz_W = ( Δz_W.^(-1) ) |> cvtDiagOp
 
        T_Δv_T = T_Δx_T * T_Δy_T * T_Δz_T


        # ===== [ END grid length, area and volume ] =====

        # ===== [ BEGIN interpolation ] =====
        function selfDivision(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse
            for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg:_end] ./= wgt
                end
            end
            
            return transpose(m_t) |> sparse
        end
        
        ones_T  = ones(Float64, bmo.T_pts)
        ones_U  = ones(Float64, bmo.U_pts)
        ones_V  = ones(Float64, bmo.V_pts)
        ones_W  = ones(Float64, bmo.W_pts)
 
        U_interp_T = (bmo.U_W_T + bmo.U_E_T) * T_mask_T
        U_interp_T = selfDivision(U_interp_T, ones_T)
        
        V_interp_T = (bmo.V_S_T + bmo.V_N_T) * T_mask_T
        V_interp_T = selfDivision(V_interp_T, ones_T)
        
        W_interp_T = (bmo.W_DN_T + bmo.W_UP_T) * T_mask_T
        W_interp_T = selfDivision(W_interp_T, ones_T)
 
        T_interp_U = (bmo.T_E_U + bmo.T_W_U) * U_mask_U
        T_interp_U = selfDivision(T_interp_U, ones_U)

        T_interp_V = (bmo.T_N_V + bmo.T_S_V) * V_mask_V
        T_interp_V = selfDivision(T_interp_V, ones_V)
 
        T_interp_W = (bmo.T_DN_W + bmo.T_UP_W) * W_mask_W
        T_interp_W = selfDivision(T_interp_W, ones_W)
 
        V_interp_U = (bmo.V_N_T + bmo.V_S_T) * (bmo.T_E_U + bmo.T_W_U) * U_mask_U
        V_interp_U = selfDivision(V_interp_U, ones_U)
        
        U_interp_V = (bmo.U_E_T + bmo.U_W_T) * (bmo.T_N_V + bmo.T_S_V) * V_mask_V
        U_interp_V = selfDivision(U_interp_V, ones_V)

        W_interp_U = (bmo.W_UP_T + bmo.W_DN_T) * (bmo.T_E_U + bmo.T_W_U) * U_mask_U
        W_interp_U = selfDivision(W_interp_U, ones_U)

        W_interp_V = (bmo.W_UP_T + bmo.W_DN_T) * (bmo.T_N_V + bmo.T_S_V) * V_mask_V
        W_interp_V = selfDivision(W_interp_V, ones_V)

        U_interp_W = (bmo.U_E_T + bmo.U_W_T) * (bmo.T_UP_W + bmo.T_DN_W) * W_mask_W
        U_interp_W = selfDivision(U_interp_W, ones_W)

        V_interp_W = (bmo.V_N_T + bmo.V_S_T) * (bmo.T_UP_W + bmo.T_DN_W) * W_mask_W
        V_interp_W = selfDivision(V_interp_W, ones_W)

        # Sending EW then NS will be the same as sending NS then EW
        T_NEWSavg_T = T_interp_V * T_interp_U #(bmo.T_N_T + bmo.T_S_T + bmo.T_E_T + bmo.T_W_T) / 4 
        #T_NEWSavg_T = selfDivision(T_NEWSavg_T, ones_T)





        # ===== [ END interpolation ] =====

        # MAGIC!!
        T_DIVx_U  = T_mask_T   * T_invΔv_T   * ( bmo.T_W_U  - bmo.T_E_U  ) * U_Δax_U  ; dropzeros!(T_DIVx_U);
        T_DIVy_V  = T_mask_T   * T_invΔv_T   * ( bmo.T_S_V  - bmo.T_N_V  ) * V_Δay_V  ; dropzeros!(T_DIVy_V);
        T_DIVz_W  = T_mask_T   * T_invΔv_T   * ( bmo.T_DN_W - bmo.T_UP_W ) * W_Δaz_W  ; dropzeros!(T_DIVz_W);
        

        U_∂x_T = U_mask_U * U_invΔx_U * (bmo.U_W_T  - bmo.U_E_T)                 ; dropzeros!(U_∂x_T);
        V_∂y_T = V_mask_V * V_invΔy_V * (bmo.V_S_T  - bmo.V_N_T)                 ; dropzeros!(V_∂y_T);
        W_∂z_T = W_mask_W * W_invΔz_W * (bmo.W_DN_T - bmo.W_UP_T)                ; dropzeros!(W_∂z_T);

        T_∂x_U = T_mask_T * T_invΔx_T * ( bmo.T_W_U - bmo.T_E_U )               ; dropzeros!(T_∂x_U);
        T_∂y_V = T_mask_T * T_invΔy_T * ( bmo.T_S_V - bmo.T_N_V )               ; dropzeros!(T_∂y_V);
        T_∂z_W = T_mask_T * T_invΔz_T * ( bmo.T_DN_W - bmo.T_UP_W )             ; dropzeros!(T_∂z_W);
        
        U_∂x_U = U_mask_U * U_invΔx_U * (bmo.U_W_U  - bmo.U_E_U) / 2            ; dropzeros!(U_∂x_U);
        V_∂y_V = V_mask_V * V_invΔy_V * (bmo.V_S_V  - bmo.V_N_V) / 2            ; dropzeros!(V_∂y_V);
        W_∂z_W = W_mask_W * W_invΔz_W * (bmo.W_DN_W - bmo.W_UP_W) / 2          ; dropzeros!(W_∂z_W);
        
        T_LAPx_T   =  T_DIVx_U * U_∂x_T
        T_LAPy_T   =  T_DIVy_V * V_∂y_T
      
        return new(
            bmo,
            
            T_DIVx_U,
            T_DIVy_V,
            T_DIVz_W,
         
            U_∂x_T,
            V_∂y_T,
            W_∂z_T,

            T_∂x_U,
            T_∂y_V,
            T_∂z_W,

            U_∂x_U,
            V_∂y_V,
            W_∂z_W, 

            T_LAPx_T,
            T_LAPy_T,

            U_interp_T,
            V_interp_T,
            W_interp_T,

            T_interp_U,
            T_interp_V,
            T_interp_W,

            V_interp_U,
            U_interp_V,

            W_interp_U,
            W_interp_V,

            U_interp_W,
            V_interp_W,

            T_NEWSavg_T,

            T_mask_T,
            U_mask_U,
            V_mask_V,
            W_mask_W,

            T_Δx_T,
            U_Δx_U,
            V_Δx_V,
            W_Δx_W,

            T_Δy_T,
            U_Δy_U,
            V_Δy_V,
            W_Δy_W,

            T_Δz_T,
            U_Δz_U,
            V_Δz_V,
            W_Δz_W,

            T_Δv_T,
        )
    end
end
