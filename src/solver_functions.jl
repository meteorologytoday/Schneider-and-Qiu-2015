
# E_0_W is the zero-th order E, while E0 is a constant
function getE_0(pp::PhyParams, z)
    tmp =  (z .+ pp.s0) / pp.γ0
#    return tmp .* E0 .* exp.( 1.0 .- tmp ) * 0 .+ E0
    return tmp .* pp.E0 .* exp.( 1.0 .- tmp )
end

function getE_0_ideal(pp::PhyParams, z)
    return getE_0(pp, z) * 0 .+ E0
end

function getE_1(pp::PhyParams, z, δ)
    return δ .* z ./ pp.γ0 .* pp.dlnγdδ .* getE_0(pp, z) 
end


function solve_order0!(
    m :: Model;
)
    st = m.st
    ev = m.ev
    pp = ev.pp
    gd_col = ev.gd_col

    Nz = gd_col.Nz

    amo_col = m.co.amo_col

    E_0_cW = getE_0(pp, gd_col.z_W[:])
    cW_E_0_cW = spdiagm(0=>E_0_cW)

    # This matrix rotate the 2-D vector 90 degress clockwise
    zhat_cross = spdiagm(Nz => - ones(Float64, Nz), - Nz => ones(Float64, Nz) )

    bottom_approx_op = spzeros(Nz+1, Nz)
    bottom_approx_op[1, 1] = - getE_0(pp, pp.s0) / pp.s0

    # This is just for u or v only
    tmp = - ( - amo_col.T_DIVz_W ) * ( - cW_E_0_cW * amo_col.W_mask_W * amo_col.W_∂z_T + bottom_approx_op )

    # Need to use blockdiag to apply for u and v together
    op = blockdiag(tmp, tmp) + pp.f0 * zhat_cross
    
    
    RHS = pp.f0 * zhat_cross * st.VEL_g
    
    println("Solving 0th order")
    @time F = lu(op)
    @time st.VEL_0[:] = F \ RHS


    #=
    dUds_0 = zeros(Float64, (Ns+1) * 2)
    dUds_0_x = view(dUds_0, 1:(Ns+1))
    dUds_0_y = view(dUds_0, (Ns+2):(2*(Ns+1)))

    # ===== constants =====

    reduced_g0 = g0 * ΔΘ / Θ0


    # ===== Solving for 0st order =====
    println("Solving for 0st order")
    # s0 = 0.05
    s0 = s_T[1]
    γ0 = 0.3


    =#

end

function solve_order1_thermal!(m :: Model; wn_rad :: Float64 = Inf)

    st = m.st
    ev = m.ev
    pp = ev.pp
    gd = ev.gd
    gd_col = ev.gd_col

    amo_slb = m.co.amo_slb
    bmo_slb = amo_slb.bmo

    Nx = gd.Nx
    Ny = gd.Ny
    Nz = gd.Nz

    Lx = gd.x_U[end, 1, 1] - gd.x_U[1, 1, 1]
    Ly = gd.y_V[1, end, 1] - gd.y_V[1, 1, 1]

    k_map, l_map = wavenumber_map2d(Nx, Ny)

    wnk_map = k_map * 2π/Lx
    wnl_map = l_map * 2π/Ly

    compute_map = (k_map.^2 + l_map.^2) .< wn_rad^2

    # Compute avg(u)
    st.u_0_mean = sum(st.u_0 .* gd_col.Δz_T[:]) / sum(gd_col.Δz_T[:])
    st.v_0_mean = sum(st.v_0 .* gd_col.Δz_T[:]) / sum(gd_col.Δz_T[:])

    VEL_0_mean = zeros(Float64, amo_slb.bmo.U_pts * 2)

    st.SST_1_coe[:, :] = AbstractFFTs.fft(reshape(st.SST_1, Nx, Ny))

    for i=1:Nx, j=1:Ny

        _k = k_map[i, j]
        _l = l_map[i, j]
     
        _wnk = wnk_map[i, j]
        _wnl = wnl_map[i, j]
        
        if compute_map[i, j] == false
            continue
        end
        
        op = sparse( [ ( 1 + ( im * (st.u_0_mean * _wnk + st.v_0_mean * _wnl) + pp.A_h * (_wnk^2 + _wnl^2) ) /  (pp.γ_Θ/pp.h_0) ) ])
        RHS = [ st.SST_1_coe[i, j]]
        
        sol = op \ RHS
        st.Θ_1_coe[i, j] = sol[1]
    end

    st.Θ_1[:, :] = real.(ifft(st.Θ_1_coe))
end


function solve_order1_momentum!(m :: Model; wn_rad :: Float64 = Inf)
 
    st = m.st
    ev = m.ev
    co = m.co
    pp = ev.pp
    gd = ev.gd
    gd_col = ev.gd_col

    Nx = gd.Nx
    Ny = gd.Ny
    Nz = gd.Nz
    
    reduced_g0 = pp.g0 * pp.ΔΘ / pp.Θ0

    Lx = gd.x_U[end, 1, 1] - gd.x_U[1, 1, 1]
    Ly = gd.y_V[1, end, 1] - gd.y_V[1, 1, 1]

    k_map, l_map = wavenumber_map2d(Nx, Ny)

    wnk_map = k_map * 2π/Lx
    wnl_map = l_map * 2π/Ly

    compute_map = (k_map.^2 + l_map.^2) .< wn_rad^2

    amo = co.amo_col
    bmo = amo.bmo

    s_T = ev.gd_col.z_T[:]
    s_W = ev.gd_col.z_W[:]

    T_pts = bmo.T_pts
    U_pts = bmo.U_pts
    V_pts = bmo.V_pts
    W_pts = bmo.W_pts

    mask_eff_W = reshape(amo.W_mask_W * ones(Float64, W_pts), amo.bmo.W_dim...)

    # Create coversion matrix and its inverse

    # w: the top and bottom w^*(1) are not free variables
    W_num = reshape(collect(1:length(mask_eff_W)), size(mask_eff_W)...)
    active_num_eff_W = W_num[ mask_eff_W .==1 ]
    eW_send_W = amo.bmo.W_I_W[active_num_eff_W, :]; dropzeros!(eW_send_W)
    W_send_eW = sparse(eW_send_W')

    full_send_eff = blockdiag(
        sparse(I, T_pts, T_pts),
        sparse(I, T_pts, T_pts),
        W_send_eW,
        sparse(I, 1, 1),
    )

    u_0_T = copy(st.u_0)[:]
    v_0_T = copy(st.v_0)[:]
    T_u_0_T = spdiagm(0 => u_0_T )
    T_v_0_T = spdiagm(0 => v_0_T )

    # building vertical advection operator
    duds_0_W = amo.W_∂z_T * u_0_T
    duds_0_W[1] = u_0_T[1] / s_T[1]
    duds_0_W[end] = 0.0
    duds_0_T = amo.T_interp_W * duds_0_W
    T_duds_0_T = spdiagm(0=>duds_0_T)

    dvds_0_W = amo.W_∂z_T * v_0_T
    dvds_0_W[1] = v_0_T[1] / s_T[1]
    dvds_0_W[end] = 0.0
    dvds_0_T = amo.T_interp_W * duds_0_W
    T_dvds_0_T = spdiagm(0=>dvds_0_T)

    E_0_W = getE_0(pp, gd_col.z_W[:])
    W_E_0_W = spdiagm(0=>E_0_W)

    δ_1_coe = st.SST_1_coe - st.Θ_1_coe

    # Nx and Ny are in spectral sense.
    # Notice the arrangement: the fastest running index
    # runs within the same choice of wavenumber pair (k, l)
    coe_vec = zeros(ComplexF64, Nz * 3 + 2, Nx, Ny)
        
    for i=1:Nx, j=1:Ny

        _k = k_map[i, j]
        _l = l_map[i, j]
     
        _wnk = wnk_map[i, j]
        _wnl = wnl_map[i, j]
        
        if compute_map[i, j] == false
            continue
        end

        zero_wn = ( _k == _l == 0 )
        
        u_0_κ_T = u_0_T * _wnk + v_0_T * _wnl

        op_cont = [
            (sparse(I, Nz, Nz) * (im*pp.h_0*_wnk))  (sparse(I, Nz, Nz) * (im*pp.h_0*_wnl))  (pp.h_0 * amo.T_DIVz_W)  reshape( im*u_0_κ_T, :, 1)
        ]
        
        #println("Size of op_cont: ",    size(op_cont))

        # ========= Creating momentum operator =========
        # Coriolis operator, equivalent to `f0 ẑ ×`
        op_mom_cori = pp.f0 * sparse([
            spzeros(T_pts, T_pts)  (- bmo.T_I_T)           spzeros(T_pts, W_pts + 1) ;
            bmo.T_I_T              spzeros(T_pts, T_pts)   spzeros(T_pts, W_pts + 1) ; 
        ])

        op_mom_hadv= [
            spdiagm(0=>im*u_0_κ_T)  spzeros(T_pts, T_pts)   spzeros(T_pts, W_pts)  spzeros(T_pts, 1)  ;
            spzeros(T_pts, T_pts)   spdiagm(0=>im*u_0_κ_T)  spzeros(T_pts, W_pts)  spzeros(T_pts, 1)  ;
        ]

        op_mom_vadv= [
            spzeros(T_pts, 2*T_pts)  T_duds_0_T*amo.T_interp_W  spzeros(T_pts, 1)  ;
            spzeros(T_pts, 2*T_pts)  T_dvds_0_T*amo.T_interp_W  spzeros(T_pts, 1)  ;
        ]

        op_mom_h = [
            spzeros(T_pts, 2*T_pts + W_pts)  ((reduced_g0 * im * _wnk)*ones(T_pts, 1)) ;
            spzeros(T_pts, 2*T_pts + W_pts)  ((reduced_g0 * im * _wnl)*ones(T_pts, 1)) ;
        ]

        rmbot_W = ones(Float64, bmo.W_dim...)
        rmtop_W = ones(Float64, bmo.W_dim...)

        rmbot_W[:, :, 1]   .= 0.0
        rmtop_W[:, :, end] .= 0.0

        W_rmbot_W = spdiagm(0=>view(rmbot_W, :))
        W_rmtop_W = spdiagm(0=>view(rmtop_W, :))
        W_mask_W = W_rmbot_W * W_rmtop_W

        W_bottom_approx_op_x_T = spzeros(W_pts, T_pts)
        W_bottom_approx_op_y_T = spzeros(W_pts, T_pts)
       
         
        E_0_s0 = getE_0(pp, s_T[1])
        W_bottom_approx_op_x_T[1, 1] = - E_0_s0 / s_T[1]
        W_bottom_approx_op_y_T[1, 1] = - E_0_s0 / s_T[1]

        op_mom_vdif = - [
            (- amo.T_DIVz_W) * ( W_mask_W * (- W_E_0_W) * amo.W_∂z_T + W_bottom_approx_op_x_T )   spzeros(T_pts, T_pts)                                     spzeros(T_pts, W_pts + 1) ;
            spzeros(T_pts, T_pts)                                     (- amo.T_DIVz_W) * ( W_mask_W * (- W_E_0_W) * amo.W_∂z_T + W_bottom_approx_op_y_T )   spzeros(T_pts, W_pts + 1) ;
        ]

        op_mom = op_mom_hadv + op_mom_vadv + op_mom_h + op_mom_vdif + op_mom_cori
       
        # Construct right-hand-side

        hgt_factor = 1.0 .- s_T
        RHS_mom_pt = [
            (reduced_g0 * pp.h_0 / pp.ΔΘ * im * _wnk) * hgt_factor .* ones(T_pts, 1) * st.Θ_1_coe[i, j] ;
            (reduced_g0 * pp.h_0 / pp.ΔΘ * im * _wnl) * hgt_factor .* ones(T_pts, 1) * st.Θ_1_coe[i, j] ;
        ]

        W_E_1_W = spdiagm(0 => getE_1(pp, s_W, δ_1_coe[i, j]))

        RHS_mom_vdif = ( 
            (- blockdiag(amo.T_DIVz_W, amo.T_DIVz_W))
          * blockdiag(- W_E_1_W, - W_E_1_W)
          * [ duds_0_W ; dvds_0_W ; ]
        )

        if zero_wn == true
            op_full = op_mom[:, 1:(2*Nz)]
            RHS = (RHS_mom_pt + RHS_mom_vdif)[1:(2*Nz)]
            sol = solveComplexMatrix(op_full, RHS)
            coe_vec[1:(2*Nz), i, j] = sol
        else 
            op_full = [ op_mom ; op_cont ; ] * full_send_eff
            RHS = [ (RHS_mom_pt + RHS_mom_vdif) ; zeros(ComplexF64, T_pts) ; ]
            sol = full_send_eff * solveComplexMatrix(op_full, RHS)
            coe_vec[:, i, j] = sol
        end
    end

    # Inverse for each layer
    idx=0;
    u_1_coe = view(coe_vec, (idx+1):(idx+Nz), :, :);    idx += Nz;
    v_1_coe = view(coe_vec, (idx+1):(idx+Nz), :, :);    idx += Nz;
    w_1_coe = view(coe_vec, (idx+1):(idx+Nz+1), :, :);  idx += Nz+1;
    h_1_coe = view(coe_vec, (idx+1):(idx+1), :, :);     idx += 1;


    for k = 1:Nz
        st.u_1[:, :, k] = real.(ifft(u_1_coe[k, :, :])) 
        st.v_1[:, :, k] = real.(ifft(v_1_coe[k, :, :])) 
    end

    for k = 1:Nz+1
        st.w_1[:, :, k] = real.(ifft(w_1_coe[k, :, :]))
    end

    st.h_1[:, :] = real.(ifft(h_1_coe[1, :, :]))


end
