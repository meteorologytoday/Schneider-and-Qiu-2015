mutable struct State

    VEL_g :: AbstractArray{Float64, 1}
    u_g   :: AbstractArray{Float64, 1}
    v_g   :: AbstractArray{Float64, 1}

    VEL_0 :: AbstractArray{Float64, 1}
    u_0   :: AbstractArray{Float64, 1}
    v_0   :: AbstractArray{Float64, 1}

    dVELds_0 :: AbstractArray{Float64, 1}
    duds_0 :: AbstractArray{Float64, 1}
    dvds_0 :: AbstractArray{Float64, 1}

    u_0_mean   :: Float64
    v_0_mean   :: Float64



    SST_1 :: AbstractArray{Float64, 2}
    Θ_1 :: AbstractArray{Float64, 2}
    
    SST_1_coe :: AbstractArray{ComplexF64, 2}
    Θ_1_coe   :: AbstractArray{ComplexF64, 2}

    u_1 :: AbstractArray{Float64, 3}
    v_1 :: AbstractArray{Float64, 3}
    w_1 :: AbstractArray{Float64, 3}
    h_1 :: AbstractArray{Float64, 2}

    function State(
        ev :: Env,
    )

        Nx = ev.gd.Nx
        Ny = ev.gd.Ny
        Nz = ev.gd.Nz
        
        VEL_g= zeros(Float64, Nz * 2)
        u_g = view(VEL_g, 1:Nz)
        v_g = view(VEL_g, (Nz+1):(2*Nz))

        VEL_0 = zeros(Float64, Nz * 2)
        u_0 = view(VEL_0, 1:Nz)
        v_0 = view(VEL_0, (Nz+1):(2*Nz))

        dVELds_0 = zeros(Float64, (Nz+1) * 2)
        duds_0 = view(dVELds_0, 1:(Nz+1))
        dvds_0 = view(dVELds_0, (Nz+2):(2*(Nz+1)))
 
        SST_1 = zeros( Float64, Nx, Ny) 
        Θ_1 = zeros( Float64, Nx, Ny) 
        
        SST_1_coe = zeros(ComplexF64, Nx, Ny) 
        Θ_1_coe = zeros(ComplexF64, Nx, Ny) 
 
        u_1   = zeros(Float64, Nx, Ny, Nz)
        v_1   = zeros(Float64, Nx, Ny, Nz)
        w_1   = zeros(Float64, Nx, Ny, Nz+1)
        h_1   = zeros(Float64, Nx, Ny)
           
        return new(

            VEL_g,
            u_g,
            v_g,

            VEL_0,
            u_0,
            v_0,

            dVELds_0,
            duds_0,
            dvds_0,

            0.0,
            0.0,

            SST_1,
            Θ_1,

            SST_1_coe,
            Θ_1_coe,
        
            u_1,    
            v_1,    
            w_1,    
            h_1,    
        )    
    end
end

