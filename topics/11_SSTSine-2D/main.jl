using Distributed

@everywhere include("../../src/BLM.jl")
@everywhere using Formatting





ΔΘs = [0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
Ugs = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
Vg = 0.0

wvlens = [20, 40, 60, 80, 100, 100, 120, 140, 160, 180, 200] * 1.0  # km
ΔSSTs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, ]



@everywhere output_dir = format("output_strong_E0")
mkpath(output_dir)

#@distributed for Ug in Ugs, ΔSST in ΔSSTs
@sync @distributed for (ΔΘ, Ug, wvlen, ΔSST) in collect(Iterators.product(ΔΘs, Ugs, wvlens, ΔSSTs))
    
    output_filename = format("DTheta_{:.1f}-Ug_{:.1f}-wvlen_{:03d}-dSST_{:.2f}.nc", ΔΘ, Ug, wvlen, ΔSST)
    full_filename = joinpath(output_dir, output_filename) 
    
    printfmtln("Doing the case: (ΔΘ, Ug, wvlen, ΔSST) = ({:f}, {:f}, {:d}, {:f}). Output: {:s}", ΔΘ, Ug, wvlen, ΔSST, output_filename)
    
    if isfile(full_filename)
        
        println("File $full_filename already exist. Skip solver.")
       
    else

        println("File $full_filename does not exist. Solve the model now...")
        println("Instantiating a boundary layer model...")

        Δx = 1.0 # km
        Δy = Δx  # km
        Nx = Int64(wvlen / Δx)
        Ny = 4

        m = BLM.Model(BLM.Env(
            Δx = Δx * 1e3,
            Δy = Δy * 1e3,
            Nx = Nx,
            Ny = Ny,
            Nz = 10,
            f0 = 1e-4,
            ΔΘ = ΔΘ,
            Θ0 = 290.0,
            g0 = 9.81,
            h_0 = 1.5e3,
            A_h = 1e5,
            γ_Θ = 0.25,
            E0  = 5.7e-6,
            dlnγdδ = 0.5,
            γ0     = 0.3,
        ))
        println("Done.")

        st = m.st
        ev = m.ev

        gd_slb = ev.gd_slb
        s_T = ev.gd_col.z_T[:]

        Lx = gd_slb.x_U[end, 1, 1] - gd_slb.x_U[1, 1, 1]
        Ly = gd_slb.y_V[1, end, 1] - gd_slb.y_V[1, 1, 1]

        wn_rad = gd_slb.Nx / 2
        
        println("# Physical parameters:")
        BLM.printInfo(m.ev.pp)

        # Setup geostrophic wind (background pressure gradient)
        st.u_g .= Ug
        st.v_g .= Vg


        # Setup SST_1
        x_c = Lx / 2.0
        y_c = Ly / 2.0
        #st.SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 + (gd_slb.y_T .- y_c).^2 ) / 50e3^2 )
        #st.SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 ) / 50e3^2 )
        
        st.SST_1[:] = ΔSST * sin.(2 * π * gd_slb.x_T / Lx)



        println("Now solve the 0-th order")
        BLM.solve_order0!(m)
        println("Done solving the 0-th order")

        println("Now solve the 1-th order thermal")
        BLM.solve_order1_thermal!(m; wn_rad=wn_rad)
        println("Done solving the 1-th order thermal")

        println("Now solve the 1-th order momentum")
        BLM.solve_order1_momentum!(m; wn_rad=wn_rad)
        println("Done solving the 1-th order momentum")

        BLM.savedata(m, full_filename, overwrite=true, props=Dict(
            "Ug" => Ug,
            "wvlen" => wvlen,
            "dSST" => ΔSST,
        ))

        m = Nothing

    end

end







