mutable struct Env

    pp ::PhyParams
    gd     :: Grid
    gd_col :: Grid
    gd_slb :: Grid

    function Env(;
        Δx :: Float64,
        Δy :: Float64,
        Nx :: Int64,
        Ny :: Int64,
        Nz :: Int64,
        f0 :: Float64,
        ΔΘ :: Float64,
        Θ0     :: Float64,
        g0     :: Float64,
        h_0    :: Float64,
        A_h    :: Float64,
        γ_Θ    :: Float64,
        E0     :: Float64,
        dlnγdδ :: Float64,
        γ0     :: Float64,
    )
        
        gd = Grid(
            Δx = Δx,
            Δy = Δy,
            Nx = Nx,
            Ny = Ny,
            Nz = Nz,
        )

        gd_col = Grid(
            Δx = Δx,
            Δy = Δy,
            Nx = 1,
            Ny = 1,
            Nz = Nz,
        )

        gd_slb = Grid(
            Δx = Δx,
            Δy = Δy,
            Nx = Nx,
            Ny = Ny,
            Nz = 1,
        )


        s0 = gd_col.z_T[1, 1, 1]

        pp = PhyParams(
            f0,
            ΔΘ,
            Θ0,
            g0,
            h_0,
            A_h,
            γ_Θ,
            E0,
            dlnγdδ,
            γ0,
            s0,
        )



        return new(
            pp, gd, gd_col, gd_slb,
        )
        
    end
end
