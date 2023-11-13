mutable struct Env

    f0     :: Float64
    gd     :: Grid
    gd_col :: Grid
    gd_slb :: Grid

    function Env(;
        f0 :: Float64,
        Δx :: Float64,
        Δy :: Float64,
        Nx :: Int64,
        Ny :: Int64,
        Nz :: Int64,
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



        return new(
            f0, gd, gd_col, gd_slb,
        )
        
    end
end
