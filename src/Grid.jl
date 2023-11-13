struct Grid

    Nx :: Int64
    Ny :: Int64
    Nz :: Int64
 
    mask :: AbstractArray{Int64, 3}

    Δx_T  :: AbstractArray{Float64, 3}
    Δy_T  :: AbstractArray{Float64, 3}
    Δz_T  :: AbstractArray{Float64, 3}

    Δx_U  :: AbstractArray{Float64, 3}
    Δy_U  :: AbstractArray{Float64, 3}
    Δz_U  :: AbstractArray{Float64, 3}
 
    Δx_V  :: AbstractArray{Float64, 3}
    Δy_V  :: AbstractArray{Float64, 3}
    Δz_V  :: AbstractArray{Float64, 3}
 
    Δx_W  :: AbstractArray{Float64, 3}
    Δy_W  :: AbstractArray{Float64, 3}
    Δz_W  :: AbstractArray{Float64, 3}

    x_T   :: AbstractArray{Float64, 3}
    x_U   :: AbstractArray{Float64, 3}
    x_V   :: AbstractArray{Float64, 3}
    x_W   :: AbstractArray{Float64, 3}

    y_T   :: AbstractArray{Float64, 3}
    y_U   :: AbstractArray{Float64, 3}
    y_V   :: AbstractArray{Float64, 3}
    y_W   :: AbstractArray{Float64, 3}

    z_T   :: AbstractArray{Float64, 3}
    z_U   :: AbstractArray{Float64, 3}
    z_V   :: AbstractArray{Float64, 3}
    z_W   :: AbstractArray{Float64, 3}

    function Grid(;
        Δx    :: Float64,
        Δy    :: Float64,
        Nx    :: Int64,
        Ny    :: Int64,
        Nz    :: Int64,
    )

        mask = ones(Int64, Nx, Ny, Nz)
        
        _z_T, _z_W, _Δz_T, _Δz_W = genVerticalGrid(Nz=Nz)
        _x_T, _x_U, _Δx_T, _Δx_U, _y_T, _y_V, _Δy_T, _Δy_V = genHorizontalGrid(Nx=Nx, Ny=Ny, Δx=Δx, Δy=Δy)

        z_makeMesh = (a, nx, ny) -> repeat( reshape(a, 1, 1, :), outer=(nx, ny,  1) )
        y_makeMesh = (a, nx, nz) -> repeat( reshape(a, 1, :, 1), outer=(nx, 1,  nz) )
        x_makeMesh = (a, ny, nz) -> repeat( reshape(a, :, 1, 1), outer=(1,  ny, nz) )

        z_T   = z_makeMesh(_z_T,  Nx,   Ny)
        z_U   = z_makeMesh(_z_T,  Nx,   Ny)
        z_V   = z_makeMesh(_z_T,  Nx,   Ny)
        z_W   = z_makeMesh(_z_W,  Nx,   Ny)
    
        Δz_T  = z_makeMesh(_Δz_T, Nx,   Ny)
        Δz_U  = z_makeMesh(_Δz_T, Nx,   Ny)
        Δz_V  = z_makeMesh(_Δz_T, Nx,   Ny)
        Δz_W  = z_makeMesh(_Δz_W, Nx,   Ny)

        y_T   = y_makeMesh(_y_T,  Nx,   Nz  )
        y_U   = y_makeMesh(_y_T,  Nx,   Nz  )
        y_V   = y_makeMesh(_y_V,  Nx,   Nz  )
        y_W   = y_makeMesh(_y_T,  Nx,   Nz+1)

        Δy_T   = y_makeMesh(_Δy_T, Nx,   Nz  )
        Δy_U   = y_makeMesh(_Δy_T, Nx,   Nz  )
        Δy_V   = y_makeMesh(_Δy_V, Nx,   Nz  )
        Δy_W   = y_makeMesh(_Δy_T, Nx,   Nz+1)

        x_T    = x_makeMesh(_x_T,  Ny,   Nz  )
        x_U    = x_makeMesh(_x_U,  Ny,   Nz  )
        x_V    = x_makeMesh(_x_T,  Ny,   Nz  )
        x_W    = x_makeMesh(_x_T,  Ny,   Nz+1)

        Δx_T   = x_makeMesh(_Δx_T, Ny,   Nz  )
        Δx_U   = x_makeMesh(_Δx_U, Ny,   Nz  )
        Δx_V   = x_makeMesh(_Δx_T, Ny,   Nz  )
        Δx_W   = x_makeMesh(_Δx_T, Ny,   Nz+1)

        return new(

            Nx,
            Ny,
            Nz,
         
            mask,

            Δx_T,
            Δy_T,
            Δz_T,

            Δx_U,
            Δy_U,
            Δz_U,

            Δx_V,
            Δy_V,
            Δz_V,
 
            Δx_W,
            Δy_W,
            Δz_W,
 
            x_T,
            x_U,
            x_V,
            x_W,

            y_T,
            y_U,
            y_V,
            y_W,

            z_T,
            z_U,
            z_V,
            z_W,

        ) 
        
    end
end


function genHorizontalGrid(;
    Nx :: Int64,
    Ny :: Int64,
    Δx :: Float64,
    Δy :: Float64,
)

    y_V = collect(Float64, 0:Ny) * Δy
    y_T = (y_V[1:end-1] + y_V[2:end]) / 2.0

    Δy_T = ones(Float64, Ny) * Δy
    Δy_V = ones(Float64, Ny+1) * Δy
 

    x_U = collect(Float64, 0:Nx) * Δx
    x_T = (x_U[1:end-1] + x_U[2:end]) / 2.0

    Δx_T = ones(Float64, Nx) * Δx
    Δx_U = ones(Float64, Nx+1) * Δx
    
    y_V = y_V[1:end-1]
    Δy_V = Δy_V[1:end-1]

    x_U = x_U[1:end-1]
    Δx_U = Δx_U[1:end-1]


    return x_T, x_U, Δx_T, Δx_U, y_T, y_V, Δy_T, Δy_V

end


function genVerticalGrid(;
    Nz    :: Int64,
)

    z_W = collect(Float64, range(0, 1, length=Nz+1))

    Δz_W = similar(z_W)

    Δz_T = z_W[2:end] - z_W[1:end-1]
    Δz_W[2:end-1] = ( Δz_T[2:end] + Δz_T[1:end-1] ) / 2.0
    Δz_W[1] = Δz_W[2]
    Δz_W[end] = Δz_W[end-1]

    z_T = (z_W[1:end-1] + z_W[2:end]) / 2.0

    return z_T, z_W, Δz_T, Δz_W
end
