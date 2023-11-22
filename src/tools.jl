using FFTW
using AbstractFFTs

function wavenumber_map(N)

    if N < 3
        throw(ErrorException("N must be greater than 3"))
    end

    if mod(N, 2) != 0
        throw(ErrorException("N must be an even number"))
    end


    k_map = zeros(Float64, N)

    half_N = Int(N/2)

    k_map[1]        = 0.0
    k_map[half_N+1] = half_N

    for i=2:half_N
        k_map[i] = i-1
        k_map[N-i+2] = - (i-1)
    end
    
    return k_map

end

function wavenumber_map2d(Nx, Ny)

    k_map = wavenumber_map(Nx)
    l_map = wavenumber_map(Ny)

    k_map = repeat(reshape(k_map, :, 1), outer=(1, Ny))
    l_map = repeat(reshape(l_map, 1, :), outer=(Nx, 1))

    return k_map, l_map
end



function solveComplexMatrix(A, b)

    Ar = real.(A)
    Ac = imag.(A)

    br = real.(b)
    bc = imag.(b)

    M = [ Ar  (-Ac) ;
          Ac    Ar  ;  ]

    B = [ br ; bc ; ]

    F = lu(M)
    X = F \ B

    N = length(b)
    xr = view(X, 1:N)
    xi = view(X, (N+1):(2*N))

    return xr + im * xi

end


