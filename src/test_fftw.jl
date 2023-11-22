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



Lx = 30.0
Ly = 25.0
k = 8.0
l = 2.0
Nx = 100
Ny = 100

f(x, y) = cos.(k * x * π / Lx) .* sin.(l * y * π / Ly)

dfdx(x, y) = - k*π/Lx * sin.(k * x * π / Lx) .* sin.(l * y * π / Ly)
dfdy(x, y) =   l*π/Ly * cos.(k * x * π / Lx) .* cos.(l * y * π / Ly)


x_U_vec = collect(range(0.0, Lx, length=Nx+1))
y_V_vec = collect(range(0.0, Ly, length=Ny+1))

x_T_vec = (x_U_vec[1:end-1] + x_U_vec[2:end] ) / 2
y_T_vec = (y_V_vec[1:end-1] + y_V_vec[2:end] ) / 2


x_T = repeat(reshape(x_T_vec, :, 1), outer=(1, Ny))
y_T = repeat(reshape(y_T_vec, 1, :), outer=(Nx, 1))


F = f(x_T, y_T)

dFdx = dfdx(x_T, y_T)
dFdy = dfdy(x_T, y_T)
println("size(F) = ", size(F))

F_coe = AbstractFFTs.fft(F)

dFdx_coe = copy(F_coe)
dFdy_coe = copy(F_coe)


k_map, l_map = wavenumber_map2d(Nx, Ny)
for i = 1:Nx, j = 1:Ny 
    _k = k_map[i, j]
    _l = l_map[i, j]

    dFdx_coe[i, j] = F_coe[i, j] * ( 2*π/Lx * _k * im)
    dFdy_coe[i, j] = F_coe[i, j] * ( 2*π/Ly * _l * im)
end

dFdx_num = real.(ifft(dFdx_coe))
dFdy_num = real.(ifft(dFdy_coe))

println("dFdx[:, 1] = ",dFdx[:, 1])
println("dFdx_num[:, 1] = ", dFdx_num[:, 1])

println("dFdy[:, 1] = ",dFdy[:, 1])
println("dFdy_num[:, 1] = ", dFdy_num[:, 1])


diff = ( (dFdx - dFdx_num).^2 + (dFdy - dFdy_num).^2 ).^0.5

println("Max gradient difference : ", maximum(diff))


println("Loading PyPlot...")
using PyPlot
plt = PyPlot
println("Done.")

fig, ax = plt.subplots(1, 3)

cs = ax[1].contour(x_T_vec, y_T_vec, F', 10, colors="black")
plt.clabel(cs)

ax[1].streamplot(x_T_vec, y_T_vec, dFdx', dFdy')


cs = ax[2].contour(x_T_vec, y_T_vec, F', 10, colors="black")
plt.clabel(cs)

ax[2].streamplot(x_T_vec, y_T_vec, dFdx_num', dFdy_num')


cs = ax[3].contour(x_T_vec, y_T_vec, diff', 10, colors="black")
plt.clabel(cs)




plt.show()










