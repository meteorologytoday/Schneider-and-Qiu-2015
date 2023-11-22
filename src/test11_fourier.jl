using SparseArrays
using LinearAlgebra, Arpack
using YAXArrays, DimensionalData, NetCDF
include("BLM.jl")
include("tools.jl")

plot_temperature = false

m = BLM.Model(BLM.Env(
    Δx = 5e3,
    Δy = 5e3,
    Nx = 200,
    Ny = 200,
    Nz = 10,
    f0 = 1e-4,
))


gd = m.ev.gd
Lx = gd.x_U[end, 1, 1] - gd.x_U[1, 1, 1]
Ly = gd.y_V[1, end, 1] - gd.y_V[1, 1, 1]
gd_col = m.ev.gd_col
gd_slb = m.ev.gd_slb
Ns = gd_col.Nz
s_T = gd_col.z_T[:]
s_W = gd_col.z_W[:]

wn_rad = 100 

k_map, l_map = wavenumber_map2d(gd.Nx, gd.Ny)

wnk_map = k_map * 2π/Lx
wnl_map = l_map * 2π/Ly

compute_map = (k_map.^2 + l_map.^2) .< wn_rad^2
display(compute_map)



VEL_0 = zeros(Float64, Ns * 2)
u_0 = view(VEL_0, 1:Ns)
v_0 = view(VEL_0, (Ns+1):(2*Ns))

VEL_g = copy(VEL_0)
ug = view(VEL_g, 1:Ns)
vg = view(VEL_g, (Ns+1):(2*Ns))


dUds_0 = zeros(Float64, (Ns+1) * 2)
dUds_0_x = view(dUds_0, 1:(Ns+1))
dUds_0_y = view(dUds_0, (Ns+2):(2*(Ns+1)))

# ===== constants =====

ΔΘ = 10.0
Θ0 = 290.0
g0 = 10.0
reduced_g0 = g0 * ΔΘ / Θ0


# ===== Solving for 0st order =====
println("Solving for 0st order")
# Construct operator
f0 = m.ev.f0
E0 = 2.5e-5
dlnγdδ = 0.5
Ug = 11.0 

ug .= Ug
vg .= 0.0

# s0 = 0.05
s0 = s_T[1]
γ0 = 0.3


# E_0_W is the zero-th order E, while E0 is a constant
function getE_0(z)
    tmp =  (z .+ s0) / γ0
#    return tmp .* E0 .* exp.( 1.0 .- tmp ) * 0 .+ E0
    return tmp .* E0 .* exp.( 1.0 .- tmp )
end

function getE_0_ideal(z)
    return getE_0(z) * 0 .+ E0
end

function getE_1(z, δ)
    return δ .* z ./ γ0 .* dlnγdδ .* getE_0(z) 
end

amo_col = m.co.amo_col

E_0_cW = getE_0(gd_col.z_W[:])
cW_E_0_cW = spdiagm(0=>E_0_cW)


# This matrix rotate the 2-D vector 90 degress clockwise
zhat_cross = spdiagm(Ns => - ones(Float64, Ns), - Ns => ones(Float64, Ns) )

keepsfc_W = zeros(Float64, Ns+1)
keepsfc_W[1] = 1.0
W_keepsfc_W = spdiagm(0=>keepsfc_W)
bottom_approx_op = (- getE_0(s0) / s0) * W_keepsfc_W * amo_col.bmo.W_DN_W * amo_col.W_interp_T

tmp = - ( - amo_col.T_DIVz_W ) * ( - cW_E_0_cW * amo_col.W_mask_W * amo_col.W_∂z_T + bottom_approx_op )
op_0th = blockdiag(tmp, tmp) + f0 * zhat_cross
RHS_0th = f0 * zhat_cross * VEL_g


println("Solving 0th order")
@time F = lu(op_0th)
@time VEL_0[:] = F \ RHS_0th



# Compute duds_0 and such







println("Compute idealized solution (assuming constant E = E0). E0 = ", E0)

k = sqrt( im * f0 / E0 )
A = [ [ exp(2*k)              (-1)                      ]
      [ (k - exp(k*s0)/s0)    (- k - exp(-k*s0)/s0)   ] ]

RHS = [0.0 ; ( Ug / s0) ;]
coefficients = A \ RHS

println("Solved coefficients: ", coefficients)
vec_u_0_ana = Ug .+ ( coefficients[1] * exp.(k*s_T) + coefficients[2] * exp.(-k*s_T) )
u_0_ana = real.(vec_u_0_ana)
v_0_ana = imag.(vec_u_0_ana)

#=
println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

println("Plot solved 0-th order solution")
fig, ax = plt.subplots(1, 3)
ax[1].plot(u_0, s_T, "k-", label="u_0")
ax[2].plot(v_0, s_T, "k-", label="v_0")

ax[1].plot(u_0_ana, s_T, "r-", label="u_0_ana")
ax[2].plot(v_0_ana, s_T, "r-", label="v_0_ana")


ax[1].plot(ug, s_T, "k--", label="u_g")
ax[2].plot(vg, s_T, "k--", label="v_g")


ax[1].set_title("zonal wind")
ax[2].set_title("meridional wind")

ax[1].legend()
ax[2].legend()

ax[3].plot(getE_0(s_W), s_W, "k-", label="E_0")
ax[3].plot(getE_1(s_W, 0.5), s_W, "r--", label="E_1")
ax[3].legend()

#fig.savefig("fig_zero.png", dpi=200)

plt.show()
=#

#@time VEL_0[:] = F \ (f0 * ROTCW * VEL_g - blockdiag(amo_col.T_DIVz_W, amo_col.T_DIVz_W) * blockdiag(W_E_0_W, W_E_0_W) * dUds_0)


# ===== Solving for 1st order =====
println("Solving for 1st order")
amo_slb = m.co.amo_slb
bmo_slb = amo_slb.bmo
h_0 = 1e3 # m
A_h = 1e5 # m^2 / s

# γ_Θ / (h_0 * f0) = 0.6
#γ_Θ = 0.25 * h_0
γ_Θ = h_0 * abs(f0) * 0.6
ΔSST  = 1.0 # K
# Compute avg(u)
u_0_mean = sum(u_0 .* gd_col.Δz_T[:])
v_0_mean = sum(v_0 .* gd_col.Δz_T[:])

VEL_0_mean = zeros(Float64, amo_slb.bmo.U_pts * 2)

SST_1 = zeros(Float64, amo_slb.bmo.T_pts)


x_c = Lx / 2.0
y_c = Ly / 2.0

SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 + (gd_slb.y_T .- y_c).^2 ) / 50e3^2 )
#SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 ) / 50e3^2 )

SST_1_coe = AbstractFFTs.fft(reshape(SST_1, gd.Nx, gd.Ny))
Θ_1_coe = copy(SST_1_coe)


println("Solve for 1st order temperature")
for i=1:gd.Nx, j=1:gd.Ny

    _k = k_map[i, j]
    _l = l_map[i, j]
 
    _wnk = wnk_map[i, j]
    _wnl = wnl_map[i, j]
    
    #println("Solve for k, l = (", _k, ", ", _l, ")")
   
    if compute_map[i, j] == false
        #println("Not in the computed radius. Skip.")
        continue
    end
    

    op_thermal = sparse( [ ( 1 + ( im * (u_0_mean * _wnk + v_0_mean * _wnl) + A_h * (_wnk^2 + _wnl^2) ) /  (γ_Θ/h_0) ) ])
    RHS = [ SST_1_coe[i, j]]
    
    sol = op_thermal \ RHS #solveComplexMatrix(op_thermal, RHS)    
    Θ_1_coe[i, j] = sol[1]
end
println("Done.")

δ_1_coe = SST_1_coe - Θ_1_coe




# Inverse into physical space!
slb_size = (gd_slb.Nx, gd_slb.Ny)
Θ_1 = real.(ifft(Θ_1_coe))
SST_1 = reshape(SST_1, slb_size...)

if plot_temperature 
    println("Loading PyPlot")
    using PyPlot
    plt = PyPlot
    println("Done")
    println("Plot solved 1-th order thermal solution (2D)")

    fig, ax = plt.subplots(1, 2)


    levels = collect(range(-2, 2, length=41))

    x_T = gd_slb.x_T[:, 1, 1] / 1e3
    y_T = gd_slb.y_T[1, :, 1] / 1e3
    cs = ax[1].contour(x_T, y_T, SST_1', levels, colors="black")
    plt.clabel(cs)

    cs = ax[2].contour(x_T, y_T, Θ_1', levels, colors="black")
    plt.clabel(cs)

    for _ax in ax
        _ax.set_xlabel("x [km]")
        _ax.set_xlabel("y [km]")
        _ax.grid()
    end

    ax[1].set_title("\$  \\mathrm{SST}^{\\left(1\\right)} \$")
    ax[2].set_title("\$  \\Theta^{\\left(1\\right)} \$")

    plt.show()
end
#δ_1 = SST_1 - Θ_1 # sT grid


# convert δ_1 to δ_1_W

#=
δ_1_W = repeat(δ_1, gd_col.Nz+1)
E_1_W = δ_1_W * dlnEdδ * E_0_W
=#

# 2: constrcut effective grid and sender matrices
# find effective grid of vertical velocities
# Vertical velocity: top and bottom are zeros

# From now on it is only 1 dimensional

amo = amo_col
bmo = amo.bmo

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

println("Size of full_send_eff: ", size(full_send_eff))

# 3: formulate operator
u_0_T = copy(u_0)[:]
v_0_T = copy(v_0)[:]
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

E_0_W = getE_0(gd_col.z_W[:])
W_E_0_W = spdiagm(0=>E_0_W)



println("T_pts  : ", T_pts)
println("U_pts  : ", U_pts)
println("V_pts  : ", V_pts)
println("W_pts  : ", W_pts)

# Nx and Ny are in spectral sense.
# Notice the arrangement: the fastest running index
# runs within the same choice of wavenumber pair (k, l)
coe_vec = zeros(ComplexF64, Ns * 3 + 2, gd.Nx, gd.Ny)
    
println("Solving 1st order momentum and height")
for i=1:gd.Nx, j=1:gd.Ny

    _k = k_map[i, j]
    _l = l_map[i, j]
 
    _wnk = wnk_map[i, j]
    _wnl = wnl_map[i, j]
    
    #println("Solve for k, l = (", _k, ", ", _l, ")")
   
    if compute_map[i, j] == false
        #println("Not in the computed radius. Skip.")
        continue
    end

    zero_wn = ( _k == _l == 0 )
    #u_1_coe = view(coe_vec, (i-1)Ns)
    #v_1_coe = zeros(ComplexF64, Ns)
    #w_1_coe = zeros(ComplexF64, Ns)
    
    u_0_κ_T = u_0_T * _wnk + v_0_T * _wnl

    op_cont = [
        (sparse(I, Ns, Ns) * (im*h_0*_wnk))  (sparse(I, Ns, Ns) * (im*h_0*_wnl))  (h_0 * amo.T_DIVz_W)  reshape( im*u_0_κ_T, :, 1)
    ]
    
    #println("Size of op_cont: ",    size(op_cont))

    # ========= Creating momentum operator =========
    # Coriolis operator, equivalent to `f0 ẑ ×`
    op_mom_cori = m.ev.f0 * sparse([
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
   
     
    E_0_s0 = getE_0(s_T[1])
    W_bottom_approx_op_x_T[1, 1] = - E_0_s0 / s_T[1]
    W_bottom_approx_op_y_T[1, 1] = - E_0_s0 / s_T[1]

    op_mom_vdif = - [
        (- amo.T_DIVz_W) * ( W_mask_W * (- W_E_0_W) * amo.W_∂z_T + W_bottom_approx_op_x_T )   spzeros(T_pts, T_pts)                                     spzeros(T_pts, W_pts + 1) ;
        spzeros(T_pts, T_pts)                                     (- amo.T_DIVz_W) * ( W_mask_W * (- W_E_0_W) * amo.W_∂z_T + W_bottom_approx_op_y_T )   spzeros(T_pts, W_pts + 1) ;
    ]

    #println("Size of op_mom_hadv: ", size(op_mom_hadv))
    #println("Size of op_mom_vadv: ", size(op_mom_vadv))
    #println("Size of op_mom_cori: ", size(op_mom_cori))
    #println("Size of op_mom_vdif: ", size(op_mom_vdif))
    #println("Size of op_mom_h: ",    size(op_mom_h))

    op_mom = op_mom_hadv + op_mom_vadv + op_mom_h + op_mom_vdif + op_mom_cori
   
    # Construct right-hand-side

    hgt_factor = 1.0 .- s_T
    RHS_mom_pt = [
        (reduced_g0 * h_0 / ΔΘ * im * _wnk) * hgt_factor .* ones(T_pts, 1) * Θ_1_coe[i, j] ;
        (reduced_g0 * h_0 / ΔΘ * im * _wnl) * hgt_factor .* ones(T_pts, 1) * Θ_1_coe[i, j] ;
    ]

    W_E_1_W = spdiagm(0 => getE_1(s_W, δ_1_coe[i, j]))

    RHS_mom_vdif = ( 
        (- blockdiag(amo.T_DIVz_W, amo.T_DIVz_W))
      * blockdiag(- W_E_1_W, - W_E_1_W)
      * [ duds_0_W ; dvds_0_W ; ]
    )

    #=
    τx_1_s0_W[1, 1, 1] = - E_0_s0 * u_0_T[1] / s_T[1]
    τy_1_s0_W[1, 1, 1] = - E_0_s0 * v_0_T[1] / s_T[1]

    RHS_mom_vdif_bc = (
        blockdiag(- amo.T_DIVz_W, - amo.T_DIVz_W) * [ τx_1_E_1_part_W[:] ; τy_1_E_1_part_W[:] ; ]
    )
    =#


    if zero_wn == true
        op_full = op_mom[:, 1:(2*Ns)]
        RHS = (RHS_mom_pt + RHS_mom_vdif)[1:(2*Ns)]
        sol = solveComplexMatrix(op_full, RHS)
        coe_vec[1:(2*Ns), i, j] = sol
    else 
        op_full = [ op_mom ; op_cont ; ] * full_send_eff
        RHS = [ (RHS_mom_pt + RHS_mom_vdif) ; zeros(ComplexF64, T_pts) ; ]
        sol = full_send_eff * solveComplexMatrix(op_full, RHS)
        coe_vec[:, i, j] = sol
    end
end
println("Done solving 1-st order.")

println("Now ready to inverse.")

# Inverse for each layer
idx=0;
u_1_coe = view(coe_vec, (idx+1):(idx+Ns), :, :);    idx += Ns;
v_1_coe = view(coe_vec, (idx+1):(idx+Ns), :, :);    idx += Ns;
w_1_coe = view(coe_vec, (idx+1):(idx+Ns+1), :, :);  idx += Ns+1;
h_1_coe = view(coe_vec, (idx+1):(idx+1), :, :);     idx += 1;

u_1   = zeros(Float64, gd.Nx, gd.Ny, Ns)
v_1   = zeros(Float64, gd.Nx, gd.Ny, Ns)
w_1   = zeros(Float64, gd.Nx, gd.Ny, Ns+1)
h_1   = zeros(Float64, gd.Nx, gd.Ny)

for k = 1:Ns
    u_1[:, :, k] = real.(ifft(u_1_coe[k, :, :])) 
    v_1[:, :, k] = real.(ifft(v_1_coe[k, :, :])) 
end

for k = 1:Ns+1
    w_1[:, :, k] = real.(ifft(w_1_coe[k, :, :]))
end

h_1[:, :] = real.(ifft(h_1_coe[1, :, :]))

w_1_T = (w_1[:, :, 2:end] + w_1[:, :, 1:end-1]) / 2.0

println("Inverse done.")
println("Output data...")

# Output data
axlist_T = (
    Dim{:x_T}(gd.x_T[:, 1, 1]),
    Dim{:y_T}(gd.y_T[1, :, 1]),
    Dim{:z_T}(gd.z_T[1, 1, :]),
)

axlist_U = (
    Dim{:x_U}(gd.x_U[:, 1, 1]),
    Dim{:y_T}(gd.y_T[1, :, 1]),
    Dim{:z_T}(gd.z_T[1, 1, :]),
)

axlist_V = (
    Dim{:x_T}(gd.x_T[:, 1, 1]),
    Dim{:y_V}(gd.y_V[1, :, 1]),
    Dim{:z_T}(gd.z_T[1, 1, :]),
)


axlist_W = (
    Dim{:x_T}(gd.x_T[:, 1, 1]),
    Dim{:y_T}(gd.y_T[1, :, 1]),
    Dim{:z_W}(gd.z_W[1, 1, :]),
)

axlist_sT = (
    Dim{:x_T}(gd.x_T[:, 1, 1]),
    Dim{:y_T}(gd.y_T[1, :, 1]),
)

axlist_cT = (
    Dim{:z_T}(gd.z_T[1, 1, :]),
)


data = Dict(
    :u_0   => YAXArray(axlist_cT,  u_0),
    :v_0   => YAXArray(axlist_cT,  v_0),
    :u_1   => YAXArray(axlist_U,   u_1),
    :v_1   => YAXArray(axlist_V,   v_1),
    :w_1   => YAXArray(axlist_W,   w_1),
    :h_1   => YAXArray(axlist_sT,  h_1),
    :pt_1  => YAXArray(axlist_sT,  Θ_1),
    :sst_1 => YAXArray(axlist_sT,  SST_1),
)

ds = Dataset(;data...)

filename = "test_solution_fourier.nc"
savedataset(ds,path = filename, driver=:netcdf, overwrite=true)

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

println("Plot solved 0-th order solution")
fig, ax = plt.subplots(1, 3)
ax[1].plot(u_0, s_T, "k-", label="u_0")
ax[2].plot(v_0, s_T, "k-", label="v_0")

ax[1].set_title("u_0")
ax[2].set_title("v_0")

ax[1].legend()
ax[2].legend()

ax[3].plot(u_0 / Ug, v_0 / Ug, "k-")
ax[3].scatter(u_0[1] / Ug, v_0[1] / Ug, s=10, c="red")
ax[3].set_title("Vector")

ax[3].grid(true)

fig.savefig("fig_zero.png", dpi=200)

println("Plot solved 1-th order solution (1D)")
fig, ax = plt.subplots(2, 1)

slb_size = (gd_slb.Nx, gd_slb.Ny)
SST_1_mid = reshape(SST_1, slb_size...)[:, Integer(round(gd_slb.Ny/2))]
Θ_1_mid = reshape(Θ_1, slb_size...)[:, Integer(round(gd_slb.Ny/2))]
h1_mid = reshape(h_1, slb_size...)[:, Integer(round(gd_slb.Ny/2))]

∇∇Θ_1_mid = reshape(amo_slb.T_LAPx_T * Θ_1[:], slb_size...)[:, Integer(round(gd_slb.Ny/2))]

x_T = gd_slb.x_T[:, 1, 1] / 1e3
y_T = gd_slb.y_T[1, :, 1] / 1e3

ax[1].plot(x_T, Θ_1_mid, "r-", label="\$\\Theta^{\\left(1\\right)}\$")
ax[1].plot(x_T, SST_1_mid, "b-", label="\$\\mathrm{SST}^{\\left(1\\right)}\$")

ax1_twinx = ax[1].twinx()
ax1_twinx.plot(x_T, h1_mid, "k--", label="\$ h^{\\left(1\\right)}\$")


ax1_twinx.set_ylabel("\$ h \$ [m]")
ax[1].set_ylabel("[K]")

ax[1].set_xlabel("x [km]")
ax[1].grid()

ax[1].legend(loc="upper left")
ax1_twinx.legend(loc="upper right")

ax[2].plot(x_T, SST_1_mid - Θ_1_mid, "b--", label="\$\\mathrm{SST}^{\\left(1\\right)} - \\Theta^{\\left(1\\right)} \$")
#ax[2].plot(x_T, ∇∇Θ_1_mid, "r-", label="\$ \\nabla^2 \\Theta^{\\left(1\\right)}\$")
ax[2].set_xlabel("x [km]")
ax[2].grid()

ax[2].legend()

fig.savefig("fig_mom_cx.png", dpi=200)


println("Plot solved 1-th order thermal solution (2D)")

fig, ax = plt.subplots(1, 2)

slb_size = (gd_slb.Nx, gd_slb.Ny)
levels = collect(range(-2, 2, length=41))

x_T = gd_slb.x_T[:, 1, 1] / 1e3
y_T = gd_slb.y_T[1, :, 1] / 1e3
cs = ax[1].contour(x_T, y_T, reshape(SST_1, slb_size...)', levels, colors="black")
plt.clabel(cs)

ax[1].quiver(x_T, y_T, u_1[:, :, 1]', v_1[:, :, 1]')

cs = ax[2].contour(x_T, y_T, reshape(Θ_1, slb_size...)', levels, colors="black")
plt.clabel(cs)


for _ax in ax
    _ax.set_xlabel("x [km]")
    _ax.set_xlabel("y [km]")
    _ax.grid()
end






ax[1].set_title("\$  \\mathrm{SST}^{\\left(1\\right)} \$")
ax[2].set_title("\$  \\Theta^{\\left(1\\right)} \$")

fig.savefig("fig_thermal.png", dpi=200)

println("Plot solved 1st order momentum solution (2D)")

fig, ax = plt.subplots(1, 2)

slb_size = (gd_slb.Nx, gd_slb.Ny)
levels = collect(range(-2, 2, length=41))

x_T = gd_slb.x_T[:, 1, 1] / 1e3
z_W = gd.z_W[1, 1, :]
z_T = gd.z_T[1, 1, :]
#cs = ax[1].contour(x_T, z_W, w_1[:, Integer(round(gd_slb.Ny/2)), :]', 20, colors="black")
cs = ax[1].contour(x_T, z_T, w_1_T[:, Integer(round(gd_slb.Ny/2)), :]', 20, colors="black")
plt.clabel(cs)

ax[2].plot(x_T, h_1[:, Integer(round(gd_slb.Ny/2))], color="black")

for _ax in ax
    _ax.set_xlabel("x [km]")
    _ax.grid()
end
    
ax[1].set_ylabel("s")
ax[2].set_ylabel("h_1 [m]")

ax[1].set_title("\$  w^{*\\left(1\\right)} \$")
ax[2].set_title("\$  h^{\\left(1\\right)} \$")


println("Save figure...")
fig.savefig("fig_mom.png", dpi=200)

println("Show plots...")
plt.show()
