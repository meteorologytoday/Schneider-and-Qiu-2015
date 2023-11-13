using SparseArrays
using LinearAlgebra, Arpack
using YAXArrays, DimensionalData, NetCDF
include("BLM.jl")



m = BLM.Model(BLM.Env(
    f0 = 1e-4,
    Δx = 5e3,
    Δy = 5e3,
    Nx = 100,
    Ny = 100,
    Nz = 10,
))


gd = m.ev.gd
gd_col = m.ev.gd_col
gd_slb = m.ev.gd_slb
Ns = gd_col.Nz
s_T = gd_col.z_T[:]

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

h_diffusivity = 1e2
ΔΘ = 1.0
Θ0 = 300.0
g0 = 10.0
reduced_g0 = g0 * ΔΘ / Θ0


# ===== Solving for 0st order =====
println("Solving for 0st order")
# Construct operator
f0 = m.ev.f0
E0 = 1e-5
dlnEdδ = 0.6
ug .= 20.0
vg .= 0.0
du0ds_0 = 0.1
dv0ds_0 = 0.00
s0 = 0.05
γ0 = 0.3

# E_0_W is the zero-th order E, while E0 is a constant
function getE_0(z)
    tmp =  (z .+ s0) / γ0
    return tmp .* E0 .* exp.( 1.0 .- tmp ) #* 0 .+ E0
end

function getE_1(z, δ)
#    E_0 = getE_0(z)
#    print(size(δ))
#    print(size(E_0))
    return δ .* dlnEdδ .* getE_0(z) 
end



E_0_W = getE_0(gd_col.z_W[:])
W_E_0_W = spdiagm(0=>E_0_W)


# This matrix rotate the 2-D vector 90 degress clockwise
ROTCW = spdiagm(Ns => ones(Float64, Ns), - Ns => - ones(Float64, Ns) )

dUds_0_x[1] = du0ds_0
dUds_0_y[1] = dv0ds_0


amo_col = m.co.amo_col
op_1var = amo_col.T_DIVz_W * W_E_0_W * amo_col.W_mask_W * amo_col.W_∂z_T
op = blockdiag(op_1var, op_1var) + f0 * ROTCW

println("Solving 0th order")
@time F = lu(op)
@time VEL_0[:] = F \ (f0 * ROTCW * VEL_g - blockdiag(amo_col.T_DIVz_W, amo_col.T_DIVz_W) * blockdiag(W_E_0_W, W_E_0_W) * dUds_0)

# Idealized solution
k = sqrt(im * f0 / E0)
factor = exp( - 2 * k )
B = (du0ds_0 + dv0ds_0 * im) / ( k * ( factor - 1 ) )
A = B * factor
VEL_0_ana = A * exp.(k * s_T) + B * exp.(- k * s_T)

u_0_ana = real(VEL_0_ana) + ug
v_0_ana = imag(VEL_0_ana) + vg
              



# ===== Solving for 1st order =====
println("Solving for 1st order")
amo_slb = m.co.amo_slb
bmo_slb = amo_slb.bmo
h_0 = 1e3 # m
A_h = 1e5 # m^2 / s
γ_Θ = 0.1
ΔSST  = 1.0 # K
# Compute avg(u)
u_0_mean = sum(u_0 .* gd_col.Δz_T[:])
v_0_mean = sum(v_0 .* gd_col.Δz_T[:])

VEL_0_mean = zeros(Float64, amo_slb.bmo.U_pts * 2)

SST_1 = zeros(Float64, amo_slb.bmo.T_pts)


x_c = ( gd_slb.x_U[end, 1, 1] - gd_slb.x_U[1, 1, 1] ) / 2.0
y_c = ( gd_slb.y_V[1, end, 1] - gd_slb.y_V[1, 1, 1] ) / 2.0

SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 + (gd_slb.y_T .- y_c).^2 ) / 50e3^2 )
#SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 ) / 50e3^2 )
#u_0_mean = view(VEL_0, 1:Ns)
#v_0_mean = view(VEL_0, (Ns+1):(2*Ns))
# Construct operator
op = (
       amo_slb.T_interp_U * u_0_mean * amo_slb.U_∂x_T
    .+ γ_Θ / h_0 * amo_slb.bmo.T_I_T
    - (amo_slb.T_DIVx_U * A_h * amo_slb.U_mask_U * amo_slb.U_∂x_T)
    - (amo_slb.T_DIVy_V * A_h * amo_slb.V_mask_V * amo_slb.V_∂y_T)
)

# Solve with boundary condition ∇Θ_1 = 0
println("Solving 1th order thermal")
@time F = lu(op)
@time Θ_1 = F \ ( γ_Θ / h_0 * SST_1)

# Now I can try to solve for u, w, and h


# 1: construct E(1) profile
δ_1 = SST_1 - Θ_1 # sT grid


# convert δ_1 to δ_1_W

#=
δ_1_W = repeat(δ_1, gd_col.Nz+1)
E_1_W = δ_1_W * dlnEdδ * E_0_W
=#

# 2: constrcut effective grid and sender matrices
# find effective grid of vertical velocities
# Vertical velocity: top and bottom are zeros

amo = m.co.amo
bmo = amo.bmo

T_pts = bmo.T_pts
U_pts = bmo.U_pts
V_pts = bmo.V_pts
W_pts = bmo.W_pts
sT_pts = amo_slb.bmo.T_pts
sU_pts = amo_slb.bmo.U_pts
sV_pts = amo_slb.bmo.V_pts

mask_eff_W = reshape(amo.W_mask_W * ones(Float64, W_pts), amo.bmo.W_dim...)
# Create coversion matrix and its inverse

# w: the top and bottom w^*(1) are not free variables
W_num = reshape(collect(1:length(mask_eff_W)), size(mask_eff_W)...)
active_num_eff_W = W_num[ mask_eff_W .==1 ]
eW_send_W = amo.bmo.W_I_W[active_num_eff_W, :]; dropzeros!(eW_send_W)
W_send_eW = sparse(eW_send_W')

# h: pick one grid to be the balance
sT_h_1_bc = - ones(Float64, sT_pts-1)
sT_h_1_bc[end] = 0.0

sT_send_h_1_esT = [   sparse(I, sT_pts - 1, sT_pts-1) ;
                    - ones(Float64, 1,    sT_pts-1) ; ]


#println(Matrix(sT_send_h_1_esT))

#eff_send_full = blockdiag(sparse(I, U_pts, U_pts), sparse(I, V_pts, V_pts), eW_send_W, sparse(I, sT_pts, sT_pts) )
#full_send_eff = eff_send_full'


full_send_eff = blockdiag(sparse(I, T_pts, T_pts), sparse(I, T_pts, T_pts), W_send_eW, sT_send_h_1_esT)

println("Size of full_send_eff: ", size(full_send_eff))

# 3: formulate operator

u_0_T = repeat(reshape(u_0, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))
v_0_T = repeat(reshape(v_0, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))

T_u_0_T = spdiagm(0 => u_0_T[:] )
T_v_0_T = spdiagm(0 => v_0_T[:] )

u_0_U = repeat(reshape(u_0, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))
v_0_V = repeat(reshape(v_0, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))

U_u_0_U = spdiagm(0 => u_0_U[:] )
V_v_0_V = spdiagm(0 => v_0_V[:] )



println("1: SIZE of amo_slb.U_∂x_T ", size(repeat(amo_slb.U_∂x_T, outer=(gd.Nz, 1))))
println("2: SIZE of amo.U_∂x_T ",     size(amo.U_∂x_T))

println("sT_pts : ", sT_pts)
println("T_pts  : ", T_pts)
println("U_pts  : ", U_pts)
println("V_pts  : ", V_pts)
println("W_pts  : ", W_pts)

op_cont = sparse((
      h_0 * [ (amo.T_DIVx_U * amo.U_interp_T)  (amo.T_DIVy_V * amo.V_interp_T)  spzeros(Float64, T_pts, W_pts + sT_pts) ]
      + [ (amo.T_interp_U * U_u_0_U)  (amo.T_interp_V * V_v_0_V) ]  
      * [ spzeros(Float64, T_pts, T_pts + T_pts + W_pts)   repeat(amo_slb.U_∂x_T, outer=(gd.Nz, 1)) ; 
          spzeros(Float64, T_pts, T_pts + T_pts + W_pts)   repeat(amo_slb.V_∂y_T, outer=(gd.Nz, 1)) ; ]
      + [ spzeros(Float64, T_pts, T_pts + T_pts)  amo.T_DIVz_W  spzeros(Float64, T_pts, sT_pts) ]
      - [ spzeros(Float64, T_pts, T_pts + T_pts + W_pts)   h_diffusivity * ( repeat(amo_slb.T_DIVx_U * amo_slb.U_∂x_T + amo_slb.T_DIVy_V * amo_slb.V_∂y_T, outer=(gd.Nz, 1)) ) ]
))

println("Size of op_cont: ", size(op_cont))


println("Size of bmo.T_I_T: ", size(bmo.T_I_T))

op_mom_cori = m.ev.f0 * sparse([
    spzeros(Float64, T_pts, T_pts)  (- bmo.T_I_T)                   spzeros(Float64, T_pts, W_pts + sT_pts) ;
    bmo.T_I_T                       spzeros(Float64, T_pts, T_pts)  spzeros(Float64, T_pts, W_pts + sT_pts) ; 
])

println("Size of op_mom_cori: ", size(op_mom_cori))


#=
op_mom_adv   = sparse(
    [ bmo.U_I_U  amo.U_interp_V  spzeros(Float64, U_pts, U_pts)  spzeros(Float64, U_pts, V_pts) ;
      spzeros(Float64, V_pts, U_pts)  spzeros(Float64, V_pts, V_pts)  amo.V_interp_U  bmo.V_I_V ; ] * 
    [ U_u_0_U * amo.U_∂x_U                   spzeros(Float64, U_pts, V_pts)   spzeros(Float64, U_pts, W_pts + sT_pts) ;
      V_v_0_V * amo.V_∂y_T * amo.T_interp_U  spzeros(Float64, U_pts, V_pts)   spzeros(Float64, V_pts, W_pts + sT_pts) ;
      spzeros(Float64, U_pts, U_pts)   U_u_0_U * amo.U_∂x_T * amo.T_interp_V  spzeros(Float64, U_pts, W_pts + sT_pts) ;
      spzeros(Float64, U_pts, U_pts)   V_v_0_V * amo.V_∂y_V                   spzeros(Float64, U_pts, W_pts + sT_pts) ; ]
)
=#
    
tmp = (T_u_0_T * amo.T_interp_U * amo.U_∂x_T)  +  (T_v_0_T * amo.T_interp_V * amo.V_∂y_T)
println("Size of tmp: ", size(tmp))
op_mom_adv = [ blockdiag(tmp, tmp)  spzeros(Float64, 2*T_pts, W_pts + sT_pts) ]
println("Size of op_mom_adv: ", size(op_mom_adv))

# Derive bottom Neumann boundary conditions

τx_1_0 = -1.0
τy_1_0 = 1.0
u_0_W = amo_col.W_interp_T * u_0[:]
v_0_W = amo_col.W_interp_T * v_0[:]

du_1ds_bc_W = zeros(Float64, amo.bmo.W_dim...)
dv_1ds_bc_W = zeros(Float64, amo.bmo.W_dim...)

# E_0_0 = E_0 at s=0
E_0_0 = getE_0(zeros(Float64, gd.Nx, gd.Ny))
E_1_0 = getE_1(zeros(Float64, gd.Nx, gd.Ny), reshape(δ_1, gd.Nx, gd.Ny))


du_1ds_bc_W[:, :, 1] = ( τx_1_0 .- E_1_0 * du0ds_0 ) ./ E_0_0
dv_1ds_bc_W[:, :, 1] = ( τy_1_0 .- E_1_0 * dv0ds_0 ) ./ E_0_0

rmbot_W = ones(Float64, bmo.W_dim...)
rmtop_W = ones(Float64, bmo.W_dim...)

rmbot_W[:, :, 1]   .= 0.0
rmtop_W[:, :, end] .= 0.0

W_rmbot_W = spdiagm(0=>view(rmbot_W, :))
W_rmtop_W = spdiagm(0=>view(rmtop_W, :))

op_mom_vdif = - ( blockdiag(amo.T_DIVz_W, amo.T_DIVz_W)
     * blockdiag(W_rmbot_W * amo.W_∂z_W * W_rmtop_W, 
                 W_rmbot_W * amo.W_∂z_W * W_rmtop_W,)
     * [ blockdiag(amo.W_interp_T, amo.W_interp_T)  spzeros(Float64, 2 * W_pts, W_pts + sT_pts) ]
)


T_∂x_sT = repeat(sparse(I, sT_pts, sT_pts), outer=(gd.Nz, 1)) * amo_slb.T_interp_U * amo_slb.U_∂x_T
T_∂y_sT = repeat(sparse(I, sT_pts, sT_pts), outer=(gd.Nz, 1)) * amo_slb.T_interp_V * amo_slb.V_∂y_T
op_mom_reduced_gravity = reduced_g0 * [ spzeros(Float64, 2*T_pts, 2*T_pts + W_pts)  [ T_∂x_sT ; T_∂y_sT ; ] ]

# building vertical advection operator
duds_0_cW = amo_col.W_∂z_T * u_0[:]
duds_0_cW[1] = du0ds_0
duds_0_cW[end] = 0.0
duds_0_cT = amo_col.T_interp_W * duds_0_cW

duds_0_T = repeat( reshape(duds_0_cT, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))[:]
T_duds_0_T = spdiagm(0=>duds_0_T)

dvds_0_cW = amo_col.W_∂z_T * v_0[:]
dvds_0_cW[1] = dv0ds_0
dvds_0_cW[end] = 0.0
dvds_0_cT = amo_col.T_interp_W * dvds_0_cW

dvds_0_T = repeat( reshape(dvds_0_cT, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))[:]
T_dvds_0_T = spdiagm(0=>dvds_0_T)

op_mom_vadv = blockdiag(T_duds_0_T, T_dvds_0_T) * [ 
    spzeros(Float64, 2 * T_pts, 2 * T_pts)  [ amo.T_interp_W ; amo.T_interp_W ; ] spzeros(2 * T_pts, sT_pts)
] 

println("Size of op_cont: ", size(op_cont))
println("Size of op_mom_adv: ", size(op_mom_adv))
println("Size of op_mom_vadv: ", size(op_mom_vadv))
println("Size of op_mom_cori: ", size(op_mom_cori))
println("Size of op_mom_vdif: ", size(op_mom_vdif))
println("Size of op_mom_reduced_gravity: ", size(op_mom_reduced_gravity))



#op_mom = op_mom_adv + op_mom_vadv + op_mom_cori + op_mom_vdif + op_mom_reduced_gravity
#op_mom = op_mom_cori + op_mom_vdif + op_mom_reduced_gravity
op_mom = op_mom_cori + op_mom_reduced_gravity
#op_mom = op_mom_cori# + op_mom_vdif )
println("Size of op_mom: ", size(op_mom))

# Construct right-hand-side

hgt_factor = spdiagm(0 => (1.0 .- gd.z_T[:]))
RHS_mom_pt = (reduced_g0 * h_0 / ΔΘ) * [ (hgt_factor * T_∂x_sT) ; (hgt_factor * T_∂y_sT) ; ] * Θ_1

duds_0_W = repeat( reshape(duds_0_cW, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))[:]
dvds_0_W = repeat( reshape(dvds_0_cW, 1, 1, :), outer=(gd.Nx, gd.Ny, 1))[:]
W_E_1_W = spdiagm(0=> getE_1(gd.z_W, repeat( reshape(δ_1, gd.Nx, gd.Ny, 1), outer=(1, 1, gd.Nz+1) ))[:] )

RHS_mom_vdif = ( 
    blockdiag(amo.T_DIVz_W, amo.T_DIVz_W)
  * blockdiag(W_E_1_W, W_E_1_W)
  * [ duds_0_W ; dvds_0_W ; ]
)

RHS_mom_vdif_bc = 0*(
    blockdiag(amo.T_DIVz_W, amo.T_DIVz_W) * [ du_1ds_bc_W[:] ; dv_1ds_bc_W[:] ; ]
)

RHS = [ (RHS_mom_pt + RHS_mom_vdif + RHS_mom_vdif_bc) ; zeros(Float64, T_pts) ; ]

RHS = RHS[:]

# 4: solve for solution

op_cont = op_cont[1:end-1, :]
RHS = RHS[1:end-1]

op = [ op_mom ; op_cont ; ] * full_send_eff
println("Size of op: ", size(op))
println("Size of RHS: ", size(RHS))

println("Solving 1th order momentum and height")



@time F = lu(op)
@time sol = full_send_eff * (F \ RHS)

# Remove unstable modes
#egvals, egvecs = eigs(op, nev=5, ncv=100, check=1)
#println("Eigenvalues: ", egvals)


# 5: project back to original grid 



# Output data

idx=1

u_1 = reshape(view(sol, idx:(idx+U_pts-1)),  bmo.U_dim...);     idx += U_pts
v_1 = reshape(view(sol, idx:(idx+V_pts-1)),  bmo.V_dim...);     idx += V_pts
w_1 = reshape(view(sol, idx:(idx+W_pts-1)),  bmo.W_dim...);     idx += W_pts
#h_1 = reshape(amo_slb.T_NEWSavg_T * view(sol, idx:(idx+sT_pts-1)), bmo_slb.T_dim[1:2]...); idx += sT_pts
h_1 = reshape(view(sol, idx:(idx+sT_pts-1)), bmo_slb.T_dim[1:2]...); idx += sT_pts
Θ_1 = reshape(Θ_1, bmo_slb.T_dim[1:2]...)
SST_1 = reshape(SST_1, bmo_slb.T_dim[1:2]...)
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
    :u_1   => YAXArray(axlist_U,  u_1),
    :v_1   => YAXArray(axlist_V,  v_1),
    :w_1   => YAXArray(axlist_W,  w_1),
    :h_1   => YAXArray(axlist_sT, h_1),
    :pt_1  => YAXArray(axlist_sT, Θ_1),
    :sst_1 => YAXArray(axlist_sT, SST_1),
)

ds = Dataset(;data...)

filename = "test_solution.nc"
savedataset(ds,path = filename, driver=:netcdf, overwrite=true)

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

println("Plot solved 0-th order solution")
fig, ax = plt.subplots(2, 1)
ax[1].plot(u_0     - ug, s_T, "k-", label="u_0")
ax[2].plot(v_0     - vg, s_T, "k-", label="v_0")
ax[1].plot(u_0_ana - ug, s_T, "r--", label="u_0_ana")
ax[2].plot(v_0_ana - vg, s_T, "r--", label="v_0_ana")

ax[1].set_title("u_0 - ug")
ax[2].set_title("v_0 - vg")

ax[1].legend()
ax[2].legend()

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
cs = ax[1].contour(x_T, z_W, w_1[:, Integer(round(gd_slb.Ny/2)), :]', 20, colors="black")
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


print("Save figure...")
fig.savefig("fig_mom.png", dpi=200)

print("Show plots...")
plt.show()
