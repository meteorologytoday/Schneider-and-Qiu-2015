using SparseArrays
using LinearAlgebra

include("BLM.jl")



m = BLM.Model(BLM.Env(
    f0 = 1e-4,
    Δx = 5e3,
    Δy = 5e3,
    Nx = 100,
    Ny = 100,
    Nz = 2,
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



# ===== Solving for 0st order =====
println("Solving for 0st order")
# Construct operator
f0 = m.ev.f0
E0 = 1e-4
ug .= 5.0
vg .= 0.0
du0ds_0 = 1.0
dv0ds_0 = 0.0

# This matrix rotate the 2-D vector 90 degress clockwise
ROTCW = spdiagm(Ns => ones(Float64, Ns), - Ns => - ones(Float64, Ns) )

dUds_0_x[1] = du0ds_0
dUds_0_y[1] = dv0ds_0

amo_col = m.co.amo_col
op_1var = amo_col.T_DIVz_W * E0 * amo_col.W_mask_W * amo_col.W_∂z_T
op = blockdiag(op_1var, op_1var) + f0 * ROTCW

F = lu(op)
VEL_0[:] = F \ (f0 * ROTCW * VEL_g - blockdiag(amo_col.T_DIVz_W, amo_col.T_DIVz_W) * E0 * dUds_0)

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

SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 + (gd_slb.y_T .- y_c).^2 ) / 100e3^2 )
#u_0_mean = view(VEL_0, 1:Ns)
#v_0_mean = view(VEL_0, (Ns+1):(2*Ns))
# Construct operator
op = (
       amo_slb.T_interp_U * u_0_mean * amo_slb.U_∂x_T
    .+ γ_Θ / h_0 * amo_slb.bmo.T_I_T
    - (amo_slb.T_DIVx_U * A_h * amo_slb.U_mask_U * amo_slb.U_∂x_T)
    #- (amo_slb.T_DIVy_V * A_h * amo_slb.V_mask_V * amo_slb.V_∂y_T)
)

# Solve with boundary condition ∇Θ_1 = 0
F = lu(op)
Θ_1 = F \ ( γ_Θ / h_0 * SST_1)


println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

#=
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
=#

println("Plot solved 1-th order solution (1D)")
fig, ax = plt.subplots(2, 1)

slb_size = (gd_slb.Nx, gd_slb.Ny)
SST_1_mid = reshape(SST_1, slb_size...)[:, Integer(round(gd_slb.Ny/2))]
Θ_1_mid = reshape(Θ_1, slb_size...)[:, Integer(round(gd_slb.Ny/2))]

∇∇Θ_1_mid = reshape(amo_slb.T_LAPx_T * Θ_1, slb_size...)[:, Integer(round(gd_slb.Ny/2))]

x_T = gd_slb.x_T[:, 1, 1] / 1e3
y_T = gd_slb.y_T[1, :, 1] / 1e3

ax[1].plot(x_T, Θ_1_mid, "r-", label="\$\\Theta^{\\left(1\\right)}\$")
ax[1].plot(x_T, SST_1_mid, "b-", label="\$\\mathrm{SST}^{\\left(1\\right)}\$")

ax[1].set_xlabel("x [km]")
ax[1].grid()

ax[1].legend()

ax[2].plot(x_T, SST_1_mid - Θ_1_mid, "b--", label="\$\\mathrm{SST}^{\\left(1\\right)} - \\Theta^{\\left(1\\right)} \$")
#ax[2].plot(x_T, ∇∇Θ_1_mid, "r-", label="\$ \\nabla^2 \\Theta^{\\left(1\\right)}\$")
ax[2].set_xlabel("x [km]")
ax[2].grid()

ax[2].legend()


println("Plot solved 1-th order solution (2D)")

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

plt.show()
