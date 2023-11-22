using SparseArrays

include("BLM.jl")



m = BLM.Model(BLM.Env(
    f0 = 1e-4,
    Δx = 2e3,
    Δy = 2e3,
    Nx = 50,
    Ny = 1,
    Nz = 1000,
))


gd_col = m.ev.gd_col
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



# Construct operator
f0 = m.ev.f0
E0 = 1e-4
#τ_t = 0.0
#τ_b = 0.0
ug .= 5.0
vg .= 0.0
du0ds_0 = 1.0
dv0ds_0 = 0.5


# This matrix rotate the 2-D vector 90 degress clockwise
ROTCW = spdiagm(Ns => ones(Float64, Ns), - Ns => - ones(Float64, Ns) )

dUds_0_x[1] = du0ds_0
dUds_0_y[1] = dv0ds_0

amo = m.co.amo_col
op_1var = amo.T_DIVz_W * E0 * amo.W_mask_W * amo.W_∂z_T
op = blockdiag(op_1var, op_1var) + f0 * ROTCW

VEL_0[:] = op \ (f0 * ROTCW * VEL_g - blockdiag(amo.T_DIVz_W, amo.T_DIVz_W) * E0 * dUds_0)
#println(u_0 - ug)
#println(v_0 - vg)

# Idealized solution
k = sqrt(im * f0 / E0)
factor = exp( - 2 * k )
B = (du0ds_0 + dv0ds_0 * im) / ( k * ( factor - 1 ) )
A = B * factor
VEL_0_ana = A * exp.(k * s_T) + B * exp.(- k * s_T)

u_0_ana = real(VEL_0_ana) + ug
v_0_ana = imag(VEL_0_ana) + vg
              



println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(2, 1)
ax[1].plot(u_0     - ug, s_T, "k-", label="u_0")
ax[2].plot(v_0     - vg, s_T, "k-", label="v_0")
ax[1].plot(u_0_ana - ug, s_T, "r--", label="u_0_ana")
ax[2].plot(v_0_ana - vg, s_T, "r--", label="v_0_ana")

ax[1].set_title("u_0 - ug")
ax[2].set_title("v_0 - vg")

ax[1].legend()
ax[2].legend()

plt.show()
