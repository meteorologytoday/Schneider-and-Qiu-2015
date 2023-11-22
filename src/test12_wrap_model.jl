using SparseArrays
using LinearAlgebra, Arpack
using NetCDF, YAXArrays, DimensionalData
include("BLM.jl")
include("tools.jl")


output_filename = "solved_model.nc"
Ug = 15.0
Vg = 0.0


if isfile(output_filename)

    println("File $output_filename already exist. Skip solver.")

else

    println("File $output_filename does not exist. Solve the model now...")
    println("Instantiating a boundary layer model...")
    m = BLM.Model(BLM.Env(
        Δx = 5e3,
        Δy = 5e3,
        Nx = 50,
        Ny = 50,
        Nz = 10,
        f0 = 1e-4,
        ΔΘ = 10.0,
        Θ0 = 290.0,
        g0 = 9.81,
        h_0 = 1e3,
        A_h = 1e5,
        γ_Θ = 0.25,
        E0  = 2.5e-5,
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
    ΔSST  = 1.0 # K
    x_c = Lx / 2.0
    y_c = Ly / 2.0
    #st.SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 + (gd_slb.y_T .- y_c).^2 ) / 50e3^2 )
    st.SST_1[:] = ΔSST * exp.( - ( (gd_slb.x_T .- x_c).^2 ) / 50e3^2 )



    println("Now solve the 0-th order")
    BLM.solve_order0!(m)
    println("Done solving the 0-th order")

    println("Now solve the 1-th order thermal")
    BLM.solve_order1_thermal!(m; wn_rad=wn_rad)
    println("Done solving the 1-th order thermal")

    println("Now solve the 1-th order momentum")
    BLM.solve_order1_momentum!(m; wn_rad=wn_rad)
    println("Done solving the 1-th order momentum")

    BLM.savedata(m, output_filename, overwrite=true)

    m = Nothing

end

# clear model

println("Loading file: ", output_filename)

ds = open_dataset(output_filename)
x_T = collect(ds.x_T) ./ 1e3
y_T = collect(ds.y_T) ./ 1e3
s_T = collect(ds.z_T) 

println(size(s_T))
println(size(ds.u_0))

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

println("Plot solved 0-th order solution")
fig, ax = plt.subplots(1, 3)
ax[1].plot(ds.u_0, s_T, "k-", label="u_0")
ax[2].plot(ds.v_0, s_T, "k-", label="v_0")

ax[1].set_title("u_0")
ax[2].set_title("v_0")

ax[1].legend()
ax[2].legend()

ax[3].plot(ds.u_0, ds.v_0, "k-")
ax[3].scatter(ds.u_0[1], ds.v_0[1], s=10, c="red")
ax[3].set_title("Vector")

ax[3].grid(true)

fig.savefig("fig_zero.png", dpi=200)


println("Plot solved 1-th order thermal solution (2D)")

fig, ax = plt.subplots(1, 2)

levels = collect(range(-2, 2, length=41))

cs = ax[1].contour(x_T, y_T, ds.sst_1', levels, colors="black")
plt.clabel(cs)

cs = ax[1].contour(x_T, y_T, ds.w_1[:, :, 2]', 10, colors="red")
plt.clabel(cs)

#u = zeros(Float64, 1, 1) .+ st.u_0_mean
#v = zeros(Float64, 1, 1) .+ st.v_0_mean
#ax[1].quiver([x_c/1e3,], [y_c/1e3,], u', v')

ax[1].quiver(x_T, y_T, ds.u_1[:, :, 1]', ds.v_1[:, :, 1]')

cs = ax[2].contour(x_T, y_T, ds.pt_1', levels, colors="black")


plt.clabel(cs)

for _ax in ax
    _ax.set_xlabel("x [km]")
    _ax.set_xlabel("y [km]")
    _ax.grid()
end


ax[1].set_title("\$  \\mathrm{SST}^{\\left(1\\right)} \$")
ax[2].set_title("\$  \\Theta^{\\left(1\\right)} \$")

println("Showing figures...")
fig.savefig("fig_2D.png", dpi=200)
plt.show()

