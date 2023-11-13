include("ZAAM.jl")

m = ZAAM.Model(ZAAM.Env(
    ϕn =  deg2rad(90),
    ϕs =  deg2rad(-90),
    npts = 51, # Need to be odd number so that V grid does not have ϕ=0
    J = 1e3,
    L = 1e5 * 0,
    K = 1e5,
    H = 10e3, 
    N = 1e-6 * 0,
    C = ZAAM.ρ0 * ZAAM.c_p * 1e-3 * 10.0,
))

pp = m.ev.pp
gd = m.ev.gd

ϕ_V = gd.ϕ_V[:]
f_V = 2 * gd.Ω * sin.(ϕ_V)


# Test if Legendre polynomial works

P = [
    (x,) -> x,
    (x,) -> (3*x.^2 .- 1) / 2,
    (x,) -> (5*x.^3 .- 3 * x) / 2,
    (x,) -> (35*x.^4 .- 30*x.^2 .+ 3) / 8,
]

N = 4

Γ    = zeros(Float64, length(m.ev.gd.ϕ_V), N)
∂B∂y = zeros(Float64, length(m.ev.gd.ϕ_V), N)
ana_Ψ = zeros(Float64, length(m.ev.gd.ϕ_V), N)
num_Ψ = zeros(Float64, length(m.ev.gd.ϕ_V), N)

for n=1:N

    # We assum Γ is a legendre polynomial Pn(sin(ϕ))

    Pn = P[n](sin.(ϕ_V))
    λ = n * (n+1)
    
    Γ[:, n] = Pn
    ∂B∂y[:, n] = - f_V.^(-1) .* ( (f_V.^2 .+ pp.J^2 * pp.k^4) .+ n*(1+n)*pp.J*pp.L*pp.k^2/gd.R^2 ) .* Pn 
    
    num_Ψ[:, n] = m.co.psi_solver.V_solveΨ_from∂B∂y_V * ∂B∂y[:, n]
    ana_Ψ[:, n] = - (pp.J * pp.k^4)^(-1) * ( ∂B∂y[:, n] + f_V .* Γ[:, n] )

end





println("Try to output data...")

println("Loading libraries...")
using YAXArrays, NetCDF
using Dates
println("Done")

axlist = [
    RangeAxis("yp1", collect(1:length(m.ev.gd.ϕ_V))),
    RangeAxis("legendre", collect(1:N)),
]
da_Γ     = YAXArray(axlist, Γ)
da_∂B∂y  = YAXArray(axlist, ∂B∂y)
da_ana_Ψ = YAXArray(axlist, ana_Ψ)
da_num_Ψ = YAXArray(axlist, num_Ψ)
ds = Dataset(Gamma=da_Γ, dbdy=da_∂B∂y, ana_Psi=da_ana_Ψ, num_Psi=da_num_Ψ,)

savedataset(ds, path="test_solve_psi.nc", driver=:netcdf, overwrite=true)


