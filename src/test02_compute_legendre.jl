include("ZAAM.jl")

m = ZAAM.Model(ZAAM.Env(
    ϕn =  deg2rad(90),
    ϕs =  deg2rad(-90),
    npts = 51,
    J = 1e3,
    L = 1e5 * 0,
    K = 1e5,
    H = 10e3, 
    N = 1e-6 * 0,
    C = ZAAM.ρ0 * ZAAM.c_p * 1e-3 * 10.0,
))


# Test if Legendre polynomial works

# From https://en.wikipedia.org/wiki/Legendre_polynomials
P = [
    (x,) -> x,
    (x,) -> (3*x.^2 .- 1) / 2,
    (x,) -> (5*x.^3 .- 3 * x) / 2,
    (x,) -> (35*x.^4 .- 30*x.^2 .+ 3) / 8,
]

N = length(P)

Γ = zeros(Float64, length(m.ev.gd.ϕ_V), N)
ana_DIVyΓ = zeros(Float64, length(m.ev.gd.ϕ_V), N)
num_DIVyΓ = zeros(Float64, length(m.ev.gd.ϕ_V), N)

for n=1:N
    Γ[:, n] = P[n](sin.(m.ev.gd.ϕ_V))
    λ = n * (n+1)
    
    expected_ratio = - λ / m.ev.gd.R^2

    ana_DIVyΓ[:, n] = expected_ratio * Γ[:, n]
    num_DIVyΓ[:, n] = m.co.amo.V_LAPy_V * Γ[:, n]
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
da_Γ  = YAXArray(axlist, Γ)
da_ana_DIVyΓ = YAXArray(axlist, ana_DIVyΓ)
da_num_DIVyΓ = YAXArray(axlist, num_DIVyΓ)
ds = Dataset(Gamma=da_Γ, ana_DIVyGamma=da_ana_DIVyΓ, num_DIVyGamma=da_num_DIVyΓ)

savedataset(ds, path="test_compute_legendre.nc", driver=:netcdf, overwrite=true)


