include("../../src/BLM.jl")
using Formatting
using Interpolations
using YAXArrays, DimensionalData, NetCDF

target_file = "wrfout_avg.nc"
output_file = "output2.nc"


selected_x_idx = 1
selected_y_idx = 1
Ug = 15.0
Vg = 0.0
f0 = 1e-4

E0_edges = range(0, 1, length=101) * 1e-5
E0_mids = ( E0_edges[2:end] + E0_edges[2:end] ) / 2

prior = ones(Float64, length(E0_mids))
loglikelihood = zeros(Float64, length(E0_mids))

prior[E0_mids .<= 0] .= 0.0

σ_vel = 1e-1

println("Loading netcdf")
ds = open_dataset(target_file)
z = collect((ds.PH + ds.PHB) ./ 9.81)[selected_x_idx, selected_y_idx, :, 1]
z = (z[2:end] + z[1:end-1])/2
h = collect(ds.PBLH)[selected_x_idx, selected_y_idx, 1]
obs_s = z ./ h

obs_U = collect(ds.U)
obs_U = (obs_U[2:end, :, :, :] + obs_U[1:end-1, :, :, :]) / 2
obs_U = obs_U[selected_x_idx, selected_y_idx, :, 1]

obs_V = collect(ds.V)[selected_x_idx, selected_y_idx, :, 1]

obs_s = convert(Array{Float64}, obs_s)
obs_U = convert(Array{Float64}, obs_U)
obs_V = convert(Array{Float64}, obs_V)

idx = obs_s .<= 1.0
obs_s = obs_s[idx]
obs_U = obs_U[idx]
obs_V = obs_V[idx]


println("size(obs_U) = ", size(obs_U))
println("size(obs_V) = ", size(obs_V))
println("size(obs_s) = ", size(obs_s))
println("s = ", obs_s)


println("Done")

println("Instantiating a boundary layer model...")
m = BLM.Model(BLM.Env(
    Δx = 5e3,
    Δy = 5e3,
    Nx = 3,
    Ny = 3,
    Nz = 101,
    f0 = f0,
    ΔΘ = 10.0,
    Θ0 = 290.0,
    g0 = 9.81,
    h_0 = 1.5e3,   # h_0 does not matter in solving Ekman balance
    A_h = 1e5,
    γ_Θ = 0.25,
    E0  = 2e-4,
    dlnγdδ = 0.5,
    γ0     = 0.3,
))
println("Done.")


function computeLoglikelihood(;
    obs_s :: AbstractArray{Float64, 1},
    obs_U :: AbstractArray{Float64, 1},
    obs_V :: AbstractArray{Float64, 1}, 
    E0    :: Float64,
    return_solution :: Bool = false,
)


    # Setup E0
    m.ev.pp.E0  = E0

    # Setup geostrophic wind (background pressure gradient)
    m.st.u_g .= Ug
    m.st.v_g .= Vg

    BLM.solve_order0!(m)

    # Interpolate data to the same grid
    # BLM grid to WRF grid
    interp_linear_U = linear_interpolation(m.ev.gd.z_T[1, 1, :], m.st.u_0)
    interp_linear_V = linear_interpolation(m.ev.gd.z_T[1, 1, :], m.st.v_0)

    u_0 = interp_linear_U.(obs_s)
    v_0 = interp_linear_V.(obs_s)

    diff_u = - (1/2) * ((u_0 - obs_U) / σ_vel).^2
    diff_v = - (1/2) * ((v_0 - obs_V) / σ_vel).^2

    loglikelihood = sum(diff_u) + sum(diff_v)

    if !return_solution
        return loglikelihood
    else
        return loglikelihood, u_0, v_0
    end

end

if isfile(output_file)
    println("File $output_file already exist. Skip solver.")
else
    
    println("File $output_file does not exist. Solve the model now...")
    
    for (i, E0) in enumerate(E0_mids)
        
        loglikelihood[i] = computeLoglikelihood(
            obs_s = obs_s,
            obs_U = obs_U,
            obs_V = obs_V,
            E0    = E0,
        )
        
    end

    loglikelihood .-= maximum(loglikelihood)
    post = exp.(loglikelihood) .* prior

    (_, idx) = findmax(loglikelihood)    
    E0_picked = E0_mids[idx]

    _, u_0, v_0 = computeLoglikelihood(
        obs_s = obs_s,
        obs_U = obs_U,
        obs_V = obs_V,
        E0    = E0_picked;
        return_solution = true,
    )




    # Save output
    axlist_E0 = (
        Dim{:E0}(E0_mids),
    )

    # Save output
    axlist_windprofile = (
        Dim{:s_T}(obs_s),
    )


    data = Dict(
        :loglikelihood => YAXArray(axlist_E0, loglikelihood),
        :post          => YAXArray(axlist_E0, post),
        :u_0           => YAXArray(axlist_windprofile, u_0),
        :v_0           => YAXArray(axlist_windprofile, v_0),
        :obs_U         => YAXArray(axlist_windprofile, obs_U),
        :obs_V         => YAXArray(axlist_windprofile, obs_V),
    ) 

    props = Dict(
        "E0_maxlikelihood" => E0_picked,
    )

    ds = Dataset(;data..., properties=props)

    println("Saving model into file: ", output_file)

    savedataset(ds, path = output_file, driver=:netcdf, overwrite=true)

end




