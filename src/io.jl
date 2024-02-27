
function savedata(m :: Model, filename :: String; overwrite=false, props :: Dict = Dict())

    st = m.st
    gd = m.ev.gd

    # Output data

    axlist_axisX = (
        Dim{:x_T}(gd.x_T[:, 1, 1]),
    )

    axlist_axisY = (
        Dim{:y_T}(gd.y_T[1, :, 1]),
    )

    axlist_axisZ = (
        Dim{:z_T}(gd.z_T[1, 1, :]),
    )

    axlist_axisX_U = (
        Dim{:x_U}(gd.x_U[:, 1, 1]),
    )

    axlist_axisY_V = (
        Dim{:y_V}(gd.y_V[1, :, 1]),
    )

    axlist_axisZ_W = (
        Dim{:z_W}(gd.z_W[1, 1, :]),
    )


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
        :u_0   => YAXArray(axlist_cT,  st.u_0),
        :v_0   => YAXArray(axlist_cT,  st.v_0),
        :u_1   => YAXArray(axlist_T,   st.u_1),
        :v_1   => YAXArray(axlist_T,   st.v_1),
        :w_1   => YAXArray(axlist_W,   st.w_1),
        :h_1   => YAXArray(axlist_sT,  st.h_1),
        :pt_1  => YAXArray(axlist_sT,  st.Θ_1),
        :sst_1 => YAXArray(axlist_sT,  st.SST_1),
        :dx_T  => YAXArray(axlist_axisX,  gd.Δx_T[:, 1, 1]),
        :dy_T  => YAXArray(axlist_axisY,  gd.Δy_T[1, :, 1]),
        :dz_T  => YAXArray(axlist_axisZ,  gd.Δz_T[1, 1, :]),
    )
    
    pp = m.ev.pp
    props = merge(
        props,
        Dict(
            "f0"             => pp.f0,
            "DeltaTheta"     => pp.ΔΘ,
            "Theta0"         => pp.Θ0,
            "g0"             => pp.g0,
            "A_h"            => pp.A_h,
            "gamma_Theta"    => pp.γ_Θ,
            "E0"             => pp.E0,
            "dlngammaddelta" => pp.dlnγdδ,
            "gamma0"         => pp.γ0,
            "s0"             => pp.s0,
            "h_0"            => pp.h_0,
        )
    )
    
    
    ds = Dataset(;data..., properties=props)

    println("Saving model into file: ", filename)
    savedataset(ds,path = filename, driver=:netcdf, overwrite=overwrite)

end
