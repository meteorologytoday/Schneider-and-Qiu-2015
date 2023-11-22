
function savedata(m :: Model, filename :: String; overwrite=false)

    st = m.st
    gd = m.ev.gd

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
        :u_0   => YAXArray(axlist_cT,  st.u_0),
        :v_0   => YAXArray(axlist_cT,  st.v_0),
        :u_1   => YAXArray(axlist_U,   st.u_1),
        :v_1   => YAXArray(axlist_V,   st.v_1),
        :w_1   => YAXArray(axlist_W,   st.w_1),
        :h_1   => YAXArray(axlist_sT,  st.h_1),
        :pt_1  => YAXArray(axlist_sT,  st.Θ_1),
        :sst_1 => YAXArray(axlist_sT,  st.SST_1),
    )

    ds = Dataset(;data...)

    println("Saving model into file: ", filename)
    savedataset(ds,path = filename, driver=:netcdf, overwrite=overwrite)

end
