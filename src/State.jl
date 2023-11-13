mutable struct State

    X    :: AbstractArray{Float64, 2}
    _X   :: AbstractArray{Float64, 1}

    Γ   :: AbstractArray{Float64, 1}
    B   :: AbstractArray{Float64, 1}

    Ψ   :: AbstractArray{Float64, 1}
    SST :: AbstractArray{Float64, 1}
    
    Ψ_wgt :: AbstractArray{Float64, 1}

    function State(
        ev :: Env,
    )
        Ny = ev.gd.Ny
 
        X   = zeros( Float64, Ny, 2)
        _X  = view(X, :)
        
        Γ   = view(X, :, 1)
        B   = view(X, :, 2)

        Ψ   = zeros( Float64, Ny+1) 
        SST = zeros( Float64, Ny ) 
        
        Ψ_wgt = zeros( Float64, Ny+1) 

        return new(

            X,
            _X,

            Γ,
            B,

            Ψ,
            SST,
            
            Ψ_wgt,

        )    
    end
end

