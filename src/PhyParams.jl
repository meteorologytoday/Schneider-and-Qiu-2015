mutable struct PhyParams

    J :: Float64   # momentum vertical diffusivity 
    K :: Float64   # buoyancy horizontal diffusivity
    N :: Float64   # Stability of the atmosphere 
    C :: Float64   
    k :: Float64   
    
    L :: AbstractArray{Float64, 1}   # momentum horizontal diffusivity

end
