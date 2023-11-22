

module BLM

    using SparseArrays, LinearAlgebra, Arpack
    using YAXArrays, DimensionalData, NetCDF

    function printInfo(x)
        for (i, sym) in enumerate(fieldnames(typeof(x)))
            println("[$i] $(string(sym)) = ", getfield(x, sym))
        end
    end

    include("constants.jl")
    include("tools.jl")


    include("Grid.jl")
    include("BasicMatrixOperators.jl")
    include("AdvancedMatrixOperators.jl")

    include("PhyParams.jl")
    include("Env.jl")
    include("Core.jl")
    include("State.jl")
    include("Model.jl")
    
    include("solver_functions.jl")
    
    include("io.jl")

#    include("dynamics.jl")
end
