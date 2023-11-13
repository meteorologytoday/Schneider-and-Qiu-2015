mutable struct Core

    amo :: AdvancedMatrixOperators
    amo_col :: AdvancedMatrixOperators
    amo_slb :: AdvancedMatrixOperators
    ops :: Dict

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd)
        amo_col = AdvancedMatrixOperators(gd=ev.gd_col)
        amo_slb = AdvancedMatrixOperators(gd=ev.gd_slb)
        ops = Dict(
        )

        return new(
            amo,
            amo_col,
            amo_slb,
            ops,
        )    
    end
end
