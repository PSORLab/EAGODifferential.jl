module EAGO_Differential

    using EAGO, MathOptInterface, LinearAlgebra

    const MOI = MathOptInterface

    import EAGO: build_evaluator!, set_current_node!, set_last_node!

    include("src/calc_ode.jl")

    export ImplicitODELowerEvaluator, ImplicitODEUpperEvaluator, build_evaluator!, set_current_node!, set_last_node!

    include("src/lower_evaluator/lower_evaluator.jl")
    include("src/upper_evaluator/upper_evaluator.jl")
    include("src/solve_ode.jl")

end # module
