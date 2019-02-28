@testset "Build Upper Evaluator, Set Node" begin
    upper_eval = ImplicitODEUpperEvaluator()

    f(x,p,t) = x[1]
    function h(H,x,p,t)
        #println("(p[1]*x[1]): $(p[1]*x[1])")
        H[1] = -p[1]*x[1]
    end
    function hj(J,x,p,t)
        J[1,1] = -p[1]
    end

    np = 1
    nx = 1
    nt = 100
    s = 2

    t_start = 0.0
    t_end = 1.0
    x0(p) = [1.0]
    method = :BDF

    pL = [-20.0]; pU = [-10.0]
    xL = [-0.12]; xU = [-0.04]

    # build the basic evaluator (w/o inequality constraints)
    build_evaluator!(upper_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU,  x0; hj = hj)

    # update the node info for the evaluator to initial node
    n1 = EAGO.NodeBB(Float64[-0.12,-20.0], Float64[-0.04,-10.0], -3.4, 2.1, 2, 1, true)
    EAGO.set_current_node!(upper_eval, n1)

    @test upper_eval.current_node.lower_variable_bounds[1] == -0.12
    @test upper_eval.current_node.lower_variable_bounds[2] == -20.0
    @test upper_eval.current_node.upper_variable_bounds[1] == -0.04
    @test upper_eval.current_node.upper_variable_bounds[2] == -10.0
    @test upper_eval.current_node.lower_bound == -3.4
    @test upper_eval.current_node.upper_bound == 2.1
    @test upper_eval.current_node.depth == 2
    @test upper_eval.current_node.last_branch == 1
    @test upper_eval.current_node.branch_direction

    n2 = EAGO.NodeBB(Float64[-0.121,-20.01], Float64[-0.041,-10.01], -3.41, 2.11, 21, 11, false)
    EAGO.set_last_node!(upper_eval, n2)

    @test upper_eval.last_node.lower_variable_bounds[1] == -0.121
    @test upper_eval.last_node.lower_variable_bounds[2] == -20.01
    @test upper_eval.last_node.upper_variable_bounds[1] == -0.041
    @test upper_eval.last_node.upper_variable_bounds[2] == -10.01
    @test upper_eval.last_node.lower_bound == -3.41
    @test upper_eval.last_node.upper_bound == 2.11
    @test upper_eval.last_node.depth == 21
    @test upper_eval.last_node.last_branch == 11
    @test ~upper_eval.last_node.branch_direction

    @test upper_eval.ivp.time_steps == 100
    @test upper_eval.ivp.step_size == 1.0/(99)
    @test upper_eval.ivp.time_start == 0.0
    @test upper_eval.ivp.time_end == 1.0

    @test upper_eval.has_ineq == false

    # build evaluator w/ equality constaints
    g(x,p,t) = [x[1][1] + p[1]; p[1]+t[1]];
    build_evaluator!(upper_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU,  x0; g = g, hj = hj)
    @test upper_eval.has_ineq == true
    @test upper_eval.ng == 2

    gval = upper_eval.constraints_fun([[2.1]], [3.4], [1.1])
    @test gval[1] == 5.5
    @test gval[2] == 4.5
end

@testset "Upper Evaluator - Relax Objective/Constraint" begin
    # soft build, then evaluate
    upper_eval = ImplicitODEUpperEvaluator()

    f(x,p,t) = x[1][1]
    g(x,p,t) = [3.0*x[1][1]]

    upper_eval.np = 2
    upper_eval.var_relax = [EAGO.IntervalType(1.0, 2.4);  EAGO.IntervalType(1.0, 2.4)]
    upper_eval.state_relax_n = fill(EAGO.IntervalType(1.1, 2.1),(2,2,))

    upper_eval.has_ineq = true
    upper_eval.objective_fun = f
    upper_eval.constraints_fun = g
    upper_eval.obj_relax = EAGO.IntervalType(1.0, 2.4)
    upper_eval.constraint_relax = [EAGO.IntervalType(1.0, 2.4)]

    upper_eval.obj_eval = true
    EAGO_Differential.relax_objective!(upper_eval)
    @test upper_eval.obj_relax.lo == 1.0
    @test upper_eval.obj_relax.hi == 2.4
    @test upper_eval.obj_eval == true

    upper_eval.obj_eval = false
    EAGO_Differential.relax_objective!(upper_eval)
    @test upper_eval.obj_relax.lo == 1.1
    @test upper_eval.obj_relax.hi == 2.1
    @test upper_eval.obj_eval == true

    upper_eval.cnstr_eval = true
    EAGO_Differential.relax_constraints!(upper_eval)
    @test upper_eval.constraint_relax[1].lo == 1.0
    @test upper_eval.constraint_relax[1].hi == 2.4
    @test upper_eval.cnstr_eval == true

    upper_eval.cnstr_eval = false
    EAGO_Differential.relax_constraints!(upper_eval)
    @test isapprox(upper_eval.constraint_relax[1].lo, 3.3, atol = 1E-6)
    @test isapprox(upper_eval.constraint_relax[1].hi, 6.3, atol = 1E-6)
    @test upper_eval.cnstr_eval == true
end

@testset "Upper Relaxation Calculation Routines (1D - Point)" begin
    # soft build, then evaluate
    upper = ImplicitODEUpperEvaluator()

    f(x,p,t) = x[1]
    function h(H,x,p,t)
        H[1] = p[1]*x[1]*(1-x[1])
    end
    function hj(J,x,p,t)
        J[1,1] = p[1]*(1.0-2.0*x[1])
    end
    x0(p) = [0.1]

    np = 1
    nx = 1
    nt = 10
    s = 2

    t_start = 0.0
    t_end = 1.0
    method = :BDF

    pL = [1.5]; pU = [1.5]
    xL = [0.1]; xU = [1.0]

    # build the basic evaluator (w/o inequality constraints)
    EAGO_Differential.build_evaluator!(upper, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)

    lower_vars = fill(xL[1], (nt-1,))
    upper_vars = fill(xU[1], (nt-1,))
    append!(lower_vars, pL)
    append!(upper_vars, pU)
    n = NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)

    y = [1.5]
    EAGO.set_current_node!(upper, n)
    EAGO_Differential.relax_ode_implicit!(upper, y)

    y1 = [1.55]
    EAGO_Differential.relax_ode_implicit!(upper, y1)

    @test isapprox(upper.state_relax_1[1].lo, 0.11725, atol=1E-5)
    @test isapprox(upper.state_relax_1[3].lo, 0.159266, atol=1E-5)
    @test isapprox(upper.state_relax_1[6].lo, 0.242828, atol=1E-5)
    @test isapprox(upper.state_relax_1[9].lo, 0.349821, atol=1E-5)
    @test isapprox(upper.state_relax_1[1].hi, 0.117251, atol=1E-5)
    @test isapprox(upper.state_relax_1[3].hi, 0.159267, atol=1E-5)
    @test isapprox(upper.state_relax_1[6].hi, 0.242829, atol=1E-5)
    @test isapprox(upper.state_relax_1[9].hi, 0.349822, atol=1E-5)
end

#=
@testset "Upper Relaxation Calculation Routines (3D)" begin
end
=#

@testset "Upper Evaluator MOI Wrapper" begin

    #=
        MOI.eval_objective(d::ImplicitODEUpperEvaluator, y)
        MOI.eval_constraint(d::ImplicitODEUpperEvaluator, g, y)
        MOI.eval_objective_gradient(d::ImplicitODEUpperEvaluator, df, y)
        MOI.jacobian_structure(d::ImplicitODEUpperEvaluator)
        MOI.hessian_lagrangian_structure(d::ImplicitODEUpperEvaluator)
        _hessian_lagrangian_structure(d::ImplicitODEUpperEvaluator)
        MOI.eval_constraint_jacobian(d::ImplicitODEUpperEvaluator,g,y)
        MOI.eval_constraint_jacobian_product(d::ImplicitODEUpperEvaluator, out, y, w)
        MOI.eval_constraint_jacobian_transpose_product(d::ImplicitODEUpperEvaluator, out, y, w)
        MOI.features_available(d::ImplicitODEUpperEvaluator)
        MOI.initialize(d::ImplicitODEUpperEvaluator, requested_features::Vector{Symbol}) end
        MOI.objective_expr(d::ImplicitODEUpperEvaluator)
        MOI.constraint_expr(d::ImplicitODEUpperEvaluator)
    =#
end
