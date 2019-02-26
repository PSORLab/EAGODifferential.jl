@testset "Build Lower Evaluator, Set Node" begin
    lower_eval = ImplicitODELowerEvaluator{1}()

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
    x0 = [1.0]
    method = :BDF

    pL = [-20.0]; pU = [-10.0]
    xL = [-0.12]; xU = [-0.04]

    # build the basic evaluator (w/o inequality constraints)
    build_evaluator!(lower_eval, f, h, np, nx, nt, s, t_start, t_end, x0, method, pL, pU, xL, xU; hj = hj)

    # update the node info for the evaluator to initial node
    n1 = EAGO.NodeBB(Float64[-0.12,-20.0], Float64[-0.04,-10.0], -3.4, 2.1, 2, 1, true)
    EAGO.set_current_node!(lower_eval, n1)

    @test lower_eval.current_node.lower_variable_bounds[1] == -0.12
    @test lower_eval.current_node.lower_variable_bounds[2] == -20.0
    @test lower_eval.current_node.upper_variable_bounds[1] == -0.04
    @test lower_eval.current_node.upper_variable_bounds[2] == -10.0
    @test lower_eval.current_node.lower_bound == -3.4
    @test lower_eval.current_node.upper_bound == 2.1
    @test lower_eval.current_node.depth == 2
    @test lower_eval.current_node.last_branch == 1
    @test lower_eval.current_node.branch_direction

    n2 = EAGO.NodeBB(Float64[-0.121,-20.01], Float64[-0.041,-10.01], -3.41, 2.11, 21, 11, false)
    EAGO.set_last_node!(lower_eval, n2)

    @test lower_eval.last_node.lower_variable_bounds[1] == -0.121
    @test lower_eval.last_node.lower_variable_bounds[2] == -20.01
    @test lower_eval.last_node.upper_variable_bounds[1] == -0.041
    @test lower_eval.last_node.upper_variable_bounds[2] == -10.01
    @test lower_eval.last_node.lower_bound == -3.41
    @test lower_eval.last_node.upper_bound == 2.11
    @test lower_eval.last_node.depth == 21
    @test lower_eval.last_node.last_branch == 11
    @test ~lower_eval.last_node.branch_direction

    @test lower_eval.ivp.time_steps == 100
    @test lower_eval.ivp.step_size == 1.0/(99)
    @test lower_eval.ivp.time_start == 0.0
    @test lower_eval.ivp.time_end == 1.0
    @test lower_eval.ivp.initial_values[1] == 1.0

    @test lower_eval.IC_relaxations[1].cc == 1.0
    @test lower_eval.IC_relaxations[1].cv == 1.0
    @test lower_eval.IC_relaxations[1].Intv.lo == 1.0
    @test lower_eval.IC_relaxations[1].Intv.hi == 1.0

    @test lower_eval.has_ineq == false

    # build evaluator w/ equality constaints
    g(x,p,t) = [x[1][1] + p[1]; p[1]+t[1]];
    build_evaluator!(lower_eval, f, h, np, nx, nt, s, t_start, t_end, x0, method, pL, pU, xL, xU; g = g, hj = hj)
    @test lower_eval.has_ineq == true
    @test lower_eval.ng == 2

    gval = lower_eval.constraints_fun([[2.1]], [3.4], [1.1])
    @test gval[1] == 5.5
    @test gval[2] == 4.5
end

@testset "Lower Evaluator - Relax Objective/Constraint" begin

    # soft build, then evaluate
    lower_eval = ImplicitODELowerEvaluator{2}()

    f(x,p,t) = x[1][1]
    g(x,p,t) = [3.0*x[1][1]]

    lower_eval.np = 2
    lower_eval.var_relax = [MC{2}(1.0, 2.4);  MC{2}(1.0, 2.4)]
    lower_eval.state_relax_n = fill(MC{2}(1.1, 2.1),(2,2,))

    lower_eval.has_ineq = true
    lower_eval.objective_fun = f
    lower_eval.constraints_fun = g
    lower_eval.obj_relax = MC{2}(1.0, 2.4)
    lower_eval.constraint_relax = [MC{2}(1.0, 2.4)]

    lower_eval.obj_eval = true
    EAGO_Differential.relax_objective!(lower_eval)
    @test lower_eval.obj_relax.cv == 1.0
    @test lower_eval.obj_relax.cc == 2.4
    @test lower_eval.obj_relax.Intv.lo == 1.0
    @test lower_eval.obj_relax.Intv.hi == 2.4
    @test lower_eval.obj_eval == true

    lower_eval.obj_eval = false
    EAGO_Differential.relax_objective!(lower_eval)
    @test lower_eval.obj_relax.cv == 1.1
    @test lower_eval.obj_relax.cc == 2.1
    @test lower_eval.obj_relax.Intv.lo == 1.1
    @test lower_eval.obj_relax.Intv.hi == 2.1
    @test lower_eval.obj_eval == true

    lower_eval.cnstr_eval = true
    EAGO_Differential.relax_constraints!(lower_eval)
    @test lower_eval.constraint_relax[1].cv == 1.0
    @test lower_eval.constraint_relax[1].cc == 2.4
    @test lower_eval.constraint_relax[1].Intv.lo == 1.0
    @test lower_eval.constraint_relax[1].Intv.hi == 2.4
    @test lower_eval.cnstr_eval == true

    lower_eval.cnstr_eval = false
    EAGO_Differential.relax_constraints!(lower_eval)
    @test isapprox(lower_eval.constraint_relax[1].cv, 3.3, atol = 1E-6)
    @test isapprox(lower_eval.constraint_relax[1].cc, 6.3, atol = 1E-6)
    @test isapprox(lower_eval.constraint_relax[1].Intv.lo, 3.3, atol = 1E-6)
    @test isapprox(lower_eval.constraint_relax[1].Intv.hi, 6.3, atol = 1E-6)
    @test lower_eval.cnstr_eval == true
end

@testset "Lower Relaxation Calculation Routines (1D)" begin
    # soft build, then evaluate
    lower_eval = ImplicitODELowerEvaluator{1}()

    f(x,p,t) = x[1][1]
    h(H,x,p,t) = [-p[1]*x[1][1]]
    hj(J,x,p,t) = [-p[1]]

    np = 1
    nx = 1
    nt = 100
    s = 2

    t_start = 0.0
    t_end = 1.0
    x0 = [1.0]
    method = :BDF

    pL = [-20.0]; pU = [-10.0]
    xL = [0.00]; xU = [1.00]

    # build the basic evaluator (w/o inequality constraints)
    EAGO_Differential.build_evaluator!(lower_eval, f, h, np, nx, nt, s, t_start, t_end, x0, method, pL, pU, xL, xU; hj = hj)

    y = [-16.0]
    EAGO_Differential.relax_ode_implicit!(lower_eval, y)
end

#=
@testset "Lower Relaxation Calculation Routines (3D)" begin
    EAGO_Differential.relax_ode_implicit!(lower_eval, y)
end

@testset "Lower Evaluator MOI Wrapper" begin

    lower_eval = ImplicitODELowerEvaluator{2}()

    f_3(x,p,t) = 3.0*x[1]
    g_3(x,p,t) = [2.0*x[2]; 3.0*x[2]]

    #= TODO
    function h_3!(H,x,xp,p,t)
    end
    function hj_3!(J,x,xp,p,t)
    end
    =#

    np = 2
    nx = 2
    nt = 100
    s = 2
    t_start = 0.0
    t_end = 1.0
    x0 = [1.5, 2.2]
    method = :BDF

    #= TODO
    pL = [-20.0; -20.0]
    pU = [-10.0; -10.0]
    xL = [-0.12; -0.12]
    xU = [-0.04; -0.04]

    build_evaluator!(lower_eval, f_3, h_3!, np, nx, nt, s, t_start, t_end,
                     x0, method, pL, pU, xL, xU; g = g_3, hj = hj_3!)
    =#

    # update the node info for the evaluator to initial node
    n1 = EAGO.NodeBB(Float64[-0.12,-20.0], Float64[-0.04,-10.0], -3.4, 2.1, 2, 1, true)
    EAGO.set_current_node!(lower_eval, n1)

    #y = [] TODO
    #w = [] TODO

    fval = MOI.eval_objective(lower_eval, y)
    #@test fval == 1 # TODO

    g = zeros[2]
    MOI.eval_constraint(lower_eval, g, y)
    #@test g[1] == 0.0 # TODO
    #@test g[2] == 0.0 # TODO

    df = zeros[2]
    MOI.eval_objective_gradient(lower_eval, df, y)
    #@test df[1] == 0.0 # TODO
    #@test df[2] == 0.0 # TODO

    jacobian_structure = MOI.jacobian_structure(lower_eval)
    @test jacobian_structure[1][1] == 1
    @test jacobian_structure[1][2] == 1
    @test jacobian_structure[2][1] == 1
    @test jacobian_structure[2][2] == 2
    @test jacobian_structure[3][1] == 2
    @test jacobian_structure[3][2] == 1
    @test jacobian_structure[4][1] == 2
    @test jacobian_structure[4][2] == 2

    @test_throws MOI.hessian_lagrangian_structure(lower_eval)
    @test_throws EAGO_Differential._hessian_lagrangian_structure(lower_eval)

    g = zeros[2,2]
    MOI.eval_constraint_jacobian(lower_eval, g, y)
    #@test g[1,1] == 0.0 TODO
    #@test g[1,2] == 0.0 TODO
    #@test g[2,1] == 0.0 TODO
    #@test g[2,2] == 0.0 TODO

    out = zeros[2]
    MOI.eval_constraint_jacobian_product(lower_eval, out, y, w)
    #@test out[1] == 0.0 TODO
    #@test out[2] == 0.0 TODO

    out = zeros[2]
    MOI.eval_constraint_jacobian_transpose_product(lower_eval, out, y, w)
    #@test out[1] == 0.0 TODO
    #@test out[2] == 0.0 TODO

    features = MOI.features_available(lower_eval)
    @test features[1] == :Grad
    @test features[2] == :Jac

    requested_features = Symbol[]
    @test_nowarn MOI.initialize(lower_eval, requested_features)

    @test_throws ErrorException MOI.objective_expr(lower_eval)
    @test_throws ErrorException MOI.constraint_expr(lower_eval)
end
=#
