@testset "Create MidPoint Node" begin

    np = 1
    nx = 1
    nt = 10

    pL = [1.0]; pU = [3.0]
    xL = [0.1]; xU = [1.0]

    lower_vars = fill(xL[1], (nt-1,))
    upper_vars = fill(xU[1], (nt-1,))
    append!(lower_vars, pL)
    append!(upper_vars, pU)
    n =  EAGO.NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)

    ymid = EAGO_Differential.create_mid_node(n, nx, np, nt)
    @test ymid.lower_variable_bounds[10] == 2.0
    @test ymid.upper_variable_bounds[10] == 2.0
end

@testset "Parametric Interval Preprocessing" begin

    # soft build, then evaluate
    upper_eval = ImplicitODEUpperEvaluator()

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

    pL = [1.49]; pU = [1.59]
    xL = [0.1]; xU = [1.0]

    # build the basic evaluator (w/o inequality constraints)
    EAGO_Differential.build_evaluator!(upper_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)

    lower_vars = fill(xL[1], (nt-1,))
    upper_vars = fill(xU[1], (nt-1,))
    append!(lower_vars, pL)
    append!(upper_vars, pU)
    n = EAGO.NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)

    x = EAGO.Optimizer()
    x.nlp_data = MOI.NLPBlockData(MOI.NLPBoundsPair[], upper_eval, true)
    EAGO_Differential.interval_preprocess_ode!(x, n)

    @test isapprox(n.lower_variable_bounds[1], 0.117112, atol = 1E-4)
    @test isapprox(n.upper_variable_bounds[1], 0.118447, atol = 1E-4)

    @test isapprox(n.lower_variable_bounds[4], 0.183487, atol = 1E-4)
    @test isapprox(n.upper_variable_bounds[4], 0.191221, atol = 1E-4)

    @test isapprox(n.lower_variable_bounds[9], 0.346909, atol = 1E-4)
    @test isapprox(n.upper_variable_bounds[9], 0.372437, atol = 1E-4)

    @test isapprox(n.lower_variable_bounds[10], 1.49, atol = 1E-4)
    @test isapprox(n.upper_variable_bounds[10], 1.59, atol = 1E-4)

    @test x.current_preprocess_info.feasibility
end

@testset "Interval MidPoint Upper Bound" begin

    # soft build, then evaluate
    upper_eval = ImplicitODEUpperEvaluator()

    f(x,x0,p,t) = x[9]
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

    pL = [1.49]; pU = [1.59]
    xL = [0.1]; xU = [1.0]

    lower_vars = fill(xL[1], (nt-1,))
    upper_vars = fill(xU[1], (nt-1,))
    append!(lower_vars, pL)
    append!(upper_vars, pU)
    n = EAGO.NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)

    # build the basic evaluator (w/o inequality constraints)
    EAGO_Differential.build_evaluator!(upper_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)

    lower_vars = fill(xL[1], (nt-1,))
    upper_vars = fill(xU[1], (nt-1,))
    append!(lower_vars, pL)
    append!(upper_vars, pU)
    n = EAGO.NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)

    x = EAGO.Optimizer()
    x.nlp_data = MOI.NLPBlockData(MOI.NLPBoundsPair[], upper_eval, true)
    x.current_upper_info.solution = fill(0.0, (10,))
    EAGO_Differential.midpoint_upper_bnd_ode!(x, n)

    @test x.current_upper_info.feasibility
    @test isapprox(x.current_upper_info.value, 0.35968613286903667, atol=1E-5)
end

# TODO: ACTUAL SOLVE ROUTINE
@testset "Global pODE Optimizer" begin

    opt = EAGO.Optimizer()

    f(x,x0,p,t) = x[1]
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

    pL = [1.49]; pU = [1.59]
    xL = [0.1]; xU = [1.0]

    solve_ode(f, h, hj, nothing, x0, xL, xU, pL, pU, t_start, t_end, nt, s, method, opt)
end
