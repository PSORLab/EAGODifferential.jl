function interval_preprocess_ode!(x::EAGO.Optimizer, y::EAGO.NodeBB)
    Eflag = false
    Iflag = false
    eDflag = false
    evaluator = x.nlp_data.evaluator
    nx = evaluator.nx
    np = evaluator.np
    nt = evaluator.ivp.time_steps

    evaluator.obj_eval = false
    evaluator.cnstr_eval = false
    evaluator.init_relax_run = false

    set_current_node!(evaluator, y)
    yval = (y.lower_variable_bounds + y.upper_variable_bounds)/2.0
    relax_ode_implicit!(evaluator)

    g = zeros(evaluator.ng)
    MOI.eval_constraint(evaluator, g, yval)
    Eflag = evaluator.exclusion_flag
    if ~Eflag
        EFlag = any(i-> (i > 0), g)
    end

    x.current_preprocess_info.feasibility = ~Eflag
    if ~Eflag
        if nx == 1
            y.lower_variable_bounds[1:(nt-1)] = lo.(evaluator.state_relax_1)
            y.upper_variable_bounds[1:(nt-1)] = hi.(evaluator.state_relax_1)
        else
            for i in 1:(nt-1)
                y.lower_variable_bounds[(1+(i-1)*nx):(i*nx)] = lo.(evaluator.state_relax_n[:,i])
                y.upper_variable_bounds[(1+(i-1)*nx):(i*nx)] = hi.(evaluator.state_relax_n[:,i])
            end
        end
    end
end

function create_mid_node(y::NodeBB, nx::Int, np::Int, nt::Int)

    lower_variable_bounds = y.lower_variable_bounds
    upper_variable_bounds = y.upper_variable_bounds

    P_interval = EAGO.IntervalType.(lower_variable_bounds[(nx*(nt-1)+1):(nx*(nt-1)+np)],
                                    upper_variable_bounds[(nx*(nt-1)+1):(nx*(nt-1)+np)])
    P_mid_interval = EAGO.IntervalType.(mid.(P_interval))

    ymid = deepcopy(y)

    ymid.lower_variable_bounds[(nx*(nt-1)+1):(nx*(nt-1)+np)] = lo.(P_mid_interval[:])
    ymid.upper_variable_bounds[(nx*(nt-1)+1):(nx*(nt-1)+np)] = hi.(P_mid_interval[:])

    return ymid
end

function midpoint_upper_bnd_ode!(x::EAGO.Optimizer, y::NodeBB)
    #println("START UPPER EVALUATOR")
    if EAGO.is_integer_feasible(x) #&& mod(x.CurrentIterationCount,x.UpperBoundingInterval) == 1

        evaluator = x.nlp_data.evaluator
        #println("evaluator.ng: $(evaluator.ng)")
        nx = evaluator.nx
        np = evaluator.np
        nt = evaluator.ivp.time_steps
        evaluator.obj_eval = false
        evaluator.cnstr_eval = false
        evaluator.init_relax_run = false

        node_ymid = create_mid_node(y, nx, np, nt)
        set_current_node!(evaluator, node_ymid)
        EAGO.set_to_mid!(x.current_upper_info.solution, y)
        Eflag = evaluator.exclusion_flag
        if evaluator.ng > 0
            g = zeros(evaluator.ng)
            MOI.eval_constraint(evaluator, g, x.current_upper_info.solution)
            #println("g: $(g)")
            g[1] = mid(evaluator.state_relax_n[20,29]) - 0.08
            result_status = any(i-> (i > 0), g) ? MOI.INFEASIBLE_POINT : MOI.FEASIBLE_POINT
            #println("result_status: $(result_status)")
            if (result_status == MOI.FEASIBLE_POINT)
                result_status = evaluator.exclusion_flag ? MOI.INFEASIBLE_POINT : MOI.FEASIBLE_POINT
            end
        else
            result_status = MOI.FEASIBLE_POINT
        end
        if (result_status == MOI.FEASIBLE_POINT)
            #println("result_status == MOI.FEASIBLE_POINT")
            #println("started eval objective")
            x.current_upper_info.feasibility = true
            #println("x.current_upper_info.solution: $(x.current_upper_info.solution)")
            val = MOI.eval_objective(evaluator, x.current_upper_info.solution)
            #println("evaluator: $(evaluator.state_relax_n[20,:])")
            #println("val: $(val)")
            x.current_upper_info.value = val
            x.current_upper_info.value = mid(evaluator.P[1])
        else
            #println("result_status == MOI.INFEASIBLE_POINT")
            x.current_upper_info.feasibility = false
            x.current_upper_info.value = Inf
        end
    else
        x.current_upper_info.feasibility = false
        x.current_upper_info.value = Inf
    end
end

function concave_upper_bnd_ode!(x::EAGO.Optimizer, y::NodeBB)
    # Copies initial model into working model (if initial model isn't dummy)
    # A dummy model is left iff all terms are relaxed
    if x.use_lower_factory
        factory = x.lower_factory(;x.lower_optimizer_options...)     # Should accept keyword arguments
        x.working_relaxed_optimizer = factory
        MOI.add_variables(x.working_relaxed_optimizer, x.variable_number)
    else
        if x.initial_relaxed_optimizer != EAGO.DummyOptimizer()
            x.working_relaxed_optimizer = deepcopy(x.initial_relaxed_optimizer)
        end
    end

    EAGO.update_lower_variable_bounds1!(x,y,x.working_relaxed_optimizer)
    ngrad = x.variable_number
    nx =  x.state_variables
    np = ngrad - nx
    var = x.upper_variables

    x.working_evaluator_block.evaluator.current_node = y
    midx = (y.upper_variable_bounds + y.lower_variable_bounds)/2.0
    grad_c = 0.0

    # Add objective
    if x.working_evaluator_block.has_objective
        # Calculates relaxation and subgradient
        df = zeros(Float64, np)
        f = MOI.eval_objective(x.working_evaluator_block.evaluator, midx)
        MOI.eval_objective_gradient(x.working_evaluator_block.evaluator, df, midx)
        fcc = x.working_evaluator_block.evaluator.obj_relax.cc # FIX ME
        dfcc = x.working_evaluator_block.evaluator.obj_relax.cc_grad # FIX ME

        # Add objective relaxation to model
        saf_const = fcc
        for i in (1+nx):ngrad
            grad_c = dfcc[i-nx]
            midx_c = midx[i]
            saf_const -= midx_c*grad_c
        end
        saf = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(fcc, var[(1+nx):ngrad]), saf_const)
        MOI.set(x.working_relaxed_optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), -saf)
    end

    # Optimizes the object
    tt = stdout
    redirect_stdout()
    MOI.optimize!(x.working_relaxed_optimizer)
    redirect_stdout(tt)


    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
    result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = EAGO.is_globally_optimal(termination_status, result_status_code)
    solution = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), x.lower_variables)

    # specifies node used in last problem
    x.working_evaluator_block.evaluator.last_node = y
    x.working_evaluator_block.evaluator.objective_ubd = x.global_upper_bound

    if valid_flag
        if feasible_flag
            x.current_upper_info.feasibility = true
            x.current_upper_info.value = -MOI.get(x.working_relaxed_optimizer, MOI.ObjectiveValue())
            x.current_upper_info.solution[1:end] = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), x.lower_variables)
            EAGO.set_dual!(x)
        else
            x.current_upper_info.feasibility = false
            x.current_upper_info.value = Inf
        end
    else
        x.current_upper_info.feasibility = true
        x.current_upper_info.value = Inf
    end
end

# Modifies functions post initial relaxation to use appropriate nlp evaluator
function ode_mod!(opt::Optimizer, args)

    # unpack args
    ImpLowerEval = args[1]
    ImpUpperEval = args[2]
    lower_bnds = args[3]
    upper_bnds = args[4]

    # load ode custom functions
    opt.preprocess! = interval_preprocess_ode!
    opt.relax_function! = EAGO.implicit_relax_model!
    #opt.upper_problem! = concave_upper_bnd_ode!
    opt.upper_problem! = midpoint_upper_bnd_ode!

    # load lower nlp data block
    lower_eval_block = MOI.NLPBlockData(lower_bnds, ImpLowerEval, true)
    opt.working_evaluator_block = deepcopy(lower_eval_block)
    if MOI.supports(opt.initial_relaxed_optimizer, MOI.NLPBlock())
        opt.initial_relaxed_optimizer.nlp_data = deepcopy(lower_eval_block)
        opt.working_relaxed_optimizer.nlp_data = deepcopy(lower_eval_block)
    end

    # load upper nlp data block
    upper_eval_block = MOI.NLPBlockData(upper_bnds, ImpUpperEval, true)
    opt.nlp_data = upper_eval_block

    # specifies that the decision variables participate in nonlinear expressions
    # and should be branched on
    num_state = num_state_variables(ImpLowerEval)
    num_decis = num_decision_variables(ImpLowerEval)
    for i in (num_state+1):(num_state+num_decis)
        opt.nonlinear_variable[i] = true
    end
end

"""
    solve_ode

Solves the optimization problem `min_{p} f(x,x0,p,t)` with respect to_indices
`dx/dt = h(x,p,t)` and `g(x,x0,p,t) <= 0` on `x in X` and `p in P` such that x_0 = x_0(p).
- f(x,x0,p,t) takes the full state vector, the decision vector, and the time vector.
The time vector may be useful for computation but the problem is solved in p only.
So for x[3,4] is x_3 at the after 4 time steps (or at time point 5). So x0 is the
initial condition.
"""
function solve_ode(f, h, hj, g, x0, xL, xU, pL, pU, t_start, t_end, nt, s, method, opt; ng = 0,
                   state_update = x -> (), user_time = [], user_xtL = [], user_xtU = [])

    # get dimensions & check for consistency
    @assert length(pL) == length(pU)
    @assert length(xL) == length(xU)
    @assert length(user_xtL) == length(user_xtU)
    np = length(pL); nx = length(xL);

    if (g == nothing)
        ng = 0
    elseif ng == 0
        ng = length(g(ones(nx,nt-1),ones(nx),ones(np),ones(nt)))
    end

    # sets most routines to default (expect bisection)
    EAGO.set_to_default!(opt)
    opt.bisection_function = EAGO.implicit_bisection

    # Add variables
    var_EAGO = MOI.add_variables(opt, nx*(nt-1)+np)
    if length(user_xtL) > 0
        for i in 1:length(user_xtL)
            MOI.add_constraint(opt, var_EAGO[i], MOI.GreaterThan(user_xtL[i]))
            MOI.add_constraint(opt, var_EAGO[i], MOI.LessThan(user_xtU[i]))
        end
    else
        count = 1
        for i in 1:(nt-1)
            for j in 1:nx
                MOI.add_constraint(opt, var_EAGO[count], MOI.GreaterThan(xL[j]))
                MOI.add_constraint(opt, var_EAGO[count], MOI.LessThan(xU[j]))
                count += 1
            end
        end
    end

    for i in 1:np
        MOI.add_constraint(opt, var_EAGO[i+nx*(nt-1)], MOI.GreaterThan(pL[i]))
        MOI.add_constraint(opt, var_EAGO[i+nx*(nt-1)], MOI.LessThan(pU[i]))
    end

    # Specify number of state variables
    opt.state_variables = nx*(nt-1)
    opt.upper_has_node = true

    # creates the appropriate lower evaluator
    lower = ImplicitODELowerEvaluator{np}()
    EAGO_Differential.build_evaluator!(lower, f, h, np, nx, nt, s, t_start, t_end,
                                       method, pL, pU, xL, xU, x0; hj = hj, g = g,
                                       state_update = state_update, user_time = user_time,
                                       user_xtL = user_xtL, user_xtU = user_xtU)
    upper = ImplicitODEUpperEvaluator()
    EAGO_Differential.build_evaluator!(upper, f, h, np, nx, nt, s, t_start,
                                       t_end, method, pL, pU, xL, xU,
                                       x0; hj = hj, g = g, user_time = user_time,
                                       user_xtL = user_xtL, user_xtU = user_xtU)

    # Add nlp data blocks ("SHOULD" BE THE LAST THING TO DO)
    bnd_pair = MOI.NLPBoundsPair(-Inf,0.0)
    lower_bnds = [bnd_pair for i=1:ng]
    upper_bnds = [bnd_pair for i=1:ng]

    # Set the objective sense
    MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    custom_mod_args = (lower, upper, lower_bnds, upper_bnds)
    MOI.optimize!(opt, custom_mod! = ode_mod!, custom_mod_args = custom_mod_args)
    return var_EAGO, opt
end
