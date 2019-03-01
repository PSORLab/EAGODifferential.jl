# Plot ideas convex/concave relaxations of 1D function versus time
# Relaxations of objective function in 2D for 2D function
# Kinetic Data Best Fit
# Convergence Plot (Semi-log & Linear)

#=
function save_history(name, var, opt, annot)

    df = DataFrame()
    history = opt.history

    df.get_lower_bound = get_lower_bound(history)
    df.get_upper_bound = get_upper_bound(history)
    df.get_lower_time = get_lower_time(history)
    df.get_upper_time = get_upper_time(history)
    df.get_preprocess_time = get_preprocess_time(history)
    df.get_postprocess_time = get_postprocess_time(history)
    df.get_solution_time = get_solution_time(history)
    df.get_iteration_number = get_iteration_number(history)
    df.get_node_count = get_node_count(history)
    df.get_absolute_gap = get_absolute_gap(history)
    df.get_relative_gap = get_relative_gap(history)

    df_sol = DataFrame()
    df_sol.solution_value = MOI.get.(opt, MOI.VariablePrimal(), var)

    solution_path = export_path*name*"_solution.csv"
    save_path_data = export_path*name*".csv"
    ann_path_data = export_path*name*".txt"

    # write csv with history
    CSV.write(save_path_data, df)

    # write csv with solution
    CSV.write(solution_path, df_sol)

    # write annotation
    annotate_file = open(ann_path_data, create = true)
    write(annotate_file, annot)
    close(annotate_file)
end

const grid_size = 20

function save_bounds(name, ann, f, h, np, nx, nt, s, t_start, t_end, method,
                     pL, pU, xL, xU, x0; g = g, hj = hj)

    @assert np < 4

    # creates storage points
    mc_lower_bounds = zeros(nx,(nt-1))
    mc_upper_bounds = zeros(nx,(nt-1))
    pi_mc_lower_bounds = zeros(nx,(nt-1))
    pi_mc_upper_bounds = zeros(nx,(nt-1))

    string_mc_lower_bounds = "mc_lower_bounds"
    string_mc_upper_bounds = "mc_upper_bounds"
    string_pi_mc_lower_bounds = "pi_mc_lower_bounds"
    string_pi_mc_upper_bounds = "pi_mc_upper_bounds"
    for i in 1:nx
        df[Symbol(string_mc_lower_bounds*"_$i")]  = -Inf*ones(nt-1)
        df[Symbol(string_mc_upper_bounds*"_$i")]  = Inf*ones(nt-1)
        df[Symbol(string_pi_mc_lower_bounds*"_$i")]  =-Inf*ones(nt-1)
        df[Symbol(string_pi_mc_upper_bounds*"_$i")]  = Inf*ones(nt-1)
    end

    # creates and build the lower evaluator
    lower_evaluator_mc = ImplicitODELowerEvaluator{np}()
    lower_evaluator_mc_pi = ImplicitODELowerEvaluator{np}()
    upper_evaluator = ImplicitODEUpperEvaluator()
    EAGO.build_evaluator!(lower_evaluator_mc, f, h, np, nx, nt, s, t_start, t_end,
                          method, pL, pU, xL, xU, x0; g = g, hj = hj)
    EAGO.build_evaluator!(lower_evaluator_mc_pi, f, h, np, nx, nt, s, t_start, t_end,
                          method, pL, pU, xL, xU, x0; g = g, hj = hj)
    EAGO.build_evaluator!(upper_evaluator, f, h, np, nx, nt, s, t_start, t_end,
                          method, pL, pU, xL, xU, x0; g = g, hj = hj)

    # Creates the initial node and set it for each evaluator
    lower_vars = Float64[]
    upper_vars = Float64[]
    for i in 1:nt-1
        append!(lower_vars, xL)
        append!(upper_vars, xU)
    end
    append!(lower_vars, pL)
    append!(upper_vars, pU)
    n_mc = EAGO.NodeBB(lower_vars, upper_vars, -3.4, 2.1, 2, 1, true)

    # perform midpoint evaluation to generate reference point (McCormick Only)
    EAGO.set_current_node!(lower_eval_mc, n_mc)
    y_mid_mc = (lower_vars + upper_vars)/2.0
    EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc, y_mid_mc)

    # runs upper evaluator to contract box
    n_mc = EAGO.NodeBB(lower_vars, upper_vars, -3.4, 2.1, 2, 1, true)
    EAGO.set_current_node!(upper_evaluator, n_mc)

    # gets post-contractor node
    lower_vars = Float64[]
    upper_vars = Float64[]
    if (nx == 1)
        append!(lower_vars, lo(upper_evaluator.state_relax_1))
        append!(upper_vars, hi(upper_evaluator.state_relax_1))
    else
        for i in 1:nt-1
            append!(lower_vars, lo(upper_evaluator.state_relax_n[:,i]))
            append!(upper_vars, hi(upper_evaluator.state_relax_n[:,i]))
        end
    end
    append!(lower_vars, lo(upper_evaluator.P))
    append!(upper_vars, hi(upper_evaluator.P))
    n_mc_pi = EAGO.NodeBB(lower_vars, upper_vars, -3.4, 2.1, 2, 1, true)

    # performs midpoint evaluation to genrate reference point (McCormick + PI)
    EAGO.set_current_node!(lower_evaluator_mc_pi, n_mc_pi)
    y_mid_mc_pi = (lower_vars + upper_vars)/2.0
    EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc, y_mid_mc_pi)

    if np == 1
        for (i,px) in range(pL[1], pU[1], length = grid_size)
            EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc, y)
            for 1:nx
                df[Symbol(string_mc_lower_bounds*"_$i")]  = -Inf*ones(nt-1)
                df[Symbol(string_mc_upper_bounds*"_$i")]  = Inf*ones(nt-1)
                df[Symbol(string_pi_mc_lower_bounds*"_$i")]  =-Inf*ones(nt-1)
                df[Symbol(string_pi_mc_upper_bounds*"_$i")]  = Inf*ones(nt-1)
            end
            EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc_pi, y)
            for 1:nx
            end
        end
    elseif n == 2
        #=
        for (i,px) in range(pL[1], pU[1], length = grid_size)
            for (j,py) in range(pL[2], pU[2], length = grid_size)
            end
        end
        =#
    elseif
        #=
        for (i,px) in range(pL[1], pU[1], length = grid_size)
            for (j,py) in range(pL[2], pU[2], length = grid_size)
                for (k,pz) in range(pL[3], pU[3], length = grid_size)
                end
            end
        end
        =#
    end

    # create initial condition data frame
    p_mc = MC{np}.((pL+pU)/2.0, EAGO.IntervalType.(pL,pU), 1:np)
    x0_mc = x0(p_mc)

    # add IC to existing data frames

    # saves trajectories

end
=#
