# Plot ideas convex/concave relaxations of 1D function versus time
# Relaxations of objective function in 2D for 2D function
# Kinetic Data Best Fit
# Convergence Plot (Semi-log & Linear)


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
    #annotate_file = open(ann_path_data, create = true)
    #write(annotate_file, annot)
    #close(annotate_file)
end

const grid_size = 10

function save_kernel!(df, ic_opts, y, nx,
                     lower_evaluator_mc,
                     lower_evaluator_mc_pi,
                     upper_evaluator_pi,
                     string_mc_lower_bounds,
                     string_mc_upper_bounds,
                     string_pi_mc_lower_bounds,
                     string_pi_mc_upper_bounds,
                     string_pi_lower_bounds,
                     string_pi_upper_bounds)

    # updates bounds from mc (mc + intv) only
    EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc, y)
    EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc_pi, y)
    for i in 1:nx
        if (nx == 1)
            df[Symbol(string_mc_lower_bounds*"_$i")] = min.(df[Symbol(string_mc_lower_bounds*"_$i")], cv.(lower_evaluator_mc.state_relax_1))
            df[Symbol(string_mc_upper_bounds*"_$i")] = max.(df[Symbol(string_mc_upper_bounds*"_$i")], cc.(lower_evaluator_mc.state_relax_1))
            df[Symbol(string_pi_mc_lower_bounds*"_$i")] = min.(df[Symbol(string_pi_mc_lower_bounds*"_$i")], cv.(lower_evaluator_mc_pi.state_relax_1))
            df[Symbol(string_pi_mc_upper_bounds*"_$i")] = max.(df[Symbol(string_pi_mc_upper_bounds*"_$i")], cc.(lower_evaluator_mc_pi.state_relax_1))
        else
            df[Symbol(string_mc_lower_bounds*"_$i")] = min.(df[Symbol(string_mc_lower_bounds*"_$i")], cv.(lower_evaluator_mc.state_relax_n[i,:]))
            df[Symbol(string_mc_upper_bounds*"_$i")] = max.(df[Symbol(string_mc_upper_bounds*"_$i")], cc.(lower_evaluator_mc.state_relax_n[i,:]))
            df[Symbol(string_pi_mc_lower_bounds*"_$i")] = min(df[Symbol(string_pi_mc_lower_bounds*"_$i")], cv.(lower_evaluator_mc_pi.state_relax_n[i,:]))
            df[Symbol(string_pi_mc_upper_bounds*"_$i")] = max(df[Symbol(string_pi_mc_upper_bounds*"_$i")], cc.(lower_evaluator_mc_pi.state_relax_n[i,:]))
        end
    end

    # updates initial condition
    for i in 1:nx
        ic_opts[Symbol(string_mc_lower_bounds*"_$i")] = min.(ic_opts[Symbol(string_mc_lower_bounds*"_$i")], lower_evaluator_mc.IC_relaxations[i].cv)
        ic_opts[Symbol(string_mc_upper_bounds*"_$i")] = max.(ic_opts[Symbol(string_mc_upper_bounds*"_$i")], lower_evaluator_mc.IC_relaxations[i].cc)
        ic_opts[Symbol(string_pi_mc_lower_bounds*"_$i")] = min.(ic_opts[Symbol(string_pi_mc_lower_bounds*"_$i")], lower_evaluator_mc_pi.IC_relaxations[i].cv)
        ic_opts[Symbol(string_pi_mc_upper_bounds*"_$i")] = max.(ic_opts[Symbol(string_pi_mc_upper_bounds*"_$i")], lower_evaluator_mc_pi.IC_relaxations[i].cc)
    end
end

function save_bounds(name, ann, f, h, np, nx, nt, s, t_start, t_end, method,
                     pL, pU, xL, xU, x0; g = nothing, hj = nothing)

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
    string_pi_lower_bounds = "pi_lower_bounds"
    string_pi_upper_bounds = "pi_upper_bounds"

    df = DataFrame()
    for i in 1:nx
        df[Symbol(string_mc_lower_bounds*"_$i")] = Inf*ones(nt-1)
        df[Symbol(string_mc_upper_bounds*"_$i")] = -Inf*ones(nt-1)
        df[Symbol(string_pi_mc_lower_bounds*"_$i")] = Inf*ones(nt-1)
        df[Symbol(string_pi_mc_upper_bounds*"_$i")] = -Inf*ones(nt-1)
    end

    # create initial condition data frame
    ic_df = DataFrame()
    for i in 1:nx
        ic_df[Symbol(string_mc_lower_bounds*"_$i")] = [Inf]
        ic_df[Symbol(string_mc_upper_bounds*"_$i")] = [-Inf]
        ic_df[Symbol(string_pi_mc_lower_bounds*"_$i")] = [Inf]
        ic_df[Symbol(string_pi_mc_upper_bounds*"_$i")] = [-Inf]
    end

    # creates and build the evaluators
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
    EAGO.set_current_node!(lower_evaluator_mc, n_mc)
    y_mid_mc = (lower_vars + upper_vars)/2.0
    EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc, y_mid_mc)

    # runs upper evaluator to contract box
    n_mc = EAGO.NodeBB(lower_vars, upper_vars, -3.4, 2.1, 2, 1, true)
    EAGO.set_current_node!(upper_evaluator, n_mc)
    EAGO_Differential.relax_ode_implicit!(upper_evaluator)

    # stores the pi bounds to the dictionary enty
    for i in 1:nx
        if (nx == 1)
            df[Symbol(string_pi_lower_bounds*"_$i")] = min.(df[Symbol(string_pi_mc_lower_bounds*"_$i")], lo.(upper_evaluator.state_relax_1))
            df[Symbol(string_pi_upper_bounds*"_$i")] = max.(df[Symbol(string_pi_mc_upper_bounds*"_$i")], hi.(upper_evaluator.state_relax_1))
        else
            df[Symbol(string_pi_lower_bounds*"_$i")] = min(df[Symbol(string_pi_mc_lower_bounds*"_$i")], lo.(upper_evaluator.state_relax_n[i,:]))
            df[Symbol(string_pi_upper_bounds*"_$i")] = max(df[Symbol(string_pi_mc_upper_bounds*"_$i")], hi.(upper_evaluator.state_relax_n[i,:]))
        end
        ic_df[Symbol(string_pi_lower_bounds*"_$i")] = [lo(lower_evaluator_mc.IC_relaxations[i])]
        ic_df[Symbol(string_pi_upper_bounds*"_$i")] = [hi(lower_evaluator_mc.IC_relaxations[i])]
    end

    # gets post-contractor node
    lower_vars = Float64[]
    upper_vars = Float64[]
    if (nx == 1)
        append!(lower_vars, lo.(upper_evaluator.state_relax_1))
        append!(upper_vars, hi.(upper_evaluator.state_relax_1))
    else
        for i in 1:nt-1
            append!(lower_vars, lo.(upper_evaluator.state_relax_n[:,i]))
            append!(upper_vars, hi.(upper_evaluator.state_relax_n[:,i]))
        end
    end
    append!(lower_vars, lo.(upper_evaluator.P))
    append!(upper_vars, hi.(upper_evaluator.P))
    n_mc_pi = EAGO.NodeBB(lower_vars, upper_vars, -3.4, 2.1, 2, 1, true)

    #performs midpoint evaluation to generate reference point (McCormick + PI)
    EAGO.set_current_node!(lower_evaluator_mc_pi, n_mc_pi)
    y_mid_mc_pi = (lower_vars + upper_vars)/2.0
    EAGO_Differential.relax_ode_implicit!(lower_evaluator_mc_pi, y_mid_mc_pi)

    # performs main grid search
    if np == 1
        for px in range(pL[1], pU[1], length = grid_size)
            y = Float64[ones((nt-1)*nx); px]
            save_kernel!(df, ic_df, y, nx,
                         lower_evaluator_mc,
                         lower_evaluator_mc_pi,
                         upper_evaluator,
                         string_mc_lower_bounds,
                         string_mc_upper_bounds,
                         string_pi_mc_lower_bounds,
                         string_pi_mc_upper_bounds,
                         string_pi_lower_bounds,
                         string_pi_upper_bounds)
        end
    elseif np == 2
        range1 = range(pL[1], pU[1], length = grid_size)
        for px in range1
            range2 = range(pL[2], pU[2], length = grid_size)
            for py in range2
                y = Float64[ones((nt-1)*nx); px; py]
                save_kernel!(df, ic_df, y, nx,
                             lower_evaluator_mc,
                             lower_evaluator_mc_pi,
                             upper_evaluator,
                             string_mc_lower_bounds,
                             string_mc_upper_bounds,
                             string_pi_mc_lower_bounds,
                             string_pi_mc_upper_bounds,
                             string_pi_lower_bounds,
                             string_pi_upper_bounds)
            end
        end
    elseif np == 3
        for px in range(pL[1], pU[1], length = grid_size)
            for py in range(pL[2], pU[2], length = grid_size)
                for pz in range(pL[3], pU[3], length = grid_size)
                    y = Float64[ones((nt-1)*nx); px; py; pz]
                    save_kernel!(df, ic_df, y, nx,
                             lower_evaluator_mc,
                             lower_evaluator_mc_pi,
                             upper_evaluator,
                             string_mc_lower_bounds,
                             string_mc_upper_bounds,
                             string_pi_mc_lower_bounds,
                             string_pi_mc_upper_bounds,
                             string_pi_lower_bounds,
                             string_pi_upper_bounds)
                end
            end
        end
    end

    # saves trajectories to csv file
    append!(ic_df, df)
    ic_df.time = [i for i in range(t_start, stop = t_end, length = nt)]
    save_path_data = export_path*name*".csv"
    CSV.write(save_path_data, ic_df)

    # saves the annotations
    #annot = "$ann: method = $method of order $s, time from $t_start to $t_end"
    #ann_path_data = export_path*name*".txt"
    #annotate_file = open(ann_path_data, create = true)
    #write(annotate_file, annot)
    #close(annotate_file)

    return lower_evaluator_mc_pi
end
