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
    annotate_file = open(ann_path_data, create = true)
    write(annotate_file, annot)
    close(annotate_file)
end

function save_bounds(name, opt)
end
