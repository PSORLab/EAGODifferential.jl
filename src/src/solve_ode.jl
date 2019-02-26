function interval_preprocess_ode(x::EAGO.Optimizer, y::EAGO.NodeBB)

    Eflag = false
    Iflag = false
    eDflag = false

    set_current_node!(x.nlp_data.evaluator, y)
    nx = x.state_variables
    np = x.variable_number - x.state_variables
    current_node = x.nlp_data.evaluator.current_node

    g = zeros(x.nlp_data.evaluator.ng)
    MOI.eval_constraint(x.nlp_data.evaluator, g, y)
    Eflag = x.nlp_data.evaluator.exclusion_flag
    if ~Eflag
        EFlag = any(i-> (i > 0), g)
    end

    x.current_preprocess_info.feasibility = ~Eflag
    if ~Eflag # TODO: Update node if feasible....
        y.lower_variable_bounds[1:nx] = lo.(X1)
        y.upper_variable_bounds[1:nx] = hi.(X1)
    end
end

function midpoint_upper_bnd_ode(x::Optimizer,y::NodeBB)
    if is_integer_feasible(x) #&& mod(x.CurrentIterationCount,x.UpperBoundingInterval) == 1

        node_ymid = # TODO: create midpoint node
        set_current_node!(x.nlp_data.evaluator, node_ymid)
        set_to_mid!(x.current_upper_info.solution, node_ymid)
        if x.nlp_data.evaluator.ng > 0
            g = zeros(x.nlp_data.evaluator.ng)
            MOI.eval_constraint(x.nlp_data.evaluator, g, x.current_upper_info.solution)
            result_status = any(i-> (i > 0), g) ? MOI.INFEASIBLE_POINT : MOI.FEASIBLE_POINT
            if (result_status == MOI.FEASIBLE_POINT)
                result_status = x.nlp_data.evaluator.exclusion_flag ? MOI.INFEASIBLE_POINT : MOI.FEASIBLE_POINT
            end
        else
            result_status = MOI.FEASIBLE_POINT
        end
        if (result_status == MOI.FEASIBLE_POINT)
            x.current_upper_info.feasibility = true
            val = MOI.eval_objective(x.nlp_data.evaluator, x.current_upper_info.solution)
            x.current_upper_info.value = val
        else
            x.current_upper_info.feasibility = false
            x.current_upper_info.value = Inf
        end
    else
        x.current_upper_info.feasibility = false
        x.current_upper_info.value = Inf
    end
end

"""
    solve_ode

Solves the optimization problem `min_{x,p} f(x,p,t)` with respect to_indices
`dx/dt = h(x,p,t)` on `x in X` and `p in P`.
"""
function solve_ode(f, h, hj, X, P, time_start, time_end, time_steps)
end
