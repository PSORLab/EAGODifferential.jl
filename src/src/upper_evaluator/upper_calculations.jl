"""
    implicit_ode_calc_functions!(d::ImplicitODEUpperEvaluator,y)

Main calculation kernel used for midpoint upper bounding calculation.
"""
function relax_ode_implicit!(d::ImplicitODEUpperEvaluator)
    # Generate new parameters for implicit relaxation if necessary
    if ~EAGO.same_box(d.current_node, d.last_node, 0.0)
        println("ran ode upper kernel")
        t = [0.0; 0.0]
        nx = d.nx
        np = d.np
        nt = d.ivp.time_steps

        kmax = d.kmax
        etol = d.etol
        rtol = d.rtol

        d.last_node = d.current_node
        d.obj_eval = false
        d.cnstr_eval = false

        lower_var_bnds::Vector{Float64} = d.current_node.lower_variable_bounds
        upper_var_bnds::Vector{Float64} = d.current_node.upper_variable_bounds
        d.P[:] = EAGO.IntervalType.(lower_var_bnds[(nx*(nt-1)+1):(nx*(nt-1)+np)],
                                    upper_var_bnds[(nx*(nt-1)+1):(nx*(nt-1)+np)])
        d.IC_relaxations[:] = d.initial_condition_fun(d.P)

        if (nx == 1)

            fill!(d.state_relax_1, EAGO.IntervalType(lower_var_bnds[1], upper_var_bnds[1]))
            t[1] = d.ivp.time[2]
            Eflg, Iflg, eDflg = EAGO.param_intv_contractor!(d.state_fun[1], d.state_jac_fun[1],
                                                       view(d.state_relax_1, 1:1), d.Ntemp, d.N, d.Xi,
                                                       d.X1, d.IC_relaxations, t,
                                                       d.Y, d.J, d.H, d.P,
                                                       d.inc, d.incLow, d.incHigh,
                                                       nx, kmax, etol, rtol)
            if Eflg
                d.exclusion_flag = Eflg
            end
            if eDflg
                d.extended_division_flag = eDflg
            end
            if Eflg || eDflg
               return
            end

            s = 1
            t = zeros(2)
            x0 = EAGO.IntervalType[]
            for i in 2:(nt-1)
                t[1] = d.ivp.time[i+1]
                t[2] = d.ivp.time[i]
                s = min(i, d.ivp.method_order)
                empty!(x0)
                if (s == i)
                    for j in 1:(s-1)
                        append!(x0, d.state_relax_1[i-j])
                    end
                    append!(x0, d.IC_relaxations)
                else
                    for j in 1:s
                        append!(x0, d.state_relax_1[i-j])
                    end
                end
                Eflg, Iflg, eDflg = EAGO.param_intv_contractor!(d.state_fun[1], d.state_jac_fun[1],
                                                          view(d.state_relax_1, i:i), d.Ntemp, d.N, d.Xi,
                                                          d.X1, x0, t, d.Y, d.J, d.H, d.P,
                                                          d.inc, d.incLow, d.incHigh,
                                                          nx, kmax, etol, rtol)
                if Eflg
                    d.exclusion_flag = Eflg
                end
                if eDflg
                    d.extended_division_flag = eDflg
                end
                if Eflg || eDflg
                    return
                end
            end
        else
            for i in 1:(nt-1)
                d.state_relax_n[:,i] .= EAGO.IntervalType.(lower_var_bnds[(nx*(i-1)+1):nx*i],
                                                           upper_var_bnds[(nx*(i-1)+1):nx*i])
            end
            t[1] = d.ivp.time[2]
            Eflg, Iflg, eDflg = EAGO.param_intv_contractor!(d.state_fun[1], d.state_jac_fun[1],
                                                            view(d.state_relax_n, 1:nx, 1:1), d.Ntemp,
                                                            d.N, d.Xi, d.X1, d.IC_relaxations, t, d.Y,
                                                            d.J, d.H, d.P, d.inc, d.incLow, d.incHigh,
                                                            nx, kmax, etol, rtol)
            if Eflg
                d.exclusion_flag = Eflg
            end
            if eDflg
                d.extended_division_flag = eDflg
            end
            if Eflg || eDflg
                return
            end

            x0 = EAGO.IntervalType[]
            for i in 2:(nt-1)
                t[1] = d.ivp.time[i+1]
                t[2] = d.ivp.time[i]
                s = min(i, d.ivp.method_order)
                empty!(x0)
                if (s == i)
                    for j in 1:(s-1)
                        append!(x0, d.state_relax_n[:, i-j])
                    end
                    append!(x0, d.IC_relaxations)
                else
                    for j in 1:s
                        append!(x0, d.state_relax_n[:, i-j])
                    end
                end
                Eflg, Iflg, eDflg = EAGO.param_intv_contractor!(d.state_fun[s], d.state_jac_fun[s],
                                                           view(d.state_relax_n, 1:nx, i:i), d.Ntemp,
                                                           d.N, d.Xi, d.X1, x0, t, d.Y, d.J, d.H, d.P,
                                                           d.inc, d.incLow, d.incHigh,
                                                           nx, kmax, etol, rtol)
                if Eflg
                    d.exclusion_flag = Eflg
                end
                if eDflg
                    d.extended_division_flag = eDflg
                end
                if Eflg || eDflg
                    return
                end
            end
        end
    end
end
# WILL USE THIS FOR PARAMETRIC INTERVAL CONTRACTOR... AND MOI FOR UPPER PROBLEM

function MOI.eval_objective(d::ImplicitODEUpperEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = 0.0
        if (d.has_nlobj)
            relax_ode_implicit!(d)
            val = d.objective_fun(d.state_relax_1, d.IC_relaxations, d.P, d.ivp.time)
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

function MOI.eval_constraint(d::ImplicitODEUpperEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if (d.ng) > 0
            relax_ode_implicit!(d)
            g[:] = d.constraints_fun(d.state_relax_1, d.IC_relaxations, d.P, d.ivp.time)
        end
    end
    return
end

MOI.features_available(d::ImplicitODEUpperEvaluator) = Symbol[]
function MOI.initialize(d::ImplicitODEUpperEvaluator, requested_features::Vector{Symbol}) end
