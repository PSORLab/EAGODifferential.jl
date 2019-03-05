function relax_ode_implicit!(d::ImplicitODELowerEvaluator, y)

    t = [0.0; 0.0]
    nx = d.nx
    np = d.np
    kmax = d.opts.kmax
    nt = d.ivp.time_steps
    shift = nx*(nt-1)+1

    lower_var_bnds = d.current_node.lower_variable_bounds
    upper_var_bnds = d.current_node.upper_variable_bounds
    # Generate new parameters for implicit relaxation if necessary

    current_node = d.current_node
    last_node = d.last_node
    new_box_flag = ~EAGO.same_box(d.current_node, d.last_node, 0.0)
    if new_box_flag
        d.last_node = d.current_node
        d.obj_eval = false
        d.cnstr_eval = false

        indx = 1
        for j in 1:(nt-1)
            d.X[1:nx,j] = EAGO.IntervalType.(lower_var_bnds[((j-1)*nx+1):(j*nx)],
                                             upper_var_bnds[((j-1)*nx+1):(j*nx)])
        end

        P = EAGO.IntervalType.(lower_var_bnds[shift:end], upper_var_bnds[shift:end])
        ref_p = (lower_var_bnds[shift:end] + upper_var_bnds[shift:end])/2.0
        pref_mc = MC{np}.(y[shift:end], d.P, 1:np)

        d.P[:] = EAGO.IntervalType.(lower_var_bnds[shift:end], upper_var_bnds[shift:end])
        d.ref_p[:] = (lower_var_bnds[shift:end] + upper_var_bnds[shift:end])/2.0
        d.pref_mc[:] = MC{np}.(y[shift:end], d.P, 1:np)
        d.IC_relaxations[:] = d.initial_condition_fun(d.pref_mc)
        d.state_update(d.IC_relaxations[:])

        s = 1
        if (nx == 1)
            t[1] = d.ivp.time[2]
            gen_expansion_params!(d.state_fun[1], d.state_jac_fun[1], d.pref_mc,
                                  d.IC_relaxations, view(d.state_relax_1, 1:1),
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                                  d.P, d.opts, view(d.state_ref_relaxation_1, 1:1, 1:kmax, 1:1),
                                  d.H, d.J, d.Y, true, t, true)

            x0 = MC{np}[]
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
                gen_expansion_params!(d.state_fun[s], d.state_jac_fun[s], d.pref_mc,
                                      x0, view(d.state_relax_1, i:i), d.xa_mc, d.xA_mc,
                                      d.z_mc, d.aff_mc, d.X[1:nx,i], d.P, d.opts,
                                      view(d.state_ref_relaxation_1, 1:1, 1:kmax, i:i),
                                      d.H, d.J, d.Y, true, t, true)
            end
        else
            t[1] = d.ivp.time[2]
            gen_expansion_params!(d.state_fun[1], d.state_jac_fun[1], d.pref_mc,
                                  d.IC_relaxations, view(d.state_relax_n, 1:nx, 1:1),
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                                  d.P, d.opts, view(d.state_ref_relaxation_n, 1:nx, 1:kmax, 1:1),
                                  d.H, d.J, d.Y, true, t, true)
            d.state_update(view(d.state_relax_n, 1:nx, 1:1))
            x0 = MC{np}[]
            for i in 2:(nt-1)
                t[1] = d.ivp.time[i+1]
                t[2] = d.ivp.time[i]
                s = min(i, d.ivp.method_order)
                empty!(x0)
                if (s == i)
                    for j in 1:(s-1)
                        append!(x0, d.state_relax_n[:,i-j])
                    end
                    append!(x0, d.IC_relaxations)
                else
                    for j in 1:s
                        append!(x0, d.state_relax_n[:,i-j])
                    end
                end
                gen_expansion_params!(d.state_fun[s], d.state_jac_fun[s], d.pref_mc,
                                      x0, view(d.state_relax_n, 1:nx, i:i), d.xa_mc, d.xA_mc,
                                      d.z_mc, d.aff_mc, d.X[1:nx,i], d.P, d.opts,
                                      view(d.state_ref_relaxation_n, 1:nx, 1:kmax, i:i),
                                      d.H, d.J, d.Y, true, t, true)
                d.state_update(view(d.state_relax_n, 1:nx, i:i))
                xref = view(d.state_ref_relaxation_n, 1:nx, 1:kmax, 1:1)
            end
        end
    end
    # Generate new value of implicit relaxation
    new_point_flag = false

    for i in 1:np
        if d.ref_p[i] != y[i]
            new_point_flag = true; break
        end
    end
    if new_point_flag && (~d.init_relax_run)
        d.obj_eval = false
        d.cnstr_eval = false
        temp = MC{np}.(y[shift:end], d.P, 1:np)
        d.var_relax[:] = MC{np}.(y[shift:end], d.P, 1:np)
        d.IC_relaxations[:] = d.initial_condition_fun(d.var_relax)
        if (nx == 1)
            t[1] = d.ivp.time[2]
            implicit_relax_h!(d.state_fun[1], d.state_jac_fun[1], d.var_relax, d.pref_mc,
                              d.IC_relaxations, view(d.state_relax_1, 1:nx, 1:1),
                              d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                              d.P, d.opts, view(d.state_ref_relaxation_1, 1:1, 1:kmax, 1:1),
                              d.H, d.J, d.Y, true, t, true)

            s = 1
            t = zeros(2)
            x0 = MC{np}[]
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
                implicit_relax_h!(d.state_fun[s], d.state_jac_fun[s], d.var_relax, d.pref_mc,
                                  x0, view(d.state_relax_1, i:i), d.xa_mc, d.xA_mc,
                                  d.z_mc, d.aff_mc, d.X[1:nx,i], d.P, d.opts,
                                  view(d.state_ref_relaxation_1, 1:1, 1:kmax, i:i),
                                  d.H, d.J, d.Y, true, t, true)
            end
        else
            t[1] = d.ivp.time[2]
            implicit_relax_h!(d.state_fun[1], d.state_jac_fun[1], d.var_relax, d.pref_mc,
                              d.IC_relaxations, view(d.state_relax_n, 1:nx, 1:1),
                              d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                              d.P, d.opts, view(d.state_ref_relaxation_n, 1:nx, :, 1:1),
                              d.H, d.J, d.Y, true, t, true)
            d.state_update(view(d.state_relax_n, 1:nx, 1:1))
            x0 = MC{np}[]
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
                implicit_relax_h!(d.state_fun[s], d.state_jac_fun[s], d.var_relax, d.pref_mc,
                                  x0, view(d.state_relax_n, 1:nx, i:i),
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,i],
                                  d.P, d.opts, view(d.state_ref_relaxation_n, 1:nx, :, i:i),
                                  d.H, d.J, d.Y, true, t, true)
                d.state_update(view(d.state_relax_n, 1:nx, i:i))
            end
        end
    end
end
