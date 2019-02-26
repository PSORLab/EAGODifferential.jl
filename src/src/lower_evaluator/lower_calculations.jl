function relax_ode_implicit!(d::ImplicitODELowerEvaluator, y)

    t = [0.0; 0.0]
    nx = d.nx
    np = d.np
    kmax = d.opts.kmax
    nt = d.ivp.time_steps

    lower_var_bnds = d.current_node.lower_variable_bounds
    upper_var_bnds = d.current_node.upper_variable_bounds
    # Generate new parameters for implicit relaxation if necessary
    if ~EAGO.same_box(d.current_node, d.last_node, 0.0)

        println("generate parameters")
        d.last_node = d.current_node
        d.obj_eval = false
        d.cnstr_eval = false

        indx = 1
        for j in 1:(nt-1)
            d.X[1:nx,j] = EAGO.IntervalType.(lower_var_bnds[((j-1)*nx+1):(j*nx)],
                                             upper_var_bnds[((j-1)*nx+1):(j*nx)])
        end

        shift = nx*(nt-1)+1
        d.P[:] = EAGO.IntervalType.(lower_var_bnds[shift:end], upper_var_bnds[shift:end])
        d.ref_p[:] = (lower_var_bnds[shift:end] + upper_var_bnds[shift:end])/2.0
        d.pref_mc[:] = MC{np}.(y, d.P, 1:np)
        d.IC_relaxations[:] = d.initial_condition_fun(d.pref_mc)

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
            println("relax time t[1] = $(t[1])")

            println("initial condition relaxations")
            for j in 1:nx
                print_out = d.IC_relaxations[j]
                grad_string = "grad ic[$j](cv, cc)[1,2,3] =  ($(round(print_out.cv_grad[1]; digits=3)), $(round(print_out.cc_grad[1]; digits=3))), ($(round(print_out.cv_grad[2]; digits=3)), $(round(print_out.cc_grad[2]; digits=3))), ($(round(print_out.cv_grad[3]; digits=3)), $(round(print_out.cc_grad[3]; digits=3))) "
                println("ic[$j](cv, cc, Intv) =  ($(print_out.cv),$(print_out.cc),[$(round(print_out.Intv.lo; digits=3)),$(round(print_out.Intv.hi; digits=3))])                     "*grad_string)
            end

            gen_expansion_params!(d.state_fun[1], d.state_jac_fun[1], d.pref_mc,
                                  d.IC_relaxations, view(d.state_relax_n, 1:nx, 1:1),
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                                  d.P, d.opts, view(d.state_ref_relaxation_n, 1:nx, 1:kmax, 1:1),
                                  d.H, d.J, d.Y, true, t, true)
            xref = view(d.state_ref_relaxation_n, 1:nx, 1:kmax, 1:1)
            for i in 1:kmax
                println("parameter #$i")
                for j in 1:nx
                    print_out = xref[j,i]
                    println("x[$j](cv, cc, Intv) =  ($(print_out.cv),$(print_out.cc),[$(round(print_out.Intv.lo; digits=3)),$(round(print_out.Intv.hi; digits=3))])")
                end
                for j in 1:nx
                    print_out = xref[j,i]
                    println("x[$j](cv, cc)[1,2,3] =  ($(round(print_out.cv_grad[1]; digits=3)), $(round(print_out.cc_grad[1]; digits=3))), ($(round(print_out.cv_grad[2]; digits=3)), $(round(print_out.cc_grad[2]; digits=3))), ($(round(print_out.cv_grad[3]; digits=3)), $(round(print_out.cc_grad[3]; digits=3))) ")
                end
            end

            x0 = MC{np}[]
            for i in 2:(nt-1)
                t[1] = d.ivp.time[i+1]
                t[2] = d.ivp.time[i]
                println("relax time t[$i] = $(t)")
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
                xref = view(d.state_ref_relaxation_n, 1:nx, 1:kmax, 1:1)
                #println("xref: $xref")
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
        temp = MC{np}.(y, d.P, 1:np)
        #println("evaluate at new MC point: $temp")
        d.var_relax[:] = MC{np}.(y, d.P, 1:np)
        d.IC_relaxations[:] = d.initial_condition_fun(d.var_relax)
        println("d.var_relax: $(d.var_relax)")
        if (nx == 1)
            t[1] = d.ivp.time[2]
            println("relax time t[1] == $(t[1])")
            implicit_relax_h!(d.state_fun[1], d.state_jac_fun[1], d.var_relax, d.pref_mc,
                              d.IC_relaxations, view(d.state_relax_1, 1:1),
                              d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                              d.P, d.opts, view(d.state_ref_relaxation_1, 1:kmax, 1:1),
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
                                  view(d.state_ref_relaxation_1, 1:kmax, i:i),
                                  d.H, d.J, d.Y, true, t, true)
            end
        else
            println("relax time 1")
            t[1] = d.ivp.time[2]
            implicit_relax_h!(d.state_fun[1], d.state_jac_fun[1], d.var_relax, d.pref_mc,
                              d.IC_relaxations, view(d.state_relax_n, 1:nx, 1:1),
                              d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,1],
                              d.P, d.opts, view(d.state_ref_relaxation_n, 1:nx, :, 1:1),
                              d.H, d.J, d.Y, true, t, true)

            x0 = MC{np}[]
            for i in 2:(nt-1)
                println("relax time $i")
                t[1] = d.ivp.time[i+1]
                t[2] = d.ivp.time[i]
                println("time 2 assigned: $(t[1]), time 1 assigned: $(t[2])")
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
                println("x0: $x0")
                implicit_relax_h!(d.state_fun[s], d.state_jac_fun[s], d.var_relax, d.pref_mc,
                                  x0, view(d.state_relax_n, 1:nx, i:i),
                                  d.xa_mc, d.xA_mc, d.z_mc, d.aff_mc, d.X[1:nx,i],
                                  d.P, d.opts, view(d.state_ref_relaxation_n, 1:nx, :, i:i),
                                  d.H, d.J, d.Y, true, t, true)
            end
        end
    end
end
