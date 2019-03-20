# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitODELowerEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = zero(eltype(y))
        if d.has_nlobj
            #try
                relax_ode_implicit!(d,y)
                relax_objective!(d)
                val = d.obj_relax.cv
            #catch
            #end
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitODELowerEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            #try
            relax_ode_implicit!(d,y)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.constraint_relax[i].cv
            end
            #catch
            #end
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitODELowerEvaluator, df, y)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            try
                relax_ode_implicit!(d,y)
                relax_objective!(d)
                for j in 1:d.np
                    df[j] = d.obj_relax.cv_grad[j]
                end
            catch
                for j in 1:d.np
                    df[j] = 0.0
                end
            end
        else
            error("No nonlinear objective.")
        end
    end
    return
end

# LOOKS GREAT!
function MOI.jacobian_structure(d::ImplicitODELowerEvaluator)
    # use user-defined sparsity pattern if possible
    if length(d.jacobian_sparsity) > 0
        return d.jacobian_sparsity
    else
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.ng for idx in 1:d.np]
        return d.jacobian_sparsity
    end
end

# LOOKS GREAT!
function MOI.hessian_lagrangian_structure(d::ImplicitODELowerEvaluator)
    error("Hessian computations not currently supported by Implicit optimizer.")
end

# LOOKS GREAT!
function _hessian_lagrangian_structure(d::ImplicitODELowerEvaluator)
    error("Hessian lagrangian structure not supported by Implicit optimizer.")
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian(d::ImplicitODELowerEvaluator,g,y)
    d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            relax_ode_implicit!(d,y)
            relax_constraints!(d)
            fill!(g, 0.0)
            if (d.ng == 1) && (d.np == 1)
                g[1] = d.constraint_relax[1].cv_grad[1]
            else
                for (i,j) in d.jacobian_sparsity
                    temp = d.constraint_relax[i].cv_grad[j]
                    g[i,j] = d.constraint_relax[i].cv_grad[j]
                end
            end
        end
    end
    return g
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian_product(d::ImplicitODELowerEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                relax_ode_implicit!(d,y)
                relax_constraints!(d)
                fill!(out, 0.0)
                for (i,j) in d.jacobian_sparsity
                    out[i] += d.constraint_relax[i].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# PROBABLY DONE
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitODELowerEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                relax_ode_implicit!(d,y)
                relax_constraints!(d)
                fill!(out, 0.0)
                for (i,j) in d.jacobian_sparsity
                    out[i] += d.constraint_relax[i].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# LOOKS GREAT
function MOI.features_available(d::ImplicitODELowerEvaluator)
    features = Symbol[]
    if !d.disable_1storder
        push!(features,:Grad)
        push!(features,:Jac)
    end
    if !d.disable_2ndorder
        push!(features,:Hess)
        push!(features,:HessVec)
    end
    return features
end

# LOOKS GREAT
function MOI.initialize(d::ImplicitODELowerEvaluator, requested_features::Vector{Symbol}) end

# LOOKS GREAT
MOI.objective_expr(d::ImplicitODELowerEvaluator) = error("EAGO.ImplicitODELowerEvaluator doesn't provide expression graphs of constraint functions.")
#LOOKS GREAT
MOI.constraint_expr(d::ImplicitODELowerEvaluator) = error("EAGO.ImplicitODELowerEvaluator doesn't provide expression graphs of constraint functions.")

#=
function eval_constraint_cc(d::ImplicitODELowerEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            relax_ode_implicit!(d,y)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.constraint_relax[i].cc
            end
        end
    end
    return
end
=#

function eval_constraint_cc(d::ImplicitODELowerEvaluator, g::Vector{Float64}, y::Vector{Float64})
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            relax_ode_implicit!(d,y)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.constraint_relax[i].cc
            end
        end
    end
    return
end

function eval_constraint_cc_grad(d::ImplicitODELowerEvaluator, g, y)
    d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            relax_ode_implicit!(d,y)
            relax_constraints!(d)
            fill!(g, 0.0)
            if (d.ng == 1) && (d.np == 1)

                g[1] = d.constraint_relax[1].cc_grad[1]
            else
                for (i,j) in d.jacobian_sparsity
                    temp = d.constraint_relax[i].cc_grad[j]
                    g[i,j] = d.constraint_relax[i].cc_grad[j]
                end
            end
        end
    end
    return g
end
