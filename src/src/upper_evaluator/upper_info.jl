# LOOKS GREAT!
function MOI.eval_objective(d::ImplicitODEUpperEvaluator, y)
    d.eval_objective_timer += @elapsed begin
        val = zero(eltype(y))
        if d.has_nlobj
            relax_ode_implicit!(d)
            relax_objective!(d)
            val = d.obj_relax.hi
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

# LOOKS GREAT!
function MOI.eval_constraint(d::ImplicitODEUpperEvaluator, g, y)
    d.eval_constraint_timer += @elapsed begin
        if d.ng > 0
            relax_ode_implicit!(d)
            relax_constraints!(d)
            for i in 1:d.ng
                g[i] = d.cnstr_relax[i].lo
            end
        end
    end
    return
end

# LOOKS GREAT!
function MOI.eval_objective_gradient(d::ImplicitODEUpperEvaluator, df, y)
    d.eval_objective_timer += @elapsed begin
        if d.has_nlobj
            fill!(df, 0.0)
        else
            error("No nonlinear objective.")
        end
    end
    return
end

# LOOKS GREAT!
function MOI.jacobian_structure(d::ImplicitODEUpperEvaluator)
    # use user-defined sparsity pattern if possible
    if length(d.jacobian_sparsity) > 0
        return d.jacobian_sparsity
    else
        d.jacobian_sparsity = Tuple{Int64,Int64}[(row, idx) for row in 1:d.ng for idx in 1:d.np]
        return d.jacobian_sparsity
    end
end

# LOOKS GREAT!
function MOI.hessian_lagrangian_structure(d::ImplicitODEUpperEvaluator)
    error("Hessian computations not currently supported by Implicit optimizer.")
end

# LOOKS GREAT!
function _hessian_lagrangian_structure(d::ImplicitODEUpperEvaluator)
    error("Hessian lagrangian structure not supported by Implicit optimizer.")
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian(d::ImplicitODEUpperEvaluator,g,y)
    d.eval_constraint_jacobian_timer += @elapsed begin
        if d.ng > 0
            fill!(g, 0.0)
        end
    end
    return g
end

# LOOKS GREAT!
function MOI.eval_constraint_jacobian_product(d::ImplicitODEUpperEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            if d.ng > 0
                fill!(out, 0.0)
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# PROBABLY DONE
function MOI.eval_constraint_jacobian_transpose_product(d::ImplicitODEUpperEvaluator, out, y, w)
    if (!d.disable_1storder)
        d.eval_constraint_jacobian_timer += @elapsed begin
            fill!(out,0.0)
        end
    else
        error("First order information unavailable.")
    end
    return
end

# LOOKS GREAT
function MOI.features_available(d::ImplicitODEUpperEvaluator)
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
function MOI.initialize(d::ImplicitODEUpperEvaluator, requested_features::Vector{Symbol}) end

# LOOKS GREAT
MOI.objective_expr(d::ImplicitODEUpperEvaluator) = error("EAGO.ImplicitODEUpperEvaluator doesn't provide expression graphs of constraint functions.")
#LOOKS GREAT
MOI.constraint_expr(d::ImplicitODEUpperEvaluator) = error("EAGO.ImplicitODEUpperEvaluator doesn't provide expression graphs of constraint functions.")
