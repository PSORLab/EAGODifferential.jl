mutable struct ImplicitODEUpperEvaluator <: MOI.AbstractNLPEvaluator
    current_node::NodeBB
    last_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool
    has_ineq::Bool
    objective_fun
    constraints_fun
    initial_condition_fun

    ivp::IVPInfo

    state_fun::Vector{Function}
    state_jac_fun::Vector{Function}

    np::Int
    ng::Int
    nx::Int
    jacobian_sparsity::Vector{Tuple{Int64,Int64}}
    obj_eval::Bool
    cnstr_eval::Bool
    init_relax_run::Bool

    P::Vector{IntervalType}
    X::Array{IntervalType,2}
    H::Vector{IntervalType}
    J::Array{IntervalType,2}
    Y::Array{Float64,2}

    var_relax::Vector{EAGO.IntervalType}
    IC_relaxations::Vector{EAGO.IntervalType}
    state_relax_1::Vector{EAGO.IntervalType}
    state_relax_n::Array{EAGO.IntervalType,2}
    constraint_relax::Vector{EAGO.IntervalType}
    obj_relax::EAGO.IntervalType

    kmax::Int
    etol::Float64
    rtol::Float64
    exclusion_flag::Bool
    inclusion_flag::Bool
    extended_division_flag::Bool

    inc::Vector{Bool}
    incLow::Vector{Bool}
    incHigh::Vector{Bool}

    N::Vector{EAGO.IntervalType}
    Ntemp::Vector{EAGO.IntervalType}
    Xi::Vector{EAGO.IntervalType}
    X1::Vector{EAGO.IntervalType}

    # timer
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    debug

    function ImplicitODEUpperEvaluator()
        d = new()

        d.objective_fun = nothing
        d.constraints_fun = nothing
        d.initial_condition_fun = nothing
        d.state_fun =  Function[]
        d.state_jac_fun = Function[]

        d.obj_eval = false
        d.cnstr_eval = false
        d.init_relax_run = false
        d.jacobian_sparsity = Tuple{Int64,Int64}[]

        d.np = 0
        d.nx = 0
        d.ng = 0
        d.has_ineq = false

        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0

        d.var_relax = fill(zero(EAGO.IntervalType),(1,))
        d.state_relax_1 = [zero(EAGO.IntervalType)]
        d.state_relax_n = fill(zero(EAGO.IntervalType), (2,2,))
        d.obj_relax = zero(EAGO.IntervalType)
        d.constraint_relax = [zero(EAGO.IntervalType)]
        d.P = fill(EAGO.IntervalType(0.0), (1,))
#        d.X = fill(EAGO.IntervalType(0.0), (2,2))

        d.H = [zero(EAGO.IntervalType)]
        d.J = fill(zero(EAGO.IntervalType), (2,2))
        d.Y = fill(0.0, (2,2))

        d.ivp = IVPInfo()
        d.IC_relaxations = [zero(EAGO.IntervalType)]

        d.disable_1storder = false
        d.disable_2ndorder = true
        d.has_nlobj = false

        d.inc = fill(false, (1,))
        d.incLow = fill(false, (1,))
        d.incHigh = fill(false, (1,))

        d.N = fill(EAGO.IntervalType(0.0), (1,))
        d.Ntemp = fill(EAGO.IntervalType(0.0), (1,))
        d.X1 = fill(EAGO.IntervalType(0.0), (1,))
        d.Xi = fill(EAGO.IntervalType(0.0), (1,))

        d.kmax = 50
        d.etol = 1E-12
        d.rtol = 1E-12
        d.exclusion_flag = false
        d.inclusion_flag = false
        d.extended_division_flag = false

        d.current_node = NodeBB()
        d.last_node = NodeBB()

        d.debug = NaN

        return d
    end
end

"""
    build_evaluator!

A function that builds the ImplicitODELowerEvaluator form user-supplied functions.
Note that by convention the objective function and constraint function are specified in
terms of discretized indices. This is done to allow for common objectives and constraints
encountered in optimal control formulations.
"""
function EAGO.build_evaluator!(d::ImplicitODEUpperEvaluator, f::Function, h::Function, np::Int, nx::Int,
                               nt::Int, s::Int, t_start::Float64, t_end::Float64, method::Symbol,
                               pL::Vector{Float64}, pU::Vector{Float64}, xL::Vector{Float64}, xU::Vector{Float64}, ic::Function;
                               g = nothing, hj = nothing)

    # setup objective and constraint functions
    d.has_nlobj = true
    d.objective_fun = f
    d.initial_condition_fun = ic

    # sets up initial value problem parameters
    d.ivp.time_start = t_start
    d.ivp.time_end = t_end
    d.ivp.time_steps = nt
    d.ivp.step_size = (t_end - t_start)/(nt - 1)
    d.ivp.time = zeros(nt)
    d.ivp.time[1] = t_start
    for i in 2:nt
        d.ivp.time[i] = d.ivp.time[i-1] + d.ivp.step_size
    end
    d.ivp.method = method
    d.ivp.method_order = s

    # add h and hj functions formated for AM and BDF methods
    if d.ivp.method == :BDF
        push!(d.state_fun, (hout,xout,x,p,t) -> bdf_kernel_1!(h, hout, d.ivp.step_size, xout, x[1:nx], p, t[1]))
        push!(d.state_fun, (hout,xout,x,p,t) -> bdf_kernel_2!(h, hout, d.ivp.step_size, xout, x[1:nx], x[(nx+1):(2*nx)], p, t[1]))
        push!(d.state_jac_fun, (jout,xout,x,p,t) -> bdf_jac_kernel_1!(hj, jout, d.ivp.step_size, xout, p, t[1]))
        push!(d.state_jac_fun, (jout,xout,x,p,t) -> bdf_jac_kernel_2!(hj, jout, d.ivp.step_size, xout, p, t[1]))
    elseif d.ivp.method == :AM
        push!(d.state_fun, (hout,xout,x,p,t) -> am_kernel_1!(h, hout, d.ivp.step_size, xout, x[1:nx], p, t[1]))
        push!(d.state_fun, (hout,xout,x,p,t) -> am_kernel_2!(h, hout, d.ivp.step_size, xout, x[1:nx], p, t[1], t[2]))
        push!(d.state_jac_fun, (jout,xout,x,p,t) -> am_jac_kernel_1!(hj, jout, d.ivp.step_size, xout, p, t[1]))
        push!(d.state_jac_fun, (jout,xout,x,p,t) -> am_jac_kernel_2!(hj, jout, d.ivp.step_size, xout, p, t[1]))
    end

    # set implicit routine information

    # preallocates the storage variables
    temp = zero(EAGO.IntervalType)

    # get dimension sizes
    d.np = np; d.nx = nx;
    if ~(g === nothing)
        d.has_ineq = true
        d.constraints_fun = g
        ng = length(g([1.0 for i=1:nx, j=1:(nt-1)], ones(nx), ones(np), d.ivp.time))
        d.ng = ng
        d.constraint_relax = fill(temp, (ng,))
    end

    d.var_relax = fill(temp, (np,))
    if (nx == 1)
        d.state_relax_1 = fill(temp, (nt-1,))
    else
        d.state_relax_n = fill(temp, (nx, nt-1))
    end

    # storage for intermediate calculations
    d.IC_relaxations = fill(temp, (nx,))
    d.H = fill(temp, (nx,))
    if (nx > 1)
        d.J = zeros(EAGO.IntervalType, nx, nx)
        d.Y = zeros(Float64, nx, nx)
    else
        d.J = zeros(EAGO.IntervalType, 2, 2)
        d.Y = zeros(Float64, 2, 2)
    end

    d.P = fill(IntervalType(0.0), (np,))
    for i in 1:np
        d.P[i] = IntervalType(pL[i], pU[i])
    end

    d.N = fill(IntervalType(0.0), (nx,))
    d.Ntemp = fill(IntervalType(0.0), (nx,))

#    d.X = fill(IntervalType(0.0), (nx, nt-1))
    d.X1 = fill(IntervalType(0.0), (nx,))
    d.Xi = fill(IntervalType(0.0), (nx,))

    d.inc = fill(false, (nx,))
    d.incLow = fill(false, (nx,))
    d.incHigh = fill(false, (nx,))

#    indx = 1
#    for i in 1:(nt-1)
#        for j in 1:nx
#            d.X[j,i] = IntervalType(xL[j], xU[j])
#        end
#        indx += 1
#    end
end

function EAGO.set_current_node!(x::ImplicitODEUpperEvaluator,n::NodeBB)
    x.current_node = copy(n)
end

function EAGO.set_last_node!(x::ImplicitODEUpperEvaluator,n::NodeBB)
    x.last_node =  copy(n)
end

# LOOKS GREAT!
function relax_objective!(d::ImplicitODEUpperEvaluator)
    if ~d.obj_eval
        if d.nx == 1
            d.obj_relax = d.objective_fun(d.state_relax_1, d.IC_relaxations, d.var_relax, d.ivp.time)
        else
            d.obj_relax =  d.objective_fun(d.state_relax_n, d.IC_relaxations, d.var_relax, d.ivp.time)
        end
        d.obj_eval = true
    end
end

# LOOKS GREAT!
function relax_constraints!(d::ImplicitODEUpperEvaluator)
    if d.has_ineq
        if ~d.cnstr_eval
            if d.nx == 1
                d.constraint_relax[:] = d.constraints_fun(d.state_relax_1, d.IC_relaxations, d.var_relax, d.ivp.time)
            else
                d.constraint_relax[:] = d.constraints_fun(d.state_relax_n, d.IC_relaxations, d.var_relax, d.ivp.time)
            end
            d.cnstr_eval = true
        end
    end
end

num_state_variables(d::ImplicitODEUpperEvaluator) = d.nx*(d.nt-1)
num_decision_variables(d::ImplicitODEUpperEvaluator) = d.np

include("upper_calculations.jl")
include("upper_info.jl")
