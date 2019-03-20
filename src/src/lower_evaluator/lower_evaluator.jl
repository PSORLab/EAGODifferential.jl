mutable struct ImplicitODELowerEvaluator{N} <: MOI.AbstractNLPEvaluator
    current_node::NodeBB
    last_node::NodeBB

    disable_1storder::Bool
    disable_2ndorder::Bool
    has_nlobj::Bool
    has_ineq::Bool
    objective_fun
    state_update
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

    xa_mc::Vector{MC{N}}
    xA_mc::Vector{MC{N}}
    z_mc::Vector{MC{N}}
    aff_mc::Vector{MC{N}}
    H::Vector{MC{N}}
    J::Array{MC{N},2}
    Y::Array{Float64,2}

    last_p::Vector{Float64}
    ref_p::Vector{Float64}
    pref_mc::Vector{MC{N}}
    var_relax::Vector{MC{N}}

    IC_relaxations::Vector{MC{N}}

    state_relax_1::Vector{MC{N}}                     # Moved to 2D
    state_ref_relaxation_1::Array{MC{N},3}           # Moved to 3D

    state_relax_n::Array{MC{N},2}
    state_ref_relaxation_n::Array{MC{N},3}

    constraint_relax::Vector{MC{N}}
    obj_relax::MC{N}

    opts::mc_opts

    objective_ubd::Float64

    # timer
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    debug

    function ImplicitODELowerEvaluator{N}() where N
        d = new()

        d.objective_fun = nothing
        d.constraints_fun = nothing
        d.state_update = nothing
        d.initial_condition_fun = nothing
        d.state_fun =  Function[]
        d.state_jac_fun = Function[]

        d.obj_eval = false
        d.cnstr_eval = false
        d.init_relax_run = false
        d.jacobian_sparsity = Tuple{Int64,Int64}[]

        d.np = N
        d.nx = 0
        d.ng = 0
        d.has_ineq = false

        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0

        d.var_relax = fill(zero(MC{N}),(N,))
        d.state_relax_1 = [zero(MC{N})]
        d.state_ref_relaxation_1 = fill(zero(MC{N}), (2,2,2,))
        d.state_relax_n = fill(zero(MC{N}), (2,2,))
        d.state_ref_relaxation_n = fill(zero(MC{N}), (2,2,2,))
        d.obj_relax = zero(MC{N})
        d.constraint_relax = [zero(MC{N})]
        d.P = fill(IntervalType(0.0), (1,))
        d.X = fill(IntervalType(0.0), (2,2))

        d.xa_mc = [zero(MC{N})]
        d.xA_mc = [zero(MC{N})]
        d.z_mc = [zero(MC{N})]
        d.aff_mc = [zero(MC{N})]
        d.H = [zero(MC{N})]
        d.J = fill(zero(MC{N}), (2,2))
        d.Y = fill(0.0, (2,2))

        d.last_p = Float64[0.0]
        d.ref_p = Float64[0.0]
        d.pref_mc = fill(zero(MC{N}),(N,))

        d.ivp = IVPInfo()
        d.IC_relaxations = [zero(MC{N})]

        d.disable_1storder = false
        d.disable_2ndorder = true
        d.has_nlobj = false

        d.current_node = NodeBB()
        d.last_node = NodeBB()

        d.opts = mc_opts()

        d.objective_ubd = Inf

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
function EAGO.build_evaluator!(d::ImplicitODELowerEvaluator, f::Function, h::Function, np::Int, nx::Int,
                               nt::Int, s::Int, t_start::Float64, t_end::Float64, method::Symbol,
                               pL::Vector{Float64}, pU::Vector{Float64}, xL::Vector{Float64}, xU::Vector{Float64}, ic::Function;
                               g = nothing, hj = nothing, state_update = x -> (), user_xtL = [], user_xtU = [], user_time = [])

    # setup objective and constraint functions
    d.has_nlobj = true
    d.objective_fun = f
    d.state_update = state_update
    d.initial_condition_fun = ic

    # sets up initial value problem parameters
    d.ivp.time_start = t_start
    d.ivp.time_end = t_end
    d.ivp.time_steps = nt
    d.ivp.step_size = (t_end - t_start)/(nt - 1)
    if length(user_time) <= 0
        d.ivp.time = zeros(nt)
        d.ivp.time[1] = t_start
        for i in 2:nt
            d.ivp.time[i] = d.ivp.time[i-1] + d.ivp.step_size
        end
    else
        d.ivp.time = user_time
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
    d.opts.nx = nx
    d.opts.np = np

    # preallocates the storage variables
    temp = zero(MC{np})

    # get dimension sizes
    d.np = np; d.nx = nx;
    if ~(g === nothing)
        d.has_ineq = true
        d.constraints_fun = g
        #ng = length(g([ones(nx) for j in 1:nt], ones(nx), ones(np), d.ivp.time))
        ng = 1
        d.ng = ng
        d.constraint_relax = fill(temp, (ng,))
    end

    d.var_relax = fill(temp, (np,))
    if (nx == 1)
        d.state_relax_1 = fill(temp, (nt-1,))
        d.state_ref_relaxation_1 = fill(temp, (1, d.opts.kmax, (nt-1)))
    else
        d.state_relax_n = fill(temp, (nx, nt-1))
        d.state_ref_relaxation_n = fill(temp, (nx, d.opts.kmax, (nt-1)))
    end

    d.pref_mc = fill(temp, (np,))

    # storage for intermediate calculations
    d.IC_relaxations = fill(temp, (nx,))
    d.xa_mc = fill(temp, (nx,))
    d.xA_mc = fill(temp, (nx,))
    d.z_mc = fill(temp, (nx,))
    d.aff_mc = fill(temp, (nx,))
    d.H = fill(temp, (nx,))
    if (nx > 1)
        d.J = zeros(MC{np}, nx, nx)
        d.Y = zeros(Float64, nx, nx)
    else
        d.J = zeros(MC{np}, 2, 2)
        d.Y = zeros(Float64, 2, 2)
    end

    d.P = fill(IntervalType(0.0), (np,))
    for i in 1:np
        d.P[i] = IntervalType(pL[i], pU[i])
    end

    d.X = fill(IntervalType(0.0), (nx, nt-1))
    indx = 1
    for i in 1:(nt-1)
        for j in 1:nx
            d.X[j,i] = IntervalType(xL[j], xU[j])
        end
        indx += 1
    end

    # allocates the reference points
    d.last_p = zeros(Float64, np); fill!(d.last_p, NaN)
    d.ref_p = zeros(Float64, np)
end

function EAGO.set_current_node!(x::ImplicitODELowerEvaluator,n::NodeBB)
    x.current_node = n
end

function EAGO.set_last_node!(x::ImplicitODELowerEvaluator,n::NodeBB)
    x.last_node = n
end

# LOOKS GREAT!
function relax_objective!(d::ImplicitODELowerEvaluator)
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
function relax_constraints!(d::ImplicitODELowerEvaluator)
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

num_state_variables(d::ImplicitODELowerEvaluator) = d.nx*(d.ivp.time_steps-1)
num_decision_variables(d::ImplicitODELowerEvaluator) = d.np

include("lower_calculations.jl")
include("lower_info.jl")
