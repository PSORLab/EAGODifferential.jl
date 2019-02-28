mutable struct IVPInfo
    time_start::Float64
    time_end::Float64
    time_steps::Int
    step_size::Float64
    time::Vector{Float64}
    method::Symbol
    method_order::Int
end
IVPInfo() = IVPInfo(0.0, 0.0, 0, 0.0, Float64[], :bdf, 2)

# Defines kernel operator for BDF integrator (1st degree, implicit euler)
@inline function bdf_kernel_1!(h!::Function, hout, delT::Float64, xout, x1, p, tout)
    h!(hout, xout, p, tout)
    hout[:] *= -delT
    hout[:] += xout
    hout[:] -= x1
end

@inline function bdf_jac_kernel_1!(hj!::Function, Jout, delT::Float64, xout, p, tout)
    #println("delT: $delT")
    #println("xout: $xout")
    #println("p: $p")
    #println("tout: $tout")
    hj!(Jout, xout, p, tout)
    Jout[:,:] *= -delT
    Jout[:,:] += I
    #println("Jout: $Jout")
end

# second order BDF kernel
@inline function bdf_kernel_2!(h!::Function, hout, delT::Float64, xout, x1, x2, p, tout)
    h!(hout, xout, p, tout)
    hout[:] *= -delT*2.0/3.0
    hout[:] += xout
    hout[:] -= (4.0/3.0)*x1
    hout[:] += (1.0/3.0)*x2
end

@inline function bdf_jac_kernel_2!(hj!::Function, Jout, delT::Float64, xout, p, tout)
    hj!(Jout, xout, p, tout)
    Jout[:,:] *= -delT*2.0/3.0
    Jout[:,:] += I
end

# Defines kernel operator for AM integrator (1st degree, implicit euler)
@inline function am_kernel_1!(h!::Function, hout, delT::Float64, xout, x1, p, tout)
    h!(hout, xout, p, tout)
    hout[:] *= -delT
    hout[:] += xout
    hout[:] -= x1
end

@inline function am_jac_kernel_1!(hj!::Function, Jout, delT::Float64, xout, p, tout)
    hj!(Jout, xout, p, tout)
    Jout[:,:] *= -delT
    Jout[:,:] += I
end

# Defines kernel operator for AM integrator (2nd degree, Trapezoidal)
@inline function am_kernel_2!(h!::Function, hout, delT::Float64, xout, x1, p, tout, t1)
    htemp = copy(hout)
    h!(hout, xout, p, tout)
    h!(htemp, x1, p, t1)
    hout[:] += htemp
    hout[:] *= -0.5*delT
    hout[:] += xout
    hout[:] -= x1
end

@inline function am_jac_kernel_2!(hj!::Function, Jout, delT::Float64, xout, p, tout)
    hj!(Jout, xout, p, tout)
    Jout[:,:] *= - 0.5*delT
    Jout[:,:] += I
end
