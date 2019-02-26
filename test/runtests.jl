#!/usr/bin/env julia

#using Compat
using Compat.Test
using EAGO, EAGO_Differential, MathOptInterface
const MOI = MathOptInterface

@testset "Kernel and Jacobians" begin

    function h!(h, x, p, t)
        h[1] = x[1]*p[2] + t[1]*p[1]*x[1]^2 + x[2]^3
        h[2] = p[2]*x[1] + p[1]*x[2]^2
    end

    function hj!(J, x, p, t)
        J[1,1] = 2.0*t[1]*p[1]*x[1] + p[2]
        J[1,2] = 3.0*x[2]^2
        J[2,1] = p[2]
        J[2,2] = 2.0*p[1]*x[2]
    end

    nx = 2
    delT = 0.01
    x = [2.0; 1.0]
    xpast1 = [2.3; 1.1]
    xpast2 = [2.5; 1.7]
    p = [2.3; 6.7]
    tpast = [1.9]
    t = [2.0]

    hout = zeros(2)
    Jout = zeros(2,2)

    EAGO_Differential.bdf_kernel_1!(h!, hout, delT, x, xpast1, p, t)
    @test isapprox(hout[1], -0.62799999, atol=1E-6)
    @test isapprox(hout[2], -0.25700000, atol=1E-6)


    EAGO_Differential.bdf_jac_kernel_1!(hj!, Jout, delT, x, p, t)
    @test Jout[1,1] == 0.749
    @test Jout[1,2] == -0.03
    @test Jout[2,1] == -0.067
    @test Jout[2,2] == 0.954

    EAGO_Differential.bdf_kernel_2!(h!, hout, delT, x, xpast1, xpast2, p, t)
    @test isapprox(hout[1], -0.45199999999, atol = 1E-6)
    @test isapprox(hout[2], -0.0046666666, atol = 1E-6)

    EAGO_Differential.bdf_jac_kernel_2!(hj!, Jout, delT, x, p, t)
    @test isapprox(Jout[1,1], 0.83266666666, atol = 1E-6)
    @test Jout[1,2] == -0.02
    @test isapprox(Jout[2,1], -0.04466666666, atol = 1E-6)
    @test isapprox(Jout[2,2], 0.96933333333, atol = 1E-6)

    EAGO_Differential.am_kernel_1!(h!, hout, delT, x, xpast1, p, t)
    @test isapprox(hout[1], -0.627999999, atol = 1E-6)
    @test isapprox(hout[2], -0.257000000, atol = 1E-6)

    EAGO_Differential.am_jac_kernel_1!(hj!, Jout, delT, x, p, t)
    @test Jout[1,1] == 0.749
    @test Jout[1,2] == -0.03
    @test Jout[2,1] == -0.067
    @test Jout[2,2] == 0.954

    EAGO_Differential.am_kernel_2!(h!, hout, delT, x, xpast1, p, t, tpast)
    @test isapprox(hout[1], -0.6632914999, atol = 1E-6)
    @test isapprox(hout[2], -0.2694650000, atol = 1E-6)

    EAGO_Differential.am_jac_kernel_2!(hj!, Jout, delT, x, p ,t)
    @test Jout[1,1] == 0.8745
    @test Jout[1,2] == -0.015
    @test Jout[2,1] == -0.0335
    @test Jout[2,2] == 0.977
end

include("lower_evaluator_tests.jl")
include("upper_evaluator_tests.jl")
include("core_optimizer_test.jl")
