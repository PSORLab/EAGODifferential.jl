tau = 20.0
k2 = 0.4

pL = [0.9; 0.8; 10.0]
pU = [1.1; 1.0; 50.0]

# Need to determine upper bounds on xB = X, xD = Y
#
X = 1.0
Y = 1.0
xL = [0.0; 0.0; 0.0; 0.0]
xU = [0.2; X; 0.5; Y]

f(out, x, x0, p, t) = ()
g(out, x, x0, p, t) = ()

function h(out, x, p, t)
    out[1] = -p[3]*x[1]*x[2] - k2*x[1]*x[3] + (p[1] - 2.0*x[1])/tau
    out[2] = -p[3]*x[1]*x[2] + (p[2] - 2.0*x[2])/tau
    out[3] = p[3]*x[1]*x[2] - k2*x[1]*x[3] - (2.0/tau)*x[3]
    out[4] = k2*x[1]*x[3]  - (2.0/tau)*x[4]
end

function hj(out, x, p, t)

    fill!(out, zero(p[1]))

    out[1,1] = -p[3]*x[2] - k2*x[3] - 2.0/tau
    out[1,2] = -p[3]*x[1]
    out[1,3] =  -k2*x[1]

    out[2,1] = -p[3]*x[2]
    out[2,2] = -p[3]*x[1] - 2.0/tau

    out[3,1] = p[3]*x[2] - k2*x[3]
    out[3,2] = p[3]*x[1]
    out[3,3] =  -k2*x[1] - (2.0/tau)

    out[4,1] = k2*x[3]
    out[4,3] = k2*x[1]
    out[4,4] = -(2.0/tau)
end

x0(p) = [0.0; 0.0; 0.0; 0.0]

time_start = 0.0
time_end = 10.0
nx = 4
np = 3
nt1 = 21

lower_vars = Float64[]
upper_vars = Float64[]
for i in 1:(nt1-1)
    append!(lower_vars, xL)
    append!(upper_vars, xU)
end
append!(lower_vars, pL)
append!(upper_vars, pU)
n = NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)
y = [xL]; (pL+pU)/2]

lower_eval1 = ImplicitODELowerEvaluator{np}()
build_evaluator!(lower_eval1, f, h, np, nx, nt1, 1, t_start, t_end, :AM, pL, pU, xL, xU, x0; g = g, hj = hj)
EAGO.set_current_node!(lower_eval1, n)
EAGO_Differential.relax_ode_implicit!(lower_eval1, y)

lower_eval2 = ImplicitODELowerEvaluator{np}()
build_evaluator!(lower_eval2, f, h, np, nx, nt, 2, t_start, t_end, :AM, pL, pU, xL, xU, x0; g = g, hj = hj)
EAGO.set_current_node!(lower_eval2, n)
EAGO_Differential.relax_ode_implicit!(lower_eval2, y)

lower_eval3 = ImplicitODELowerEvaluator{np}()
build_evaluator!(lower_eval3, f, h, np, nx, nt, 2, t_start, t_end, :BDF, pL, pU, xL, xU, x0; g = g, hj = hj)
EAGO.set_current_node!(lower_eval3, n)
EAGO_Differential.relax_ode_implicit!(lower_eval3, y)
