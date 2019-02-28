using EAGO_Differential, EAGO, DataFrames

# soft build, then evaluate
lower_eval = ImplicitODELowerEvaluator{1}()

f(x,x0,p,t) = x[1]
function h(H,x,p,t)
    H[1] = p[1]*x[1]*(1-x[1])
end
function hj(J,x,p,t)
    J[1,1] = p[1]*(1.0-2.0*x[1])
end
x0(p) = [0.1]

np = 1
nx = 1
nt = 10
s = 1

t_start = 0.0
t_end = 0.333333333
method = :BDF

pL = [1.2]; pU = [1.4]
xL = [0.1]; xU = [0.4]

# build the basic evaluator (w/o inequality constraints)
EAGO_Differential.build_evaluator!(lower_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)

lower_vars = fill(xL[1], (nt-1,))
upper_vars = fill(xU[1], (nt-1,))
append!(lower_vars, pL)
append!(upper_vars, pU)
n = NodeBB(lower_vars, upper_vars, -Inf, Inf, 0, -1, false)

y = [0.2*ones(9); 1.3]
EAGO.set_current_node!(lower_eval, n)
EAGO_Differential.relax_ode_implicit!(lower_eval, y)

println(" Eval Box ")
println("       ")
upper_eval = ImplicitODEUpperEvaluator()
EAGO_Differential.build_evaluator!(upper_eval, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)
EAGO.set_current_node!(upper_eval, n)
EAGO_Differential.relax_ode_implicit!(upper_eval)

println("       ")
println(" Eval Point ")
println("       ")
upper_eval1 = ImplicitODEUpperEvaluator()
EAGO_Differential.build_evaluator!(upper_eval1, f, h, np, nx, nt, s, t_start, t_end, method, pL, pU, xL, xU, x0; hj = hj)
n_new = NodeBB([lo.(upper_eval.state_relax_1); 1.3], [hi.(upper_eval.state_relax_1); 1.3], -Inf, Inf, 0, -1, false)
EAGO.set_current_node!(upper_eval1, n_new)
EAGO_Differential.relax_ode_implicit!(upper_eval1)

state_relax_1_lo = lo.(upper_eval1.state_relax_1)
state_relax_1_hi = hi.(upper_eval1.state_relax_1)
