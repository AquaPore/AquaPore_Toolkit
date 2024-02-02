using Optim
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
lower= [0.0 ,0.0 ]
upper= [0.9, 1.0] 
initial_x = [-1.0, -1.0]
res = optimize(x->f(x),   lower, upper, NelderMead(), Optim.Options(g_tol = 1e-12,
                             iterations = 10,
                             store_trace = true,
                             show_trace = true))

X = Optim.minimizer(res)

println(X)

Of = Optim.minimum(res)

println(Of)