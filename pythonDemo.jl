using SciPy

jlFunc(x) = SciPy.optimize.rosen(x)
# jlFunc = SciPy.optimize.rosen

x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
res = SciPy.optimize.minimize(jlFunc, x0, method="Nelder-Mead", tol=1e-6)
println(res)