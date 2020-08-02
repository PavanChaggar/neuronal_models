using Plots
using DifferentialEquations
using DiffEqUncertainty
#----------------------------------------------
# Specify and Solve model:

function fitzhugh_nagumo!(du, u, p, t)

    v, w = u
    a, b, I, R, τ = p


    du[1] = v - (v^3)/3 - w  + R * I

    du[2] = (v + a - (b * w)) / τ
end

u0 = [-2.0, 1.0]
p = [0.7, 0.8, 0.34, 1.0, 10]
tspan = (0.0, 100.0)
prob = ODEProblem(fitzhugh_nagumo!,u0,tspan,p)
sol = solve(prob)

plot(sol, vars =(1))


cb = ProbIntsUncertainty(0.1,1)
cb2 = ProbIntsUncertainty(0.05,1)
ensemble_prob = EnsembleProblem(prob)
sim = solve(ensemble_prob,Euler(),trajectories=100,callback=cb,dt=1/10)
sim2 = solve(ensemble_prob,Euler(),trajectories=50,callback=cb2,dt=1/100)

plot(sim, vars=(0,1), color="grey")
plot!(sim2, vars=(0,1), color="lightsalmon2")
plot!(sol, vars =(1), color="red3", linewidth =3)
