using Plots
using DifferentialEquations

#----------------------------------------------
# Set Parameters:

v = -2.0
w = 1.0

a = 0.7
b = 0.8
τ = 10
I = 0.35


tspan = (0.0, 100.0)

#----------------------------------------------
# Specify and Solve model:

function fitzhugh_nagumo!(du, u, p, t)

    v, w = u
    a, b, I = p


    du[1] = v - w  * I - (w^3)/3

    du[2] = (v + a - (b * w)) / τ
end

u0 = [v, w]
p = [a, b, I]
prob = ODEProblem(fitzhugh_nagumo!,u0,tspan,p)
sol = solve(prob)

plot(sol, vars =1)