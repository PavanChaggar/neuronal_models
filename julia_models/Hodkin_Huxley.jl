using Pkg
Pkg.status()
using Plots
using DifferentialEquations

#----------------------------------------------
# Set Parameters:
# Membrane capacitance and potential
C_m = 1.0
V_m = -65.0

# Max sodium conductance and membrane potential
g_Na = 120.0
V_Na = 56.0

# Max potassium conductance and membrane potential
g_K = 36.0
V_K = -77.0

# Max leak conductance and membrane potential
g_l = 0.3 
V_l = -60.0

# Time to integrate over
time = 0.0:0.1:1000

# ion variables
m = 0.05
h = 0.6
n = 0.32

#----------------------------------------------
# Potassium Rate Functoins
a_n(V_m) =  0.01 * (V_m + 55) / (1 - exp(-(V_m + 55.0) / 10))

B_n(V_m) = 0.125 * exp(-(V_m + 65.0) / 80)


#----------------------------------------------
# Sodium Rate Functions
a_m(V_m) = 0.1 * (40 + V_m) / (1.0 - exp(-(V_m + 40.0) / 10.0))

a_h(V_m) = 0.07 * exp(-(V_m + 65.0) / 20.0)

B_m(V_m) = 4 * exp(-(V_m + 65.0) / 18.0)

B_h(V_m) = 1 / (exp(-(35 + V_m) / 10) + 1)


#----------------------------------------------
# Ion Currents
IK(g_K, n, V_m, V_K) = g_K * n^(4) * (V_m - V_K)

INa(g_Na, m, h, V_m, V_Na) =  g_Na * m^(3) * h * (V_m - V_Na)

LI(g_l, V_m, V_l) = g_l * (V_m - V_l)

#----------------------------------------------
# Input current
I_Inj = 20

input(t::AbstractFloat) = ifelse(0 < t < 200, zero(t), ifelse(200<t<300, one(t)*I_Inj, zero(t)))

#----------------------------------------------
# Hodking Huxley model
function hodgkin_huxley!(du, u, p, t)
    # dVdt, dndt, dmdt, dhdt = u
    #
    #
    V_m, n, m, h = u
    g_K, V_K, g_Na, V_Na, g_l, C_m = p
    
    # Total current through the membrane
    du[1] = input(t) - IK(g_K,n,V_m,V_K) - INa(g_Na, m, h, V_m, V_Na) - LI(g_l, V_m, V_l) / C_m

    # Derivative of n, potassium channel activation, w.r.t. time
    du[2] = a_n(V_m) * (1 - n) - B_n(V_m) * n

    # Derivative of m, sodium channel activion, w.r.t. time
    du[3] = a_m(V_m) * (1 - m) - B_m(V_m) * m

    # Derivative of h, leaky channel in-activion, w.r.t. time

    du[4] = a_h(V_m) * (1 - h) - B_h(V_m) * h

end

#----------------------------------------------
# Define and solve ODE problem
p = [g_K, V_K, g_Na, V_Na, g_l, C_m]
u0 = [V_m, n, m, h]
tspan = (0.0,500.0)

prob = ODEProblem(hodgkin_huxley!,u0,tspan,p)
sol = solve(prob)

#----------------------------------------------
# Plotting Solutions
plot(sol, vars=1)
#animate(sol, vars=(1), every=20)

plot(sol, vars=(2,3))

plot(sol, vars=(2,3,4))
