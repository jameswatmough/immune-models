include("staged_treiv.jl")
# sample parameter sets
#
using Plots

function param_test(
    ;infection_rate = 4.9e-5
    ,resistance_rate = 4.9e-5
    ,viral_clearance_rate = 23.61
    ,time_to_full_budding = 5
    ,budding_death_rate = 0.35
    ,eclipse_death_rate = 0.
    ,num_eclipse_stages = 3
    ,num_budding_stages = 3
    ,max_budding_rate = 1000
    ,ramp_fun = n->[i/n for i in 1:n]
  )
  progression_rate = (num_eclipse_stages+num_budding_stages)/time_to_full_budding
  Param(
    infection_rate,
    resistance_rate,
    [progression_rate for i in 1:num_eclipse_stages],
    [eclipse_death_rate for i in 1:num_eclipse_stages], 
    [budding_death_rate for i in 1:num_budding_stages], 
    [progression_rate for i in 1:num_budding_stages-1],
    max_budding_rate*ramp_fun(num_budding_stages),
    viral_clearance_rate,
  )
end


# Ensemble Runs


## set up base ode problem

stages = 16
p = param_baseline(infection_rate = 0.0,ramp = rampfun(AlgRamp,m=0,n=stages))
initial_conditions = ComponentArray(  T = 1000.,  R = 0.,  I = [[1.]; [0 for i in 2:length(p.death_rate)]],  V = 0.)
prob = ODEProblem(staged_treiv_ode!,initial_conditions,[0,20],p)

plot_frame = plot( title="Budding Rate with algebraic ramp", xlabel="Time since cell infection", ylabel="Expected Budding Rate")

## loop through ramp-up functions

for i in 0:3:stages  
  global prob = remake(prob, p = param_baseline(infection_rate = 0.0, ramp = rampfun(AlgRamp,m=i,n=stages-i)) )
  sol = solve(prob)
  t = sol.t
  I = [sol[t].I[i] for t in 1:length(sol),  i in 1:length(sol[1].I) ]
  plot!(plot_frame, t, I*prob.p.budding_rate, label=string("E=",i,"; A=",15-i) )
end

# remake the ode problem for the simple single stage TIV model

p_TIV = param_baseline(infection_rate = 0.0, ramp = rampfun(AlgRamp,m=0,n=1)) 
u0_TIV = ComponentArray(  T = 1000.,  R = 0.,  I = [1.],  V = 0.)
prob = remake(prob,u0 = u0_TIV, p = p_TIV)
sol = solve(prob)
t = sol.t
I = [sol[t].I[i] for t in 1:length(sol),  i in 1:length(sol[1].I) ]

plot!(plot_frame, t, I*prob.p.budding_rate, label=string("E=0; A=0") )
              
plot_frame


