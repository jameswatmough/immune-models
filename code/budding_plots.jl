# produces two plots 'plot_arith' and 'plot_linear'
# for runs with arithmetic and linear ramp up


include("staged_treiv.jl")

using Plots

# parameter sets
#

function param_budding(
    ;infection_rate = 0.0
    ,viral_clearance_rate = 23.61
    ,eclipse_duration = 5
    ,max_death_rate = 0.35
    ,max_budding_rate = 1000
    ,ramp = rampfun(ArithRamp,m=2,n=3)
  )
  stages = length(ramp)
  prog_rate = [stages/eclipse_duration for i in 1:stages]
  death_rate = max_death_rate*ramp
  budding_rate = max_budding_rate*ramp
  Param(
    infection_rate,
    prog_rate,
    death_rate,
    budding_rate,
    viral_clearance_rate,
  )
end

# Ensemble Runs

## Arithmetic Ramp up compared with TIV model

### set up base ode problem

stages = 16 
p = param_budding(ramp = rampfun(ArithRamp,m=0,n=stages-1))
initial_conditions = ComponentArray(
  T = 1000.,
  I = [[1.]; [0 for i in 2:length(p.death_rate)]],
  V = 0.
)
prob = ODEProblem(staged_treiv_ode!,initial_conditions,[0,20],p)

plot_arith = plot( title="Budding Rate with arithmetic ramp", xlabel="Time since cell infection", ylabel="Expected Budding Rate")

### loop through ramp-up functions

for m in 0:3:(stages-1)
  n = stages-1-m
  global prob = remake(prob, p = param_budding(ramp = rampfun(ArithRamp,m=m,n=n)))
  sol = solve(prob)
  t = sol.t
  I = [sol[t].I[i] for t in 1:length(sol),  i in 1:length(sol[1].I) ]
  plot!(plot_arith, t, I*prob.p.budding_rate, label=string("E=",m,"; A=",n) )
end

## remake the ode problem for the simple single stage TIV model

let p,u0,sol,I
  p = param_budding(ramp = rampfun(LinearRamp,start=1.0,n=0)) 
  u0 = ComponentArray( T = 1000., I = [1.],  V = 0. )
  TIVprob = remake(prob,u0 = u0, p = p)
  sol = solve(TIVprob)
  I = [sol[t].I[i] for t in 1:length(sol),  i in 1:length(sol[1].I) ]

  plot!(plot_arith, sol.t, I*TIVprob.p.budding_rate, label=string("E=0; A=0") )
end
              
plot_arith

## Linear Ramp up 

### set up base ode problem


stages = 16
p = param_budding(ramp = rampfun(LinearRamp,start=0,n=stages-1))
initial_conditions = ComponentArray(
  T = 1000.,
  I = [[1.]; [0 for i in 2:length(p.death_rate)]],
  V = 0.
 )
prob = ODEProblem(staged_treiv_ode!,initial_conditions,[0,20],p)

plot_linear = plot( title="Budding Rate with Linear Ramp", xlabel="Time since cell infection", ylabel="Expected Budding Rate")

### loop through ramp-up functions

for start in 0:.2:1  
  global prob = remake(prob, p = param_budding(ramp = rampfun(LinearRamp,start=start,n=stages-1)) )
  sol = solve(prob)
  t = sol.t
  I = [sol[t].I[i] for t in 1:length(sol),  i in 1:length(sol[1].I) ]
  plot!(plot_linear, t, I*prob.p.budding_rate, label=string("ramp start =",start) )
end

plot_linear
