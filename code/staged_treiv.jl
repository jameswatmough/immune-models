# Staged TREIV model
# really a staged-TRIV model with 
# m+n+1 stages
# where m stages have zero budding (true eclipse) 
# n stages with budding ramp up
# resistant (R) compartment
# should possibly include interferon, but keep resistance simple for first go
#
# Note there is a potential for confusion of the number of stages since a given model has eclipse stages, ramp-up stages, and infectious stages
# The convention taken here is that 
#     there are n + m + 1 stages
#     m refers to number of eclipse stages (zero budding)
#     n refers to number of ramp-up stages
#     there is one additional 'fully-productive' infected stage
# The classic TIV model has m=n=0      ( the single stage model)
# The classic TEIV model has m=1, n=0  ( a two stage model)

using DifferentialEquations, ComponentArrays, Parameters

# Model parameters used in general ode 

@with_kw struct Param
  infection_rate::Float64           = 0.000049
  resistance_rate::Float64          = 0.000049
  progression_rate::Vector{Float64} = [5/5,5/5,5/5,5/5,5/5] # m+n=5
  death_rate::Vector{Float64}       = [0,0,0,.35,.35,.35]
  budding_rate::Vector{Float64}     = [0,0,0,.1,.2,.3]
  viral_clearance_rate::Float64     = 23.61
#  @assert length(budding_death_rate) == length(budding_rate)
#  @assert length(budding_progression_rate) == length(budding_rate)
#  @assert length(eclipse_progression_rate) == length(eclipse_death_rate)
end

# ramp up functions
#
# for consistency the ramp functions should return vectors of length m+n+1

abstract type RampType end
struct GeomRamp <: RampType end
struct ArithRamp  <: RampType end
struct LinearRamp  <: RampType end

""" rampfun(::Type{ArithRamp}) 
    return budding rates for arithmetic ramp up with m+n stages with m zeros
    returns a vector of length m+n+1"""
    rampfun(::Type{ArithRamp};m=1,n=7) = [[0. for i in 1:m]; [i/(n+1) for i in 1:(n+1)]]

""" rampfun(::Type{LinearRamp}) 
    return budding rates for linear ramp up from start to 1
    returns a vector of length n+1"""
    rampfun(::Type{LinearRamp};start=-0.5,n=15) = 
        n==0 ?
          [1.0] :
          [max(0.0,1.0 + (start-1.0)*i/n) for i in n:-1:0]

""" rampfun(::Type{GeomRamp};m=1; n=7,rampfactor=0.1) 
    return budding rates for geometric ramp up with m+n stages with m zeros
    and n+1 entries starting from rampfactor and ending at one """
    rampfun(::Type{GeomRamp};m=1,n=7,rampfactor=0.1) = [[0. for i in 1:m]; [rampfactor^(1-i/n) for i in 0:n]]

# constructors for parameter sets

function param_baseline(
    ;infection_rate = 4.9e-5
    ,resistance_rate = 4.9e-5
    ,viral_clearance_rate = 23.61
    ,eclipse_duration = 5
    ,max_death_rate = 0.35
    ,max_budding_rate = 1000
    ,ramp = rampfun(AlgRamp,m=2,n=3)
  )
  stages = length(ramp)
  prog_rate = [stages/eclipse_duration for i in 1:stages]
  death_rate = max_death_rate*ramp
  budding_rate = max_budding_rate*ramp
  Param(
    infection_rate,
    resistance_rate,
    prog_rate,
    death_rate,
    budding_rate,
    viral_clearance_rate,
  )
end

# sample parameter sets
p = param_baseline(
      ramp = rampfun(GeomRamp,m=7,n=8,rampfactor=.2)
    )
        
initial_conditions = ComponentArray(
  T = 1000.,
  R = 0.,
  I = [0 for i in 1:length(p.budding_rate)],
  V = 100.
)

function staged_treiv_ode!(dx,x,p,t)

  @unpack T, R, I, V = x

  incidence = p.infection_rate.*T.*V
  budding = sum(p.budding_rate.*I)

  dx.I = -p.death_rate.*I
  dx.I[1:end-1] .-= p.progression_rate[1:end-1].*I[1:end-1] 
  dx.I[2:end] .+= p.progression_rate[1:end-1].*I[1:end-1] 
  

  dx.R = p.resistance_rate*sum(I)*T
  dx.T = -incidence - p.resistance_rate*sum(I)*T
  dx.I[1] += incidence

  dx.V = budding .- p.viral_clearance_rate.*V .- incidence

  return(dx)

end

