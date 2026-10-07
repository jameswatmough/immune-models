# Staged TEIV model
# really a staged-TIV model with 
# m+n+1 stages
# where m stages have zero budding (true eclipse) 
# n stages with budding ramp up
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

# Axes constructors for the state and event component arrays

"""
  get_state_axis(stages)

Return an Axis for state variables
"""
function get_state_axis(stages)
  return Axis(T=1,V=2,I=3:stages+2) 
end

"""
  get_event_axis(stages)

Return an Axis for transition events
"""
function get_event_axis(stages)
  return Axis(incidence = 1,
              virus_decay = 2,
              budding = 3,
              progression = 4:2+stages,
              cell_death = 3+stages:2+2*stages
             )
end

# Model parameters used in general ode 

@with_kw struct Param
  infection_rate::Float64           = 0.000049
  progression_rate::Vector{Float64} = [5/5,5/5,5/5,5/5,5/5] # m+n=5
  death_rate::Vector{Float64}       = [0,0,0,.35,.35,.35]   # m+n+1=6
  budding_rate::Vector{Float64}     = [0,0,0,.1,.2,.3]
  viral_clearance_rate::Float64     = 23.61
  state_axis                        = get_state_axis(6)
  event_axis                        = get_event_axis(6) 
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

"""
  rampfun(ArithRamp[, m, n]) 

  return budding rates for arithmetic ramp up with m+n stages with m zeros
  returns a vector of length m+n+1
"""
rampfun(::Type{ArithRamp};m=1,n=7) = [[0. for i in 1:m]; [i/(n+1) for i in 1:(n+1)]]

"""
  rampfun(LinearRamp[, start, n]) 

  return budding rates for linear ramp up from start to 1
  returns a vector of length n+1

  This is intended as a bridge between a TEIV model with an eclipse stage
  and the TIV model with no eclipse stages.
  Budding starts immediately after infection, but peaks some time after infection.
"""
rampfun(::Type{LinearRamp};start=-0.5,n=15) = 
  n==0 ?
    [1.0] :
    [max(0.0,1.0 + (start-1.0)*i/n) for i in n:-1:0]

"""
  rampfun(GeomRamp[, m, n, rampfactor) 

  return budding rates for geometric ramp up with m+n stages with m zeros
  and n+1 entries starting from rampfactor and ending at one

  This limits to a discrete delay as rampfactor -> 0
"""
rampfun(::Type{GeomRamp};m=1,n=7,rampfactor=0.1) = [[0. for i in 1:m]; [rampfactor^(1-i/n) for i in 0:n]]

# constructors for parameter sets

"""
  param_baseline([infection_rate, viral_clearance_rate, eclipse_duration, max_death_rate, max_budding_rate, ramp])

  returns a default set of parameter values

"""
function param_baseline(
    ;infection_rate = 4.9e-5
    ,viral_clearance_rate = 23.61
    ,eclipse_duration = 5
    ,max_death_rate = 0.35
    ,max_budding_rate = 1000
    ,ramp = rampfun(AlgRamp,m=2,n=3)
  )
  stages = length(ramp)
  prog_rate = [stages/eclipse_duration for i in 1:stages-1]
  death_rate = max_death_rate*ramp
  budding_rate = max_budding_rate*ramp
  state_axis = get_state_axis(stages)
  event_axis = get_event_axis(stages)
  Param(
    infection_rate,
    prog_rate,
    death_rate,
    budding_rate,
    viral_clearance_rate,
    state_axis,
    event_axis
  )
end

"""
  ic_baseline(p)

  returns a default set of initial conditions consistent with parameter values in p

  I.e., a ComponentVector based on p.state_axis

  # Arguments
  
  - `p` structure of parameters (Param)

"""
ic_baseline = function(p)
  initial_conditions = ComponentVector(zeros(lastindex(p.state_axis)),p.state_axis)
  initial_conditions.T = 1000.
  initial_conditions.V = 100.
  return(initial_conditions)
end

"""
  staged_tiv_ode!(dx,x,p,t)

  returns the vector field for the staged tiv model

  for use in constructing ODEProblems or SDEProblems

  # Arguments
  
  `dx`  the vector field (a.k.a 'drift' in the SDE model
  `x` state vector; a ComponentVector with axis p.state_axis
  `p` a Parameter structure (Param)
  `t` time
"""
function staged_tiv_ode!(dx,x,p,t)

  # get rid of any negative entries
  @. x = max(x,0.)  

  @unpack T, V, I = x

  incidence = p.infection_rate.*T.*V
  budding = sum(p.budding_rate.*I)

  dx.I = -p.death_rate.*I
  dx.I[1:end-1] .-= p.progression_rate.*I[1:end-1] 
  dx.I[2:end] .+= p.progression_rate.*I[1:end-1] 
  

  dx.T = -incidence
  dx.I[1] += incidence

  dx.V = budding .- p.viral_clearance_rate.*V .- incidence

  return(dx)

end

# SDE model

import Random, SparseArrays

# use a ComponentVector for constructing the covariance matrix
# rows of G correspond to the state variables 
# columns of G correspond to events


"""
  Gpattern(p)

  returns the sparseArray pattern for the state transition matrix passed to SDEProblem

  # Argument
  
  `p` a Parameter structure (Param)
"""
function Gpattern(p)
  # construct sparse array pattern for G, the root of the covariance matrix 

  # construct ComponentArray and sent all entries to zero
  G = ComponentMatrix( zeros(lastindex(p.state_axis),lastindex(p.event_axis)), (p.state_axis, p.event_axis))

  # set entries corresponding to transitions to 1
  stages = length(p.death_rate) # for steping along diagonals of G
  G[:T,:incidence]=1
  G[:V,:incidence]=1
  @view(G[:I,:incidence])[1]=1
  G[:V,:virus_decay] = 1
  @view(G[:I,:budding])[1] = 1
  @view(G[:I,:progression])[1:(stages+1):end] .= 1  # diagonals
  @view(G[:I,:progression])[2:(stages+1):end] .= 1  # subdiagonals
  @view(G[:I,:cell_death])[1:(stages+1):end] .= 1   # diagonals

  G = SparseArrays.sparse(G)

end


"""
  staged_tiv_noise!(dx,x,p,t)

  returns the root of the CoVariance Matrix the staged tiv model

  for use in constructing the SDEProblem

  # Arguments
  
  `dx`  the vector field (a.k.a 'drift' in the SDE model
  `x` state vector; a ComponentVector with axis p.state_axis
  `p` a Parameter structure (Param)
  `t` time
"""
function staged_tiv_noise!(dx_raw,x,p,t)


  # get rid of any negative entries
  @. x = max(x,0.)  

  @unpack T, V, I = x

  # embellish dx with row and column names
  dx = ComponentMatrix(dx_raw, (p.state_axis,p.event_axis))
  
  # target cell infection: T-=1 ; V-=1; I[1] +=1
  # target cell resistence: T-=1 ; R+=1
  # virion degradation: V-=1
  # virus production: V+=1
  # infected cell progression: I[i]-=1; I[i+1]+-1 for i in 1:(m+n)
  # infected cell death: I[i]-=1 for i in 1:(m+n+1)
  # p.budding_rate has length m+n+1; the total number of infected stages
  #
  # dx is a ComponentArray with columns indexed by events

  σ_incidence = sqrt(p.infection_rate.*T.*V)
  dx[:T,:incidence] = -σ_incidence
  dx[:V,:incidence] = -σ_incidence 
  @view(dx[:I,:incidence])[1] = σ_incidence 

  dx[:V,:budding] = sqrt(sum(p.budding_rate.*I))
  dx[:V,:virus_decay] = -sqrt(sum(p.viral_clearance_rate.*I))

  # there are m+n different independent processes for progress
  stages = length(I)
  dIprog = @view(dx[:I,:progression])
  σ_prog = sqrt.(p.progression_rate.*I[1:end-1] )
  dIprog[1:(stages+1):end] .= -σ_prog
  dIprog[2:(stages+1):end] .= σ_prog

  # there are m+n+1 different independent processes for progress
  dIdeath = @view(dx[:I,:cell_death])
  dIdeath[1:stages+1:end] = -sqrt.(p.death_rate.*I)
  
  return(dx)

end
