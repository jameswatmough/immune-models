include("staged_tiv.jl")



# test parameter set
p = param_baseline( ramp = rampfun(ArithRamp,m=1,n=2))
        
u0 = ic_baseline(p)

dx_raw = zeros(lastindex(p.state_axis),lastindex(p.event_axis))
dx = ComponentMatrix(dx_raw,(p.state_axis,p.event_axis))

G = Gpattern(p)

tspan = (0.,30.)
prob = SDEProblem(staged_tiv_ode!,staged_tiv_noise!,u0,tspan,p, noise_rate_prototype = G)

sol = solve(prob,EM(); dt=0.01)


I = [max(0.,sol[t].I[i]) for t in 1:length(sol),  i in 1:length(sol[1].I) ]
V = [max(0.,sol[t].V) for t in 1:length(sol)]

I_plot = plot(sol.t,I,title = "Infected Cells by stage of infection", ylabel = "Cells (I)", xlabel = "time (days) since infection"

V_plot = plot(sol.t,V,ylabel = "Viral Load (V)", xlabel = "time (days) since infection"

