######################### Running the Effective LV
# Set initial conditions
Ci = fill(0.0, N)
for i in 1:N
    Ci[i] = 0.1
end

tspan = (0.0, 15000.0)
prob_LV = ODEProblem(LV_dx!, Ci, t_span, p_lv)
sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
bm = sol_LV.u[length(sol_LV.t)][1:N]
length(bm[bm.>1e-7])

