function dI = emg_exp1fit(beta,dt)
dI=beta(1)*exp(-dt/beta(2));