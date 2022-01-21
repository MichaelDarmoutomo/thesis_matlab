clear
load('data/x_opt.mat')
delete(gcp('nocreate'))

%% Settings of simulation
nSim = 10000;
T = 600;
Tp = 25;
gamma = 5;
rng(1194866)

%% Generate economy
E   = GenerateEconomy(nSim,T, x_opt_struct);
rho = 1 / (1 + mean(E.r,1:2));

%%
parpool('threads')
tic
PensionFund([0.5, 0.5], E, 1, nSim, T, gamma)
toc