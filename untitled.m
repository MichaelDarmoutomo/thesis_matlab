clear
load('data/x_opt_restricted.mat')
load('data/se.MAT')
se_struct = GetParameters(abs(se));
delete(gcp('nocreate'))

%% Settings of simulation
nSim = 1000;
T = 600;
Tp = 25;
gamma = 5;
rng(1194866)

%% Generate economy
tic
E   = GenerateEconomySE(nSim,T, x_opt_struct, se_struct);
toc
rho = 1 / (1 + mean(E.r,1:2));

%%
parpool('threads')
tic
PensionFund([0.5, 0.5], E, 1, nSim, T, gamma)
toc