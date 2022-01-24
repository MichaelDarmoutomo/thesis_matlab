%%
clear
load('data/x_opt.mat')
load('data/se.MAT')
se_struct = GetParameters(abs(se));

%% Settings of simulation
nSim = 1000;
T = 600;
Tp = 25;
gamma_list = [3,5,7];
rng(1194866)

%% Generate economy
E   = GenerateEconomySE(nSim,T, x_opt_struct, se_struct);
rho = 1 / (1 + mean(E.r,1:2));

%% Optimization
options = optimoptions('fmincon', 'Display', 'iter', 'TolFun', 1e-6, ...
    'TolX', 1e-4);

delete(gcp('nocreate'))
parpool('threads')

for i=1:2
    gamma = gamma_list(i);
    
    filename = sprintf('results/gamma_%d-nSim_%d_SE.txt', gamma, nSim);

    % Base Uniform 
    disp('Base uniform')
    fun = @(x) PensionFund(x, E, 1, nSim, T, gamma);
    [x_opt_uniform, fval_uniform, exitflag_uniform, output_uniform] = fmincon(...
        fun,[0.395; 0.815],[],[],[],[],[0,0],[1,1],[], options);

    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Base uniform:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_uniform);
    fprintf(fileID, 'fval = %.4f\n\n', fval_uniform);
    fclose(fileID);
    
    % 321
    disp('321')
    fun = @(x) PensionFund(x, E, 2, nSim, T, gamma);
    [x_opt_321, fval_321, exitflag_321, output_321] = fmincon(...
        fun,[0.3626, 0.8356],[],[],[],[],[0,0],[1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, '3-2-1 adjustment factor:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_321);
    fprintf(fileID, 'fval = %.4f\n\n', fval_321);
    fclose(fileID);
    
    % Uniform future pension
    disp('Uniform future pension')
    fun = @(x) PensionFund(x, E, 3, nSim, T, gamma);
    [x_opt_future, fval_future, exitflag_future, output_future] = fmincon(...
        fun,[0.4452, 0.9257],[],[],[],[],[0,0],[1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Uniform future pension:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_future);
    fprintf(fileID, 'fval = %.4f\n\n', fval_future);
    fclose(fileID);
    
    % Optimization over life cycle
    disp('Optimization over life cycle')
    fun = @(x) PensionFund(x, E, 4, nSim, T, gamma);
    [x_opt_life, fval_life, exitflag_life, output_life] = fmincon(...
        fun,[0.4512,0.8654,0.0003,  0.0011,  0.0017],[],[],[],[],...
        [0,0,0.0001,0.0001,0.0001],[1,1,1,1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Optimization over life cycle:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_life);
    fprintf(fileID, 'fval = %.4f\n\n', fval_life);
    fclose(fileID);

    % Optimization over the states
    disp('Optimization over the states')
    fun = @(x) PensionFund(x, E, 5, nSim, T, gamma);
    [x_opt_life, fval_life, exitflag_life, output_life] = fmincon(...
        fun, [0.395,0.815,0.003,0.004,0.005,0,0,0,0,0,0],...
        [],[],[],[],[0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1],...
        [1,1,1,1,1,1,1,1,1,1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Optimization over the states:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_life);
    fprintf(fileID, 'fval = %.4f\n\n', fval_life);
    fclose(fileID);
    
    savefile = sprintf('results/gamma_%d-nSim_%d.mat', gamma, nSim);
    save(savefile)

end
