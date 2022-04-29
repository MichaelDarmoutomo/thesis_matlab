%% load data
clear

restricted = true;
se = true;

if (restricted)
   file = 'data/x_opt_restricted.mat';
else
   file = 'data/x_opt.mat';
end

load(file)

if (se)
   load('data/se_restricted.MAT')
   load('data/covariance_restricted.mat')
   se_struct = GetParameters((se_restricted./10)', restricted);
end

%% Settings of simulation
nSim = 1000;
T = 600;
Tp = 25;
gamma_list = [5];
rng(1194866)

%% Generate economy

% if (se)
%     E = GenerateEconomySE(nSim,T,  x_opt_struct,se_struct, x_opt, covariance_restricted);
% else
%     E_basic = GenerateEconomy(nSim, T, x_opt_struct);
% end
% 
% save("data/economy_cov.mat", "E", "-v7.3")
load('data/economy_cov_0275.mat')
E.delta_0pi = E.delta_0pi(:,1:nSim);
E.delta_1pi = E.delta_1pi(:,1:nSim);
E.delta_0r = E.delta_0r(:,1:nSim);
E.delta_1r = E.delta_1r(:,1:nSim);
E.K = E.K(:,:,1:nSim);
E.dZ = E.dZ(:,1:nSim,:);
E.Lambda = E.Lambda(:,1:nSim,:);
E.r = E.r(:,1:nSim);
E.pi = E.pi(:,1:nSim);
E.dPi_ = E.dPi_(:,1:nSim);
E.Pi = E.Pi(:,1:nSim);
E.dS_ = E.dS_(:,1:nSim);
E.S = E.S(:,1:nSim);
E.dX = E.dX(:,1:nSim,:);
E.X = E.X(:,1:nSim,:);
E.P = E.P(:,1:nSim,:);
E.P2 = E.P2(:,1:nSim,:);
E.w = E.w(:,1:nSim);

%% Optimization
options = optimoptions('fmincon', 'Display', 'iter', 'TolFun', 1e-6, ...
    'TolX', 1e-4);

delete(gcp('nocreate'))
parpool('threads')

for i=1:length(gamma_list)
    gamma = gamma_list(i);
    for j=10
    filename = sprintf('results/gamma_%d-nSim_%d_cov_0275.txt', gamma, nSim);

%     Base Uniform 
    disp('Base uniform')    
    fun = @(x) PensionFundSE(x, E, 1, nSim, T, gamma, j);
    [x_opt_uniform, fval_uniform, exitflag_uniform, output_uniform] = fmincon(...
        fun,[0.45; 0.95],[],[],[],[],[0,0],[1,1],[], options);

    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Base uniform:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_uniform);
    fprintf(fileID, 'fval = %.4f\n\n', fval_uniform);
    fclose(fileID);
    
%     321
    disp('321')
    fun = @(x) PensionFundSE(x, E, 2, nSim, T, gamma, j);
    [x_opt_321, fval_321, exitflag_321, output_321] = fmincon(...
        fun,[0.4995, 0.9524],[],[],[],[],[0,0],[1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, '3-2-1 adjustment factor:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_321);
    fprintf(fileID, 'fval = %.4f\n\n', fval_321);
    fclose(fileID);
    
    % Uniform future pension
    disp('Uniform future pension')
    fun = @(x) PensionFundSE(x, E, 3, nSim, T, gamma,j);
    [x_opt_future, fval_future, exitflag_future, output_future] = fmincon(...
        fun,[0.4, 0.95],[],[],[],[],[0,0],[1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Uniform future pension (x=%d):\n', j);
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_future);
    fprintf(fileID, 'fval = %.4f\n\n', fval_future);
    fclose(fileID);
    
%     Optimization over life cycle
    disp('Optimization over life cycle')
    fun = @(x) PensionFundSE(x, E, 4, nSim, T, gamma, j);
    [x_opt_life, fval_life, exitflag_life, output_life] = fmincon(...
        fun,[0.39,0.99, -0.0081, 0.0004, 0.0017],[],[],[],[],...
        [0,0,-1,-1, -1],[1,1,1,1,1],[], options);
    
    fileID = fopen(filename , 'a+');
    fprintf(fileID, 'Optimization over life cycle:\n');
    fprintf(fileID, '--------------------\n');
    fprintf(fileID, 'x_opt = %.4f\n', x_opt_life);
    fprintf(fileID, 'fval = %.4f\n\n', fval_life);
    fclose(fileID);

%     % Optimization over life cycle (fixed allocation)
%     disp('Optimization over life cycle (fixed allocation)')
%     fun = @(x) PensionFundSE(x, E, 6, nSim, T, gamma);
%     [x_opt_life, fval_life, exitflag_life, output_life] = fmincon(...
%         fun,[-0.0230, -0.0236, 0.0029],[],[],[],[],...
%         [-1,-1, -1],[1,1,1],[], options);
%     
%     fileID = fopen(filename , 'a+');
%     fprintf(fileID, 'Optimization over life cycle (fixed allocation):\n');
%     fprintf(fileID, '--------------------\n');
%     fprintf(fileID, 'x_opt = %.4f\n', x_opt_life);
%     fprintf(fileID, 'fval = %.4f\n\n', fval_life);
%     fclose(fileID);

    % Optimization over the states
%     disp('Optimization over the states')
%     fun = @(x) PensionFund(x, E, 5, nSim, T, gamma);
%     [x_opt_states, fval_states, exitflag_states, output_states] = fmincon(...
%         fun, [0.77,0.99,-0.0126,-0.0137, 0.0021,0.0001,0.0003, 0.0003,-0.0026,0.0004, 0.0001], ...
%         [],[],[],[],[0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1],...
%         [1,1,  1,1,1,  1,1,1,  1, 1, 1],[], options);
%     
%     fileID = fopen(filename , 'a+');    
%     fprintf(fileID, 'Optimization over the states:\n');
%     fprintf(fileID, '--------------------\n');
%     fprintf(fileID, 'x_opt = %.4f\n', x_opt_states);
%     fprintf(fileID, 'fval = %.4f\n\n', fval_states);
%     fclose(fileID);
    
%     savefile = sprintf('results/gamma_%d-nSim_%d_SE.mat', gamma, nSim);
%     save(savefile, "-v7.3")
    end
end
