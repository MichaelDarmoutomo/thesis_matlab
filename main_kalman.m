%% Load data and dependencies
clear
load("thesis_data_short.MAT")
% load("thesis_data_param.MAT")

%% Set initial values

delta_pi = [0.02, -0.0028, -0.0014];
delta_r = [0.02, -0.0094, -0.0024];
K = [0.05, 0.5, 1.2];
sigma_pi = [-0.0010, 0.0013, 0.0055];
sigma_s = [-0.0483, 0.0078, 0.0010, 0.1335];
eta_s = 0.04;
lambda = [0.6, -0.02];
Lambda = [0.17, 0.4, -0.5, -1];
h = [0.0038, 0.0003, 0.0003, 0.0000, 0.0008, 0.0021];

init_param = [delta_pi, delta_r, K, sigma_pi, sigma_s, eta_s, lambda, Lambda, h];
clearvars delta_pi delta_r K sigma_pi sigma_s eta_s lambda Lambda h;

%% Optimizer
clearvars options;

options = optimset('fmincon');
options = optimset(options, 'Display', 'iter');
options = optimset(options , 'MaxFunEvals' , 1e+6);
options = optimset(options , 'MaxIter' ,1e+6);
options = optimset(options , 'TolFun' , 1e-8);
options = optimset(options , 'TolX' , 1e-6);

[x_opt,fval,exitflag,output] = fminunc(@(x) wrapper(x,data), init_param, options);
x_opt_struct = GetParameters(x_opt);

%% Compute standard errors
syms x
loss = @(x) wrapper(x,data);
hessian(loss, x)

%%
wrapper(x_opt, data)

%% Helper functions
function loss=wrapper(param, data)
% number of yields
m=6;

p_ = GetParameters(param);      % Split up array of parameters to individual
s = KalmanParameters(p_, m);    % Get the parameters for Kalman filter
[V, u] = KalmanFilter(s, data); % Apply Kalman filter on parameters
loss = LogLikelihood(V, u);     % Compute loss

end
