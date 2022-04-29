%%
clear
load('data/x_opt_restricted.mat')
load('data/se_restricted.MAT')
se_struct = GetParameters(se_restricted', restricted);

%% Settings of simulation
nSim = 10000;
T = 600;
Tp = 25;
gamma = 3;
rng(1194866)

%% Generate economy
E   = GenerateEconomy(nSim,T, x_opt_struct);

%%
delete(gcp('nocreate'))
parpool('threads')

%% 
disp('1')
x = [0.7815, 0.9992];
[U_1, cec_1] = PensionFund(x, E, 1, nSim, T, gamma, 10);
% 
% disp('2')
% x= [0.6039, 1.000];
% [U_2, cec_2] = PensionFund(x, E, 2, nSim, T, gamma);
% 
% disp('3')
% x= [0.7520, 1.0000];
% [U_3, cec_3]  = PensionFund(x, E, 3, nSim, T, gamma);

disp('4')
x= [0.7718, 0.9903, -0.0126, -0.0137, 0.0021];
[U_4, cec_4]  = PensionFund(x, E, 4, nSim, T, gamma, 10);

% x= [0.3955,0.8151,0.0002,0.000,0.0002,0,0.0002,-0.0002,-0.0001,0.0001,-0.0001];
% cec_5 = PensionFund(x, E, 5, nSim, T, gamma);

% filename = sprintf('results/CEC_gamma_%d-nSim_%d.txt', gamma, nSim);
% 
% fileID = fopen(filename , 'a+');
% fprintf(fileID, 'cec_1 = %.4f\n\n', cec_1);
% fprintf(fileID, 'cec_2 = %.4f\n\n', cec_2);
% fprintf(fileID, 'cec_3 = %.4f\n\n', cec_3);
% fprintf(fileID, 'cec_4 = %.4f\n\n', cec_4);
% fprintf(fileID, 'cec_5 = %.4f\n\n', cec_5);
% fclose(fileID);

%%
rho = 1 / (1 + mean(E.r,1:2));
% SW = zeros(nSim,1);

Nbootstraps = 10000;
% gains_2 = bootstrap_welfare(U_1, U_2, Nbootstraps, nSim, rho, Tp, gamma);

% gains_3 = bootstrap_welfare(U_1, U_3, Nbootstraps, nSim, rho, Tp, gamma);

gains_4 = bootstrap_welfare(U_1, U_4, Nbootstraps, nSim, rho, Tp, gamma);
%%
function gains = bootstrap_welfare(U1, U2, Nbootstraps, nSim, rho, Tp, gamma)
    
    gains = zeros(Nbootstraps, 1);

    rng(1194866)
    
    for i = 1:Nbootstraps
        if mod(i, 100) == 0
            fprintf('Currently at %.2f percent. \n', i/Nbootstraps*100);
        end

        sim_indices = randi(nSim,nSim,1);
        
        [size_U, ~] = size(U1);
        welfare= zeros(size_U, 1);
        for idx = 1:size_U
            tmp = U1(idx,sim_indices);
            welfare(idx,:) = mean(tmp(isfinite(tmp)));
        end

%         welfare = U1(:,sim_indices);
%         welfare = mean(welfare(isfinite(welfare), 2);
        SW = sum(rho.^(100:length(welfare)) .* welfare(100:end)');
        CEC1 = -((SW*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
        
        [size_U, ~] = size(U2);
        welfare= zeros(size_U, 1);
        for idx = 1:size_U
            tmp = U2(idx,sim_indices);
            welfare(idx,:) = mean(tmp(isfinite(tmp)));
        end

%         welfare = welfare(:,sim_indices);
%         welfare = mean(welfare, 2);
        SW = sum(rho.^(100:length(welfare)) .* welfare(100:end)');
        CEC2 = -((SW*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));

        gains(i) = CEC2/CEC1 - 1;
    end


end