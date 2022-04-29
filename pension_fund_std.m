%%
clear
load('data/x_opt_restricted.mat')
load('data/se_restricted.MAT')
se_struct = GetParameters(se_restricted', restricted);

%% Settings of simulation
nSim = 10000;
T = 600;
Tp = 25;
gamma = 5;
rng(1194866)

%% Generate economy
E   = GenerateEconomy(nSim,T, x_opt_struct);
% load('data/economy_cov.mat')

%% 
delete(gcp('nocreate'))
parpool('threads')

disp('1')
x = [0.4877, 0.9995];
[U_1, cec_1] = PensionFund(x, E, 1, nSim, T, gamma, 10);

disp('2')
x= [0.4076, 1.0000];
[U_2, cec_2] = PensionFund(x, E, 2, nSim, T, gamma, 10);

disp('3')
x= [0.4052, 1.0000];
[U_3, cec_3]  = PensionFund(x, E, 3, nSim, T, gamma, 10);

disp('4')
x= [0.4865, 1.0000, 0.0001, 0.0001, 0.0008];
[U_4, cec_4]  = PensionFund(x, E, 4, nSim, T, gamma, 10);

% disp('5')
% x= [0.5217,0.8864,0.0100,0.0100,0.0006,-0.0027,-0.0007,0.0008,0.0044,0.0000,0.0005];
% [U_5, cec_5]  = PensionFundSE(x, E, 5, nSim, T, gamma);
%  
% filename = sprintf('results/CEC_gamma_%d-nSim_%d_COV.txt', gamma, nSim);
% 
% fileID = fopen(filename , 'a+');
% fprintf(fileID, 'cec_1 = %.4f\n\n', cec_1);
% fprintf(fileID, 'cec_2 = %.4f\n\n', cec_2);
% fprintf(fileID, 'cec_3 = %.4f\n\n', cec_3);
% fprintf(fileID, 'cec_4 = %.4f\n\n', cec_4);
% fprintf(fileID, 'cec_5 = %.4f\n\n', cec_5);
% fclose(fileID);
%%
idxs = find(mean(U_4, 1, 'omitnan') < -0.000001);
U1 = U_1;
U4 = U_4;
U1(:,idxs) = [];
U4(:,idxs) = [];
%% Regression results 
rho = 1 / (1 + mean(E.r,1:2));
for s=1:10000
    U1 = U_1(:,s);
    U2 = U_2(:,s);
    U3 = U_3(:,s);
    U4 = U_4(:,s);
    
    SW1(s) = sum(rho.^(100:length(U1)) .* U1(100:end)');
    SW2(s) = sum(rho.^(100:length(U2)) .* U2(100:end)');
    SW3(s) = sum(rho.^(100:length(U3)) .* U3(100:end)');
    SW4(s) = sum(rho.^(100:length(U4 )) .* U4(100:end)');
    
    CEC1(s) = -((SW1(s)*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
    CEC2(s) = -((SW2(s)*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
    CEC3(s) = -((SW3(s)*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
    CEC4(s) = -((SW4(s)*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
end

rstat1 = regstats(-CEC1',[squeeze(mean(E.P(11,:,:),3))', mean(E.r,1)',std(E.r,1)', mean(E.dS_,1)', std(E.dS_,[],1)']);
rstat2 = regstats(-CEC2',[squeeze(mean(E.P(11,:,:),3))', mean(E.r,1)',std(E.r,1)', mean(E.dS_,1)', std(E.dS_,[],1)']);
rstat3 = regstats(-CEC3',[squeeze(mean(E.P(11,:,:),3))', mean(E.r,1)',std(E.r,1)', mean(E.dS_,1)', std(E.dS_,[],1)']);
rstat4 = regstats(-CEC4',[squeeze(mean(E.P(11,:,:),3))', mean(E.r,1)',std(E.r,1)', mean(E.dS_,1)', std(E.dS_,[],1)']);

rstat1.rsquare
rstat2.rsquare
rstat3.rsquare
rstat4.rsquare

%%
rho = 1 / (1 + mean(E.r,1:2));
% SW = zeros(nSim,1);

Nbootstraps = 10000;
% gains_2 = bootstrap_welfare(U_1, U_2, Nbootstraps, nSim, rho, Tp, gamma);
% 
% gains_3 = bootstrap_welfare(U_1, U_3, Nbootstraps, nSim, rho, Tp, gamma);

gains_4 = bootstrap_welfare(U1, U4, Nbootstraps, length(U2), rho, Tp, gamma);

% gains_5 = bootstrap_welfare(U_1, U_5, Nbootstraps, nSim, rho, Tp, gamma);

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

