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
% E   = GenerateEconomy(nSim,T, x_opt_struct);
load('data/economy_cov.mat')

%% 
% delete(gcp('nocreate'))
% parpool('threads')

% x = [
%     [0.4608,1.0000];
%     [0.4487,1.0000];
%     [0.4350,1.0000];
%     [0.4331,0.9998];
%     [0.4271,1.0000];
%     [0.4288,1.0000];
%     [0.4166,1.0000];
%     [0.4126,1.0000];
%     [0.4089,1.0000]];

x = [
    [0.5026,1.0000];
    [0.5151,1.0000];
    [0.5087,0.9995];
    [0.5000,1.0000];
    [0.4966,1.0000];
    [0.4916,0.9996];
    [0.4825,0.9997];
    [0.4756,0.9998];
    [0.4617,0.9999]];
res = struct;

filename = sprintf('results/CEC_gamma_%d-nSim_%d_smeren_cov.txt', gamma, nSim);
fileID = fopen(filename , 'a+');

for j=1:9
    disp(j)
    res.(['U' int2str(j)]) = [];
    res.(['CEC' int2str(j)]) = [];
    [res.(['U' int2str(j)]), res.(['CEC' int2str(j)])] = PensionFundSE(x(j,:), E, 3, nSim, T, gamma, j);
    fprintf(fileID, 'cec_%i = %.4f\n\n', j, res.(['CEC' int2str(j)]));
end

fclose(fileID);
save('results/res_x_cov.mat', "res")
%%
[U1, CEC1] = PensionFundSE([0.4765, 0.9995], E, 1, nSim, T, gamma, 10);

%%
for j = 1:9
    disp(sprintf('%.2f',(res.(['CEC' int2str(j)])/CEC1-1)*100))
end
%%
rho = 1 / (1 + mean(E.r,1:2));

Nbootstraps = 10000;
gains = struct;

for j=1:9
    gains.(['x' int2str(j)]) = [];
    gains.(['x' int2str(j)]) = bootstrap_welfare(U1, res.(['U' int2str(j)]), ...
        Nbootstraps, nSim, rho, Tp, gamma);
    disp(std(gains.(['x' int2str(j)])));
end
% %%
% 
% % U = U_1;
% % [~, idx] = find(sum(isinf(U_4),1)>0);
% % U(:,idx) = [];
% 
% CEC_1_dist = zeros(10000,1);
% rng(1194866)
% for i = 1:10000
% %     if mod(i, 100) == 0
% %         fprintf('Currently at %.2f percent. \n', i/10000*100);
% %     end
% %     sim_indices = randi(9991,9991,1);
% %     U_ = U(:,sim_indices);
% %     [size_U, ~] = size(U_);
% %     welfare= zeros(size_U, 1);
% %     for idx = 1:size_U
% %         tmp = U_(idx,:);
% %         welfare(idx,:) = mean(tmp(isfinite(tmp)), 'omitnan');
% %     end
%     welfare = U_1(:,i);
%     SW = sum(rho.^(100:length(welfare)) .* welfare(100:end)', 'omitnan');
%     CEC_1_dist(i) = -((SW*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
% end

% %%
% % [size_U, ~] = size(U_1);
% U = zeros(size_U, 1);
% for i = 1:size_U
%     tmp = U_1(i,:);
%     U(i,:) = mean(tmp(isfinite(tmp)));
% end
% % 
%%
function gains = bootstrap_welfare(U1, U2, Nbootstraps, nSim, rho, Tp, gamma)
    
    gains = zeros(Nbootstraps, 1);

    rng(1194866)

    parfor i = 1:Nbootstraps
%         if mod(i, 1000) == 0
%             fprintf('Currently at %.2f percent. \n', i/Nbootstraps*100);
%         end

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

