function p = GetParameters(param, restricted)
%GETPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

if (~restricted)

    p.delta_pi = param(1:3);
    p.delta_r = param(4:6);
    K = param(7:9);
    p.K = [K(1), 0 ; K(2),K(3)];
    p.sigma_pi = [param(10:12), 0];
    p.sigma_s = param(13:16);
    p.eta_s = param(17);
    p.lambda = param(18:19);
    Lambda = param(20:23);
    p.Lambda = reshape(Lambda, 2, 2)';
    p.h = param(24:end).^2;
    
else
    p.delta_pi = [0, param(1:2)];
    p.delta_r = [0, param(3:4)];
    K = param(5:7);
    p.K = [K(1), 0 ; K(2),K(3)];
    p.sigma_pi = [param(8:10), 0];
    p.sigma_s = param(11:14);
    p.lambda = param(15:16);
    Lambda = param(17:20);
    p.Lambda = reshape(Lambda, 2, 2)';
    p.h = param(21:end).^2;

    % restrictions
    B_inf = inv(p.K + p.Lambda)' * p.delta_r(2:3)';
    p.delta_r(1) = log(1.021) +  p.lambda *  B_inf + 0.5 * (B_inf' * B_inf);
    p.eta_s = log(1.056) - p.delta_r(1) + 0.5 * p.sigma_s * p.sigma_s';
    p.delta_pi(1) = log(1.019) + 0.5 * p.sigma_pi * p.sigma_pi'; 
    
end

end


