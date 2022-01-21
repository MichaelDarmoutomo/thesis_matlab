function p = GetParameters(param)
%GETPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
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
end

