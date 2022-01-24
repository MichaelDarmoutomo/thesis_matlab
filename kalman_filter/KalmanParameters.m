function s=KalmanParameters(param, m)
%KALMANPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
delta_pi = param.delta_pi;
delta_r = param.delta_r;
K = param.K;
sigma_pi = param.sigma_pi;
sigma_s = param.sigma_s;
eta_s = param.eta_s;
lambda = param.lambda;
Lambda = param.Lambda;
h = param.h;

dt = 1/12;

a_ = [
    0;
    0; 
    delta_pi(1) - 0.5 * (sigma_pi * sigma_pi'); 
    delta_r(1) + eta_s - 0.5 * (sigma_s * sigma_s')];

A_ = [[-K;delta_pi(2:3);delta_r(2:3)], zeros(4,2)];

C_ = [[eye(2), zeros(2,2)]; sigma_pi; sigma_s];

[U, D] = eig(A_);
D = diag(D);
Uinv = inv(U);
F_ = diag(dt * alpha(D * dt));

phi = (U * F_ * Uinv) * a_;
Phi = expm(A_ * dt);

V = zeros(4, 4);
for (i=1:4)
    for (j=1:4)
        V_tmp = (Uinv * (C_ * C_') * Uinv');
        V(i,j) = V_tmp(i,j) * dt * alpha((D(i) + D(j)) * dt);
    end
end

Q = U * V * U';

a_tmp = fA([1,5,10,15,20,30], K, Lambda, lambda, delta_r) ./ [1,5,10,15,20,30];
% a_tmp = fA([1,5,10,15,20,30], A_, lambda, delta_r) ./ [1,5,10,15,20,30];
a = [-a_tmp, 0, 0];

B = zeros(m+2, 4);
B_tmp = fB([1,5,10,15,20,30], K, Lambda, delta_r) ./ [1,5,10,15,20,30]';
% B_tmp = fB([1,5,10,15,20,30], A _) ./ [1,5,10,15,20,30]';
B(1:m, 1:2) = -B_tmp;
B((m+1):(m+2), 3:4) = eye(2);

H = zeros(m+2);
H(1:m, 1:m) = diag(h);

s.a = a;
s.B = B;
s.H = H;
s.Q = Q;
s.phi = phi;
s.Phi = Phi;
end

function res=alpha(x)
res = (exp(x) - 1) ./ x;
res(isnan(res))=1;
end

function res=fB(tau, K, Lambda, delta_r)
M = (K + Lambda)';
Minv = inv(M);
res = zeros(length(tau), 2);
for i=1:length(tau)
    res(i,:) = Minv * (expm(-tau(i) * M) - eye(2)) * delta_r(2:3)';
end
end

function res=fAprime(tau, K, Lambda, lambda, delta_r)
res = zeros(1,length(tau));
for i=1:length(tau)
    B = fB(tau(i),  K, Lambda, delta_r);
    res(i) = -B * lambda' + 0.5 * (B * B') - delta_r(1);
end
end

function res=fA(tau, K, Lambda, lambda, delta_r)
tau = [0,tau];
p = zeros(1,length(tau)-1);
for i = 1:length(tau)-1
    p(i) = integral(@(x) fAprime(x, K, Lambda, lambda, delta_r), tau(i), tau(i+1));
end
res=cumsum(p);
end

% 
% function res=fB(tau, A_)
% res = zeros(length(tau), 2);
% for i=1:length(tau)
%     tmp = expm(A_* tau(i));
%     res(i,:) = tmp(4, 1:2);
% end
% end
% 
% function res=fAprime(tau, A_, lambda, delta_r)
% res = zeros(1,length(tau));
% for i=1:length(tau)
%     B = fB(tau(i), A_);
%     res(i) = B * lambda' + 0.5 * (B * B') - delta_r(1);
% end
% end
% 
% function res=fA(tau, A_, lambda, delta_r)
% tau = [0,tau];
% p = zeros(1,length(tau)-1);
% for i = 1:length(tau)-1
%     p(i) = integral(@(x) fAprime(x, A_, lambda, delta_r), tau(i), tau(i+1));
% end
% res=cumsum(p);
% end
