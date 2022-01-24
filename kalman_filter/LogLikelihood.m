function loss = LogLikelihood(V,u)
%LOGLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here

loss = 0;
T = size(u, 2);
for t=1:T
    loss_ = -0.5 * log(det(V(:,:,t))) - 0.5 * u(:,t)' * inv(V(:,:,t)) * u(:,t);
    loss = loss + loss_;
end

loss = -loss;
