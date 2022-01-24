function ec = GenerateEconomySE(nSim,T, param, se)
% ec = generateEconomy(nSim,T) returns economic variables for nSim simulations and
% a period of T years. Function returns a 1 x 1 structure with: 
% dZ     := 4 x nSim x T        
% Lambda := 4 x nSim x T
% r      := T x nSim
% pi     := T x nSim
% dPi_   := T x nSim
% Pi     := T x nSim
% dS_    := T x nSim
% S      := T x nSim
% dX     := 2 x nSim x T
% X      := 2 x nSim x T
% P      := nTau x nSim x T
% w      := T x nSim

%% Initialization
nSim = nSim * 10;    


dt           = 1;
ec.delta_0pi = normrnd(param.delta_pi(1), se.delta_pi(1), [1, nSim]);
ec.delta_1pi = [
    normrnd(param.delta_pi(2), se.delta_pi(2), [1, nSim]);
    normrnd(param.delta_pi(3), se.delta_pi(3), [1, nSim])];
ec.delta_0r  = normrnd(param.delta_r(1), se.delta_r(1), [1, nSim]);
ec.delta_1r  = [
    normrnd(param.delta_r(2), se.delta_r(2), [1, nSim]);
    normrnd(param.delta_r(3), se.delta_r(3), [1, nSim])];
sigma_Pi     = [
    normrnd(param.sigma_pi(1), se.sigma_pi(1), [1, nSim]);
    normrnd(param.sigma_pi(2), se.sigma_pi(2), [1, nSim]);
    normrnd(param.sigma_pi(3), se.sigma_pi(3), [1, nSim]);
    zeros(1,nSim)];
sigma_S      = [
    normrnd(param.sigma_s(1), se.sigma_s(1), [1, nSim]);
    normrnd(param.sigma_s(2), se.sigma_s(2), [1, nSim]);
    normrnd(param.sigma_s(3), se.sigma_s(3), [1, nSim]);
    normrnd(param.sigma_s(4), se.sigma_s(4), [1, nSim])];
eta_S        = normrnd(param.eta_s, se.eta_s, [1, nSim]);
% Lambda is calculated such that sigma_S' * Lambda = 0

lambda       = [
    normrnd(param.lambda(1), se.lambda(1), [1, nSim]);
    normrnd(param.lambda(2), se.lambda(2), [1, nSim])];

Lambda0      = [lambda; zeros(1,nSim); sum(-sigma_S(1:2,:) .* lambda) ./ sigma_S(4,:)];

Lambda       = [
    normrnd(param.Lambda(1,1), se.Lambda(1,1), [1, nSim]);
    normrnd(param.Lambda(1,2), se.Lambda(1,2), [1, nSim]);
    normrnd(param.Lambda(2,1), se.Lambda(2,1), [1, nSim]);
    normrnd(param.Lambda(2,2), se.Lambda(2,2), [1, nSim])];

Lambda1 = zeros(4,2,nSim);
for iSim=1:nSim
    Lambda1(:,:,iSim) = [
        Lambda(1,iSim), Lambda(2,iSim);
        Lambda(3,iSim), Lambda(4,iSim);
        0, 0;
        [sum(-sigma_S(1:3, iSim)' .* [Lambda(1,iSim), Lambda(3,iSim), 0]) / sigma_S(4, iSim),
        sum(-sigma_S(1:3, iSim)' .* [Lambda(2,iSim), Lambda(4,iSim), 0]) / sigma_S(4, iSim)]'];
end

K         = [
    normrnd(param.K(1,1), se.K(1,1), [1, nSim]);
    normrnd(param.K(2,1), se.K(2,1), [1, nSim]);
    normrnd(param.K(2,2), se.K(2,2), [1, nSim])];

ec.K = zeros(2,2,nSim);
for iSim=1:nSim
    ec.K(:,:,iSim) = [
        K(1,iSim), 0;
        K(2,iSim), K(3,iSim)];
end

w_0          = 25000;
Sigma_X      = [eye(2); zeros(2)];

counter = 0;
idx = zeros(1, nSim);
for s=1:nSim
    eigenvalues = eig(ec.K(:,:,s) + Lambda1(1:2,:,s));
    if ((sum(eigenvalues > 0.02) ~= 2) || (~isreal(eigenvalues)))
        counter = counter + 1;
        idx(:,s) = 1;
    end
end

nSim = nSim / 10;
% TODO: select only ~idx, and truncate to nSim.

ec.delta_0pi = ec.delta_0pi(:,~idx);
ec.delta_1pi = ec.delta_1pi(:,~idx);
ec.delta_0r = ec.delta_0r(:,~idx);
ec.delta_1r = ec.delta_1r(:,~idx);
sigma_Pi = sigma_Pi(:,~idx);
sigma_S = sigma_S(:,~idx);
eta_S = eta_S(:,~idx);
Lambda0 = Lambda0(:,~idx);
Lambda1 = Lambda1(:,:,~idx);
ec.K = ec.K(:,:,~idx);

ec.delta_0pi = ec.delta_0pi(:,1:nSim);
ec.delta_1pi = ec.delta_1pi(:,1:nSim);
ec.delta_0r = ec.delta_0r(:,1:nSim);
ec.delta_1r = ec.delta_1r(:,1:nSim);
sigma_Pi = sigma_Pi(:,1:nSim);
sigma_S = sigma_S(:,1:nSim);
eta_S = eta_S(:,1:nSim);
Lambda0 = Lambda0(:,1:nSim);
Lambda1 = Lambda1(:,:,1:nSim);
ec.K = ec.K(:,:,1:nSim);



%% Simulation of economy
ec.dZ          = normrnd(0, dt, [4, nSim, T]);
ec.Lambda      = zeros(4, nSim, T);
ec.r           = zeros(T,nSim);
ec.pi          = zeros(T,nSim);
ec.dPi_        = zeros(T,nSim);
ec.Pi          = zeros(nSim,T+1);
ec.dS_         = zeros(nSim,T);
ec.S           = zeros(nSim,T+1);
ec.dX          = zeros(2, nSim, T);
ec.X           = zeros(2, nSim, T+1);

ec.X(:,:,1)    = ones(2,nSim).* [2.330; 0.546];
ec.Pi(:,1)     = ones(nSim,1);
ec.S(:,1)      = ones(nSim,1);

for t = 1:T
    ec.r(t,:)          = ec.delta_0r + sum(ec.delta_1r .* ec.X(:,:,t));
    ec.pi(t,:)         = ec.delta_0pi + sum(ec.delta_1pi .* ec.X(:,:,t));
end

for t = 1:T
    for s=1:nSim
        ec.Lambda(:,s,t)   = Lambda0(:,s) + Lambda1(:,:,s) * ec.X(:,s,t);
        ec.dPi_(t,s)       = ec.pi(t,s) * dt + sigma_Pi(:,s)' * ec.dZ(:,s,t);
        ec.Pi(s,t+1)       = ec.Pi(s,t) + ec.Pi(s,t) .* ec.dPi_(t,s)';
        ec.dS_(s,t)        = (ec.r(t,s) + eta_S(:,s)) * dt + sigma_S(:,s)' * ec.dZ(:,s,t);
        ec.S(s,t+1)        = ec.S(s,t) + ec.S(s,t) .* ec.dS_(s,t);
        ec.dX(:,s,t)       = -ec.K(:,:,s) * ec.X(:,s,t) .* dt + Sigma_X' * ec.dZ(:,s,t);
        ec.X(:,s,t+1)      = ec.X(:,s,t) +  ec.dX(:,s,t);
    end
end

ec.Pi   = ec.Pi(:,1:end-1)';
ec.dS_  = ec.dS_';
ec.S    = ec.S(:,1:end-1)';
ec.X    = ec.X(:,:,1:end-1);

%% Simulate bond prices
nTau    = 64;   % Number of maturities
% c       = zeros(nTau,1);
% d       = zeros(nTau,2);
% 
% D = @(tau) (ec.K'+Lambda1'*Sigma_X) \ (expm(-tau*(ec.K'+Lambda1'*Sigma_X)) - ...
%     eye(2)) * ec.delta_1r;
% Cdot = @(s) -ec.delta_0r - Lambda0' * Sigma_X * D(s) + 0.5 * D(s)' * (Sigma_X' * Sigma_X) * D(s);
% 
% for i = 0:nTau
%     c(i+1)    = C(i, Cdot, 500);
%     d(i+1,:)  = D(i); 
% end

% c = fA(0:nTau+1, param.K, param.Lambda, param.lambda, param.delta_r)';
% d = fB(0:nTau+1, param.K, param.Lambda, param.delta_r);

c = zeros(nTau+2,nSim);
d = zeros(nTau+2, 2, nSim);
for s=1:nSim
    if mod(s, 10) == 0
        disp(s)
    end

    c(:,s) = fA(0:nTau+1, ec.K(:,:,s), Lambda1(1:2,:,s), Lambda0(1:2,s)', [ec.delta_0r(:,s);ec.delta_1r(:,s)]')';
    d(:,:,s) = fB(0:nTau+1, ec.K(:,:,s), Lambda1(1:2,:,s), [ec.delta_0r(:,s);ec.delta_1r(:,s)]');
end


% Bond prices
P = zeros(size(c,1), nSim, T);
for t = 1:T
    for s=1:nSim
        P(:,s,t) = exp(c(:,s) + d(:,:,s)*ec.X(:,s,t));
    end
end

ec.P = P(1:end-1,:,:);
ec.P2 = P(2:end,:,:);

% c2=c;
% d2=d;
% c2(1:end-1)=c(2:end);
% d2(1:end-1,:)=d(2:end,:);
% c2(end)=C(nTau+1,Cdot,50);
% d2(end,:)=D(nTau+1);
% for t = 1:T
%     ec.P2(:,:,t) = exp(c2 + d2*ec.X(:,:,t));
% end

% Wages
ec.w = w_0 * ec.Pi;


%% FUNCTIONS

function c = C(tau, Cdot, intervals)
% use riemann sum approximation
c = 0;
for rep = 1:intervals
    c = c + Cdot(rep*tau/(intervals));
end
c = c * tau / intervals;
end

function res=fB(tau, K, Lambda, delta_r)
M = (K + Lambda)';
res = zeros(length(tau), 2);
for iB=1:length(tau)
    res(iB,:) = M \ (expm(-tau(iB) * M) - eye(2)) * delta_r(2:3)';
end
end

function res=fAprime(tau, K, Lambda, lambda, delta_r)
res = zeros(1,length(tau));
for iAprime=1:length(tau)
    B = fB(tau(iAprime),  K, Lambda, delta_r);
    res(iAprime) = -B * lambda' + 0.5 * (B * B') - delta_r(1);
end
end

function res=fA(tau, K, Lambda, lambda, delta_r)
    tau = [0,tau];
    p = zeros(1,length(tau)-1);
    for iA = 1:length(tau)-1
        p(iA) = integral(@(x) fAprime(x, K, Lambda, lambda, delta_r), tau(iA), tau(iA+1));
    end
    res=cumsum(p);
end

end