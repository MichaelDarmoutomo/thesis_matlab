function ec = GenerateEconomy(nSim,T, param)
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
dt           = 1;
ec.delta_0pi = param.delta_pi(1);
ec.delta_1pi = param.delta_pi(2:3)';
ec.delta_0r  = param.delta_r(1);
ec.delta_1r  = param.delta_r(2:3)';
sigma_Pi     = param.sigma_pi';
sigma_S      = param.sigma_s';
eta_S        = param.eta_s;
% Lambda is calculated such that sigma_S' * Lambda = 0
Lambda0     = [param.lambda'; 0; -sigma_S(1:2)'*param.lambda'/sigma_S(4)];
Lambda1     = [param.Lambda; 0, 0];
Lambda1     = [param.Lambda; 0, 0; ...
    -sigma_S(1:3)'*Lambda1(1:3,1)/sigma_S(4), ...
    -sigma_S(1:3)'*Lambda1(1:3,2)/sigma_S(4)];
ec.K        = param.K;
w_0         = 25000;
Sigma_X     = [eye(2); zeros(2)];


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
    ec.r(t,:)          = ec.delta_0r + ec.delta_1r' * ec.X(:,:,t);
    ec.pi(t,:)         = ec.delta_0pi + ec.delta_1pi' * ec.X(:,:,t);
    ec.Lambda(:,:,t)   = Lambda0 + Lambda1 * ec.X(:,:,t);
    ec.dPi_(t,:)       = ec.pi(t,:) * dt + sigma_Pi' * ec.dZ(:,:,t);
    ec.Pi(:,t+1)       = ec.Pi(:,t) + ec.Pi(:,t) .* ec.dPi_(t,:)';
    ec.dS_(:,t)        = (ec.r(t,:) + eta_S) * dt + sigma_S' * ec.dZ(:,:,t);
    ec.S(:,t+1)        = ec.S(:,t) + ec.S(:,t) .* ec.dS_(:,t);
    ec.dX(:,:,t)       = -ec.K * ec.X(:,:,t) .* dt + Sigma_X' * ec.dZ(:,:,t);
    ec.X(:,:,t+1)      = ec.X(:,:,t) +  ec.dX(:,:,t);
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

c = fA(0:nTau+1, param.K, param.Lambda, param.lambda, param.delta_r)';
d = fB(0:nTau+1, param.K, param.Lambda, param.delta_r);

% Bond prices
P = zeros(length(c), nSim, T);
for t = 1:T
    P(:,:,t) = exp(c + d*ec.X(:,:,t));
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