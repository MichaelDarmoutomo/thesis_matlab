function ec = GenerateEconomySE(nSim,T, param, se, x_opt, covariance)
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


nSim0 = nSim;
nSim = nSim * 100;    

params = mvnrnd(x_opt, covariance, nSim);

dt           = 1;

% Covariance matrix method
ec.delta_1pi = params(:,1:2)';
ec.delta_1r = params(:, 3:4)';
sigma_Pi = [params(:,8:10)'; zeros(1,nSim)];
sigma_S = params(:,11:14)';

lambda = params(:,15:16)';
Lambda0 = [lambda; zeros(1,nSim); sum(-sigma_S(1:2,:) .* lambda) ./ sigma_S(4,:)];
Lambda = params(:,17:20)';

% ec.delta_1pi = [
%     normrnd(param.delta_pi(2), se.delta_pi(2), [1, nSim]);
%     normrnd(param.delta_pi(3), se.delta_pi(3), [1, nSim])];
% ec.delta_1r  = [
%     normrnd(param.delta_r(2), se.delta_r(2), [1, nSim]);
%     normrnd(param.delta_r(3), se.delta_r(3), [1, nSim])];
% sigma_Pi     = [
%     normrnd(param.sigma_pi(1), se.sigma_pi(1), [1, nSim]);
%     normrnd(param.sigma_pi(2), se.sigma_pi(2), [1, nSim]);
%     normrnd(param.sigma_pi(3), se.sigma_pi(3), [1, nSim]);
%     zeros(1,nSim)];
% sigma_S      = [
%     normrnd(param.sigma_s(1), se.sigma_s(1), [1, nSim]);
%     normrnd(param.sigma_s(2), se.sigma_s(2), [1, nSim]);
%     normrnd(param.sigma_s(3), se.sigma_s(3), [1, nSim]);
%     normrnd(param.sigma_s(4), se.sigma_s(4), [1, nSim])];
% Lambda is calculated such that sigma_S' * Lambda = 0

% lambda       = [
%     normrnd(param.lambda(1), se.lambda(1), [1, nSim]);
%     normrnd(param.lambda(2), se.lambda(2), [1, nSim])];
% 
% Lambda0      = [lambda; zeros(1,nSim); sum(-sigma_S(1:2,:) .* lambda) ./ sigma_S(4,:)];
% 
% Lambda       = [
%     normrnd(param.Lambda(1,1), se.Lambda(1,1), [1, nSim]);
%     normrnd(param.Lambda(1,2), se.Lambda(1,2), [1, nSim]);
%     normrnd(param.Lambda(2,1), se.Lambda(2,1), [1, nSim]);
%     normrnd(param.Lambda(2,2), se.Lambda(2,2), [1, nSim])];

Lambda1 = zeros(4,2,nSim);
% for iSim=1:nSim
% 
% end

% K         = [
%     normrnd(param.K(1,1), se.K(1,1), [1, nSim]);
%     normrnd(param.K(2,1), se.K(2,1), [1, nSim]);
%     normrnd(param.K(2,2), se.K(2,2), [1, nSim])];

K = params(:,5:7)';

ec.K = zeros(2,2,nSim);

for iSim=1:nSim
    Lambda1(:,:,iSim) = [
        Lambda(1,iSim), Lambda(2,iSim);
        Lambda(3,iSim), Lambda(4,iSim);
        0, 0;
        [sum(-sigma_S(1:3, iSim)' .* [Lambda(1,iSim), Lambda(3,iSim), 0]) / sigma_S(4, iSim),
        sum(-sigma_S(1:3, iSim)' .* [Lambda(2,iSim), Lambda(4,iSim), 0]) / sigma_S(4, iSim)]'];

    ec.K(:,:,iSim) = [
        K(1,iSim), 0;
        K(2,iSim), K(3,iSim)];

    B_inf                = inv(ec.K(:,:,iSim) + Lambda1(1:2,:,iSim))' * ec.delta_1r(:,iSim);
    ec.delta_0r(:,iSim)  = log(1.021) +  Lambda0(1:2,iSim)' * B_inf + 0.5 * (B_inf' * B_inf);
    eta_S(:,iSim)        = log(1.056) - ec.delta_0r(:,iSim) + 0.5 * (sigma_S(:,iSim)' * sigma_S(:,iSim));
    ec.delta_0pi(:,iSim) = log(1.019) + 0.5 * (sigma_Pi(:,iSim)' * sigma_Pi(:,iSim)); 
end

w_0          = 25000;
Sigma_X      = [eye(2); zeros(2)];

counter = 0;
idx = zeros(1, nSim);
counter2=1;
for s=1:nSim
    eigenvalues = eig(ec.K(:,:,s) + Lambda1(1:2,:,s));
    eigenvaluesK = eig(ec.K(:,:,s));
    eigenvaluesLambda = eig(Lambda1(1:2,:,s));
    if ((sum(eigenvalues > 0.0275) ~= 2 || sum(eigenvaluesK > 0.01) ~= 2) ...
            || (~isreal(eigenvalues) || ~isreal(eigenvaluesK)))
        counter = counter + 1;
        idx(:,s) = 1;        
    else
        eigenvalues_list(:,counter2) = eigenvalues;
        counter2 = counter2 + 1;
%         disp(eigenvaluesLambda)
    end

end

disp(counter)
nSim = nSim / 100;
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
    for s=1:nSim
        ec.r(t,s)          = ec.delta_0r(:,s) + sum(ec.delta_1r(:,s)' * ec.X(:,s,t));
        ec.pi(t,s)         = ec.delta_0pi(:,s) + sum(ec.delta_1pi(:,s)' * ec.X(:,s,t));
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
D = @(tau,s) (ec.K(:,:,s)'+Lambda1(:,:,s)'*Sigma_X) \ (expm(-tau*(ec.K(:,:,s)'+Lambda1(:,:,s)'*Sigma_X)) - ...
    eye(2)) * ec.delta_1r(:,s);
Cdot = @(tau, s) -ec.delta_0r(:,s) - Lambda0(:,s)' * Sigma_X * D(tau,s) + 0.5 * D(tau,s)' * (Sigma_X' * Sigma_X) * D(tau,s);
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
% indices = setdiff(1:nSim, find(ec.P(65,:,600) > 1));
% 
% ec.delta_0pi = ec.delta_0pi(:,indices);
% ec.delta_1pi = ec.delta_1pi(:,indices);
% ec.delta_0r = ec.delta_0r(:,indices);
% ec.delta_1r = ec.delta_1r(:,indices);
% ec.K = ec.K(:,:,indices);
% ec.dZ = ec.dZ(:,indices,:);
% ec.Lambda = ec.Lambda(:,indices,:);
% ec.r = ec.r(:,indices);
% ec.pi = ec.pi(:,indices);
% ec.dPi_ = ec.dPi_(:,indices);
% ec.Pi = ec.Pi(:,indices);
% ec.dS_ = ec.dS_(:,indices);
% ec.S = ec.S(:,indices);
% ec.dX = ec.dX(:,indices,:);
% ec.X = ec.X(:,indices,:);
% ec.P = ec.P(:,indices,:);
% ec.P2 = ec.P2(:,indices,:);
% 
% ec.delta_0pi = ec.delta_0pi(:,1:nSim0);
% ec.delta_1pi = ec.delta_1pi(:,1:nSim0);
% ec.delta_0r = ec.delta_0r(:,1:nSim0);
% ec.delta_1r = ec.delta_1r(:,1:nSim0);
% ec.K = ec.K(:,:,1:nSim0);
% ec.dZ = ec.dZ(:,1:nSim0,:);
% ec.Lambda = ec.Lambda(:,1:nSim0,:);
% ec.r = ec.r(:,1:nSim0);
% ec.pi = ec.pi(:,1:nSim0);
% ec.dPi_ = ec.dPi_(:,1:nSim0);
% ec.Pi = ec.Pi(:,1:nSim0);
% ec.dS_ = ec.dS_(:,1:nSim0);
% ec.S = ec.S(:,1:nSim0);
% ec.dX = ec.dX(:,1:nSim0,:);
% ec.X = ec.X(:,1:nSim0,:);
% ec.P = ec.P(:,1:nSim0,:);
% ec.P2 = ec.P2(:,1:nSim0,:);

ec.w = w_0 * ec.Pi;

%% FUNCTIONS

function c = C(tau, Cdot, intervals ,s)
% use riemann sum approximation
c = 0;
for rep = 1:intervals
    c = c + Cdot(rep*tau/(intervals),s);
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