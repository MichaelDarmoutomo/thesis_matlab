function [ smoothedX,smoothedP,smoothedPcross,x0_out, P0_out ] = KalmanSmoother(param, y)

a =  param.a;
B =  param.B;
H =  param.H;
Q =  param.Q;
phi = param.phi;
Phi = param.Phi;

T = size(y, 2);
mu = inv(eye(2) - Phi(1:2, 1:2)) * phi(1:2);
DSigma = inv(eye(4) - kron(Phi(1:2, 1:2), Phi(1:2, 1:2))) * reshape(Q(1:2, 1:2), 1, [])';

X0 = [mu', 0, 0];
P0 = diag([DSigma([1,4])', 0 , 0]);

% Run the Kalman filter
[~,~,  X, Xhat, P, Phat] = KalmanFilter(param, y);

% Initialise the smoother at the end
smoothedX(:,T)   = X(:,T);
smoothedP(:,:,T)  = P(:,:,T);
        
% Run the Kalman smoother 
for j=1:(T-1)
    % Loop backwards through time
    t = T-j;
    
    % Run the smoother
    smoothedX(:,t)      = X(:,t) + (P(:,:,t) * Phi' * inv(Phat(:,:,t+1))) * (smoothedX(:,t+1) - Xhat(:,t+1));
    smoothedP(:,:,t)    = P(:,:,t) - (P(:,:,t) * Phi' * inv(Phat(:,:,t+1)) * Phi * P(:,:,t));
end

% Calculate the xi0 and P0 separately
tmp = Phat(:,:,1) \ (smoothedX(:,1) - Xhat(:,1));
tmp(isnan(tmp)) = 0;
x0_out = X0' + P0 * Phi' * tmp ;
P0_out  = P0 -  P0 * Phi' * ((Phat(:,:,1))\(Phat(:,:,1) - smoothedP(:,:,1))) * ((Phat(:,:,1))\(Phi * P0)) ;

% Run the loop to get the cross terms where smoothedPcross(t)=P_{t,t-1|T}
for t=1:T
    if t==1
        smoothedPcross(:,:,t) = smoothedP(:,:,t) * (Phat(:,:,t)\(Phi *P0));
    else
        smoothedPcross(:,:,t) = smoothedP(:,:,t) * ((Phat(:,:,t)\(Phi *P(:,:,t-1))));
    end
end

% Close the function
end

