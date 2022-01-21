function [V, u] = KalmanFilter(param, y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

% Initialize variables
Xhat = zeros(4, T);
X = zeros(4, T);
Phat = zeros([4,4,T]);
P = zeros([4,4,T]);
V = zeros([size(H), T]);
u = zeros(size(y));
K = zeros([size(B'), T]);

for t=1:T    
    % Prediction for t == 1
    if (t == 1)
      X(:,t) = X0;
      Xhat(:,t) = X0;
      P(:,:,t) = P0;
      Phat(:,:,t) = P0;
    else 
      % Prediction step for any other t
      Xhat(:,t) = phi + Phi * X(:,t-1); % Predicted state
      Phat(:,:,t) = Phi * P(:,:,(t-1)) * Phi' + Q; % Predicted covariance (Q is variance of state eq.)
    
      % Variables needed
      V(:,:,t) = B * Phat(:,:,t) * B' + H;
  
      u(:,t) = y(:,t) - (a' + B * Xhat(:,t));
      K(:,:,t) = Phat(:,:,t) * B'* inv(V(:,:,t));
      
      % State update
      X(:,t) = Xhat(:,t) + K(:,:,t) * u(:,t);
      
      % P update
      P(:,:,t) = (eye(4) - (K(:,:,t) * B)) * Phat(:,:,t);
    end
end
V = V(:,:,2:end);
u = u(:,2:end);
end

