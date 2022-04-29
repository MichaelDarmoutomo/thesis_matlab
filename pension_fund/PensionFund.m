function [U_, CEC] = PensionFund(xx, E, afx, nSim, T, gamma, spread)
    % afx chooses adjustment factor, 1 = uniform, 2 = 3-2-1
    if afx == 4
        x = xx(1:2);
        beta = xx(3:end);
    elseif afx == 5
        x = xx(1:2);
        beta = [xx(3:5); xx(6:8); xx(9:11)]';
    else
        x = xx;
        beta= nan;
    end
    
    p       = 0.2;
    
    %% Initialization
    Tw = 40;
    Tp = 25;
    
    % Q-matrix
    Q = zeros(Tw + Tp, Tw + Tp);
    for i = 1:Tw+Tp
        for j=1:Tw+Tp
            if (i + j > Tw + 1) && (i + j < Tw + Tp + 2)
                Q(i,j) = 1;
            end
        end
    end

    rho     = 1 / (1 + mean(E.r,1:2));
    U_ = zeros(T-Tp+1, nSim);
    
    af_ = zeros(Tw+Tp,nSim);
    
    parfor (s = 1:nSim)
    
    % create B such that people receive 40% of their salary as pension
    B = zeros(Tw + Tp, Tw + Tp, T+1);
    disc = .4*[1/Tw:1/Tw:1 ones(1,Tp)];
    B(:,:,1) = E.w(1) * disc .* Q;
    
    A           = zeros(T+1, 1);
    L           = zeros(T+1, 1);
    discB       = zeros(Tw+Tp, Tw+Tp, T+1);
    Ink         = zeros(T, 1);
    Uit         = zeros(T, 1);
    Obl         = zeros(T, 1);
    CR          = ones(T, 1);
    v           = zeros(Tw+Tp, Tw+Tp, T+1);
    ind         = zeros(T, 1);
    ksi         = zeros(T, 1);
    U           = zeros(T-Tp+1, nSim);
    CRstrike    = 0;
    afag        = zeros(1, Tw+Tp);
    
    % Create alpha for gesloten smeren
    alpha=ones(Tw+Tp,1);
    for i=1:(spread-1)
        alpha(i)=i/spread;
    end
    
    af = zeros(5, Tw+Tp);

    % Uniform adjustment factor
    af(1,:) = ones(1,Tw+Tp);
    
    % 3-2-1 adjustment factor
    af(2,:) = [3*ones(1,Tw-10), 2*ones(1,10), ones(1,Tp)];
    
    % Optimized adjustment factor
    if afx == 4
        F =  @(beta) exp(beta(1)*(Tw+Tp-1:-1:0) + beta(2)*([Tw-1:-1:0 zeros(1,Tp)]) + beta(3).*([Tw-1:-1:0 zeros(1,Tp)].^2));
        af(4,:) = F(beta);
    end
    
    if afx == 5
        nu = [Tw+Tp-1:-1:0; Tw-1:-1:0 zeros(1,Tp); [Tw-1:-1:0 zeros(1,Tp)].^2];
    end
    
    % Needed for utility
    u       = @(x) x.^(1-gamma) / (1-gamma); % CRRA
    Xagg=sum(E.X,2)/nSim;
    Xagg2=zeros(2,T);
    Xagg2(:,:)=Xagg;
    rho     = 1/(1+E.delta_0r+mean(E.delta_1r'*Xagg2));
    
    B(:,:,2) = B(:,:,1);
    
    %% Run different ALMs
%     for s=1:nSim
%         if mod(s, 100) == 0
%             fprintf('Currently at %.2f percent. \n', s/nSim*100);
%         end
        discB(:,:,1)  = E.P(:,s,1) .* B(:,:,1);
        discB(:,:,2)  = E.P(:,s,2) .* B(:,:,2);
        
        if afx == 3
            % Equal changes to achievable pension
            af(3,:) = ones(1,Tw+Tp);
            W       = zeros(Tw+Tp,T);
            W(:,1)  = sum(discB(:,:,1),1)';
            H       = zeros(Tw+Tp,T);
            sum_    = zeros(2,2,Tw);
            
            sum_(:,:,1) = expm(-E.K);
            for i = 2:Tw
                sum_(:,:,i) = sum_(:,:,i-1) + expm(-i * E.K);
            end
            w_tilde = zeros(T,Tw);
            W(:,1)  = sum(discB(:,:,1),1)';
            for tau=1:Tw
                w_tilde(1,tau)  = E.w(1,s) * exp(tau * E.delta_0pi + E.delta_1pi' * sum_(:,:,tau) * E.X(:,s,1));
            end
            for j = 1:Tw
                H(j,1) =  p * E.P(1:Tw-j+1,s,1)' * w_tilde(1,1:Tw-j+1)';
            end
            af(3,:) = (W(:,1) + H(:,1))' ./ W(:,1)';
            afag = afag + af(3,:);
        end
        L(1)  = sum(discB(:,:,1), "all");
        A(1)  = L(1);
        
        for t = 2:T
            if afx == 5
                af(5,:) = exp(nu'*beta*[1;E.X(:,s,t)])';
            end
            B_year = B(:,:,t-1);
            P = E.P(:,s,t);
            P2 = E.P2(:,s,t);
            discB(:,:,t) = P .* B_year;
            L(t)         = sum(discB(:,:,t), 'all');
            A(t)=A(t-1)*(x(1)*(1+E.dS_(t,s))+(1-x(1)-x(2))*(1+E.r(t,s)))+x(2)*CR(t-1)*L(t);
            CR(t)    = A(t) / L(t);
            
            if CR(t) < 1
                CRstrike = CRstrike + 1;
            else
                CRstrike = 0;
            end
            if CRstrike == 5
                ind(t) = CR(t)-1;
                ksi(t)      = ind(t) * L(t) / sum(P .* ((alpha'.*af(afx,:)) .* B_year), 'all');
                v(:,:,t)    = ksi(t) * alpha * af(afx,:);
            elseif CR(t) > 0.9 && CR(t) < 1.2
                ind(t) = (A(t) - L(t)) / ((spread-1)*A(t) + L(t));
                ksi(t)      = ind(t) * L(t) / sum(P .* (af(afx,:) .* B_year), 'all');
                v(:,:,t)    = ksi(t) * ones(Tw+Tp,1) * af(afx,:);
            elseif CR(t) > 1.20
                ind(t) = (A(t) - L(t)) / ((spread/2-1)*A(t) + L(t));
                ksi(t)      = ind(t) * L(t) / sum(P .* (af(afx,:) .* B_year), 'all');
                v(:,:,t)    = ksi(t) * ones(Tw+Tp,1) * af(afx,:);
            else %backstop 1
                ind(t) = CR(t) / 0.9 - 1;
                ksi(t)      = ind(t) * L(t) / sum(P .* ((alpha'.*af(afx,:)) .* B_year), 'all');
                v(:,:,t)    = ksi(t) * alpha * af(afx,:);
            end
            Ink(t)      = Tw * p * E.w(t,s);
            B_year    = (1 + v(:,:,t)) .* B_year;
            
            % Nonnegativity restriction for pensions
            B_year    = (abs(B_year) + B_year) / 2;
            
            discB(:,:,t) = P .* B_year;
            L(t)         = sum(discB(:,:,t), 'all');
            CR(t)        = A(t) / L(t);
            Uit(t)       = sum(B_year(1,Tw+1:Tw+Tp));

            A(t)=A(t)-Uit(t)+Ink(t);
            B_tmp = B_year(2:Tw+Tp,1:Tw+Tp-1);
            B_year = zeros(Tw+Tp, Tw+Tp);
            B_year(1:Tw+Tp-1,2:Tw+Tp) = B_tmp;

            payout = (p * E.w(t,s)) ./ sum(P2 .* Q, 1);
            payout(Tw+1:end) = 0;
            B(:,:,t) = (B_year + payout) .* Q;
            
            discB(:,:,t)    = P2 .* B(:,:,t);
            L(t)            = sum(discB(:,:,t), 'all');
            CR(t)=A(t)/L(t);
            
            if afx == 3
                W(:,t)  = sum(discB(:,:,t),1)';
                for tau=1:Tw
                    w_tilde(t,tau) = E.w(t,s) * exp(tau * E.delta_0pi + E.delta_1pi' * sum_(:,:,tau) * E.X(:,s,t));
                end
                for jt=1:Tw
                    H(jt,t)  = p * E.P(1:Tw-jt+1,s,t)' * w_tilde(t,1:Tw-jt+1)';
                end
                af(3,:)  = (W(:,t) + H(:,t))' ./ W(:,t)';
                af_(:,s) = af(3,:);
            end
        end
        
        %% Compute utility
        rho_ = 1 / (mean(E.r(:,s)) + 1);
        
        % Loops over generations j
%         for j = -Tw:T-Tw-Tp
%             U2(j+Tw+1,s) = 0;
%             % Loops over pension period (25 years) of generation j
%             for t = j+Tw+1:j+Tw+Tp
%                 U2(j+Tw+1,s) = U2(j+Tw+1,s) + rho_^(t-j-Tw-1) * u( B(1,t-j,t)/E.Pi(t,s));
%             end
%         end
        
        % Vectorized version
        rho_sq = rho_.^(0:24);
        B_cut = squeeze(B(1, 41:end, :));
        for j = -Tw:T-Tw-Tp
            U(j+Tw+1, s) = sum(rho_sq .* u(diag(B_cut,j+Tw) ./ E.Pi(j+Tw+1:j+Tw+Tp,s))');
        end
        
        U_(:,s) = U(:,s);
    end

    if (sum(isinf(U_), 1:2) > 0)
        [size_U, ~] = size(U_);
        welfare= zeros(size_U, 1);
        for idx = 1:size_U
            tmp = U_(idx,:);
            welfare(idx,:) = mean(tmp(isfinite(tmp)), 'omitnan');
        end
    else
        welfare = mean(U_, 2, 'omitnan');
    end

%     welfare = mean(U_, 2);
    SW = sum(rho.^(100:length(welfare)) .* welfare(100:end)');
    CEC = -((SW*(1-rho)^2 * (1-gamma)) / ((1-rho^Tp)*rho^101))^((1-gamma)^(-1));
end