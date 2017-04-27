function [grad,labda,J_partial,gradNp,grad_ss] = get_adjoint_Nc(options,sol,Wp,index)

Np          = options.Np;
Nc          = options.Nc;
k           = length(sol);
labda       = zeros(k,Np);
% labdaCsum   = zeros(1,min(size(sol)));
% J_betasum   = zeros(1,min(size(sol)));
grad        = zeros(Nc,Wp.turbine.N);

tic
for n = Np:-1:1
    % Note: m = n - 1 // o = n + 1
    
    load(strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(n),'.mat'));
    Cn_xn           = sys.A;
    Cn_xm           = sys.derivatives.dAdx - sys.derivatives.dSmdx ...
                        - sys.derivatives.Q - sys.derivatives.dBc;
    Cn_betan        = - sys.derivatives.dSmdbeta;
    Jn_xm           = sys.derivatives.dJdx';
    J_betan         = sys.derivatives.dJdbeta;
    J_partial(n,:)  = J_betan;
    
    % Steady state
    Cx_ss           = sys.derivatives.dAdx + sys.A ...
                     - sys.derivatives.dSmdx - sys.derivatives.Q - sys.derivatives.dBc;
	labda_ss(:,n)   = - Cx_ss'\Jn_xm';
    grad_ss(n,:)    = J_betan + labda_ss(:,n)'*Cn_betan;
      
    if n == Np && Np == Nc  % ONLY if Np = Nc
        grad(n,:)           = J_betan;
        gradNp(n,:)         = grad(n,:);
    elseif n == Np          % Final value of grad en first of labda
        labda(:,n)          = zeros(k,1);
        grad(Nc,:)          = J_betan;
        gradNp(n,:)         = J_betan;
    elseif n >= Nc           % Everything from Nc to Np
        labda(:,n)          = - Cn_xn'\(Jn_xm'+Co_xn'*labda(:,n+1));
        grad(Nc,:)          = grad(Nc,:) + J_betan + labda(:,n)'*Cn_betan;
        gradNp(n,:)         = J_betan + labda(:,n)'*Cn_betan;
%     elseif n == Nc          % ONLY if n = Nc
%         labda(:,n)          = - Cn_xn'\(Jn_xm'+Co_xn'*labda(:,n+1));
%         grad(n,:)           = grad(n,:) + J_betan + labda(:,n)'*Cn_betan;
%         gradNp(n,:)         = J_betan + labda(:,n)'*Cn_betan;
    else                    % Everything from 0 to Nc
        labda(:,n)          = - Cn_xn'\(Jn_xm'+Co_xn'*labda(:,n+1));
        grad(n,:)           = J_betan + labda(:,n)'*Cn_betan;
        gradNp(n,:)         = grad(n,:);
    end
    
    Co_xn           = Cn_xm;
        
end
toc

end