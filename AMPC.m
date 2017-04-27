clear; clc; close all
addpath WFSim\libraries\sparse_null
addpath WFSim\bin\core                      % WFSim model directory
addpath bin\core                            % WFAMPC directory
%addpath bin\archive                        % WFAMPC directory (old needs to be replaced)

Wp.name             = 'ThreeTurbine_Ampc';
Wp.Turbulencemodel  = 'WFSim3';

%% Init
WFAMPC_initialize

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end

%% Simulate Wind Farm towards Steady State
index                   = 0;
options.startUniform    = 1;    % Start from a uniform flowfield (true) or a steady-state solution (false)
max_it                  = 50;  

RunWFSim;                          % index 1 (go to steady state)

%% First forward simulation + Backwards adjoint
load((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(options.AMPC.Nr))))

options.startUniform    = 1;    % Start from a uniform flowfield (true) or a steady-state solution (false)
max_it                  = 1;  

RunWFSim                        % index 2 (start from steady state)

[grad,labda,J_partial]   = get_adjoint_Nc(options.AMPC,x,Wp,index);

% This is temporarily, time and inputs need to be redefined in meshing
for kk=1:options.AMPC.Nc  
    beta(:,kk) = input{kk}.beta;
end

grad 	                 = gradproj(grad,beta,options.AMPC.beta_lim,eps);
dJmax                    = max(max(abs(grad)))';

%% Gradient calculation (Can this init not be done above, at least some can)
stop_ls                     = 0;
grad                        = get_adjoint_Nc(options.AMPC,x,Wp,index);
dJmaxinit                   = max(max(abs(grad)))';
alphai                      = 1/(dJmaxinit);
alphai0                     = alphai;
BETA(:,1:options.AMPC.Nr)   = beta(:,1:options.AMPC.Nr);
GRAD(:,1:options.AMPC.Nr)   = grad(1:options.AMPC.Nr,:)';
POWER(:,1:options.AMPC.Nr)  = Power(:,1:options.AMPC.Nr);
POWERTOT(:,1)               = sum(Power,2);
J(1)                        = sum(POWERTOT(:,1));
P                           = zeros(options.AMPC.Nrmax,options.AMPC.imax+1);
P2                          = zeros(options.AMPC.Nrmax,options.AMPC.imax+2);
ALPHA                       = [];
dP_normi                    = options.AMPC.dP_norm*ones(1,options.AMPC.Nrmax);
constant                    = 0;

index = 1;
    
for i = 1:options.AMPC.Nrmax
    
    % Forward simulation 1 -> Np
    load((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(options.AMPC.Nr))))
    if constant == 0
        beta        = [beta(:,options.AMPC.Nr+1:end) beta(:,end)*ones(1,options.AMPC.Nr)];
        beta0       = beta(:,options.AMPC.Nr);
    end
    
    % Write here an update of input.beta
    for kk=1:size(beta,2)
        input{kk}.beta = beta(:,kk);
    end
    
    constant    = 0;
    indexNr     = index;
    RunWFSim
    
    % Backward adjoint Np -> 1
    grad        = get_adjoint_Nc(options.AMPC,x,Wp,index);
%     grad        = gradproj(grad,beta,beta_lim,eps);
    beta_prev   = beta;
    Power_prev  = Power;
    ls          = 1;
    check1      = 0;
    check2      = 0;
    alphai      = alphai*8; % Een keer *8 proberen?
    P(i,ls)     = sum(sum(Power));
    P2(i,ls)    = P(i,ls);
    display(['i =',num2str(i),', Initial cost J =',num2str(P(i,ls),'%10.4e')])
    
    P2(i,ls+1)  = P(i,ls);
    
    if i > 1
        dP(i)       = abs((P(i,1)/P(i-1,1))-1);
    else
        dP(1)       = dP_normi(1) + 1;
    end
    
    while ls <= options.AMPC.imax && ~(check1 && check2) && dP(i) > dP_normi(i)
        
        beta_prev2  = beta;
        Power_prev2 = Power;
        check1      = 0;
        check2      = 0;
             
        % Determine new beta
%         if P(i,ls) <= P(i,1)
%             alphai  = alphai/2;
%             ALPHA(length(ALPHA)+1) = alphai;
%         elseif P2(i,ls+1) > P2(i,ls)
%             alphai  = alphai*2;
%             ALPHA(length(ALPHA)+1) = alphai;
%         end
        alphai  = alphai/2;
        ALPHA(length(ALPHA)+1) = alphai;
        
%         beta    = update_beta_abs(beta_prev,grad,beta_lim,Nc,dbeta_max,alphai);
        beta    = update_beta(beta_prev,grad,alphai,beta0,options);
        
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        load((strcat('data_WFAMpc/states/state',num2str(indexNr),'_',num2str(options.AMPC.Nr))))
        RunWFSim
        ls          = ls + 1;
        P(i,ls)     = sum(sum(Power));
        P2(i,ls+1)  = P(i,ls);
        display(['i =',num2str(i),', line search ',num2str(ls-1),', Cost J =',num2str(P(i,ls),'%10.4e')])
        
        if P(i,ls) > P(i,1) || P(i,ls-1) > P(i,1)
            check1  = 1;
        end
        if P(i,ls) < P(i,ls-1)
            check2  = 1;
        end
        
    end
    
    if dP(i) < dP_normi(i)
%         beta        = beta(:,1)*ones(1,option.AMPC.Np);
        beta(:,1:options.AMPC.Nr)    = beta(:,1)*ones(1,options.AMPC.Nr);
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        %         dP_normi(i+1)   = dP_norm*2;
        load((strcat('data_WFAMpc/states/state',num2str(indexNr),'_',num2str(options.AMPC.Nr))))
        RunWFSim
        P(i,2)      = sum(sum(Power));
        P2(i,3)     = P(i,2);
        constant    = 1;
        alphai      = alphai0;
        display(['i =',num2str(i),', constant Power, Cost J =',num2str(P(i,1),'%10.4e')]);
    elseif P(i,ls) < P(i,ls-1) && P(i,ls-1) > P(i,1)
        beta    = beta_prev2;
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        Power   = Power_prev2;
        index   = index-1;
        stopls  = 1;
    elseif P(i,ls) < P(i,1)
        beta    = beta_prev;
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        Power   = Power_prev;
        index   = indexNr;
    end
    
    % Save intermediate results
    BETA(:,1+options.AMPC.Nr*(i-1):options.AMPC.Nr*i)   = beta(:,1:options.AMPC.Nr);
    GRAD(:,1+options.AMPC.Nr*(i-1):options.AMPC.Nr*i)   = grad(1:options.AMPC.Nr,:)';
    POWER(:,1+options.AMPC.Nr*(i-1):options.AMPC.Nr*i)  = Power(:,1:options.AMPC.Nr);
    POWERTOT(:,i)                                       = sum(Power,2);
    J(i)                                                = sum(POWERTOT(:,i));
    
    if options.AMPC.ShowGrad == 1
        figure;
        plot(grad);
        figure;
        plot(beta');
    end
    
end
    
%% Figures: opmaak verbeteren!
h           = Wp.sim.h;
Nrmax       = options.AMPC.Nrmax;
Nr          = options.AMPC.Nr;

figure;hold on;plot(h:h:h*Nrmax*Nr,POWER')
plot(h:h:h*Nrmax*Nr,sum(POWER))


grid on
xlabel('Time [s]')
ylabel('Power [W]')
legend('Turbine 1','Turbine 2','Turbine 3','Total power')

figure;hold on;plot(h:h:h*Nrmax*Nr,BETA')
plot(h:h:h*Nrmax*Nr,[0.26;0.1;0.54]*ones(1,Nrmax*Nr),'--')
xlabel('Time [s]')
ylabel('\beta [-]')
grid on
legend('Turbine 1','Turbine 2','Turbine 3','Optimum T_1','Optimum T_2','Optimum T_3')

figure; plot(h:h:h*Nrmax*Nr,GRAD); hold on
plot(h:h:h*Nrmax*Nr,zeros(1,Nrmax*Nr),'--')
xlabel('Time [s]')
grid on
ylabel('\nabla_{\beta} J')

figure; plot(h:h:h*Nrmax,POWERTOT');hold on
plot(h:h:h*Nrmax,J)
xlabel('Time [s]')
ylabel('Power [W]')
legend('Turbine 1','Turbine 2','Turbine 3','Total power')
grid on

delete('data_WFAMpc/states/*.mat')
