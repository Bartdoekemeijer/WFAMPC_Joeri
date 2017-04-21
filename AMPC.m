clear; clc; close all
addpath WFSim\libraries\sparse_null
addpath WFSim\bin\core                      % WFSim model directory
addpath bin\core                            % WFAMPC directory
%addpath bin\archive                        % WFAMPC directory (old needs to be replaced)

Wp.name             = 'ThreeTurbine_Ampc';
Wp.Turbulencemodel  = 'WFSim3';

%% Init
WFAMPC_initialize


%% Control inputs (has to go in initialize)
b0              = 0.5;%[0.26;0.1;0.54];
Phi             = zeros(Wp.N,Np);       % Yaw angles in degrees (-90 < Phi < 90 degrees)
if length(b0) > 1
    beta        = b0*ones(1,Np);        % Scaled axial induction 0 < beta < 1
    beta0       = b0;
else
    beta        = b0*ones(Wp.N,Np);     % Scaled axial induction 0 < beta < 1
    beta0       = b0*ones(Wp.N,1);
end
dPhi            = zeros(Wp.N,Np);       % Yaw angles in degrees (-90 < Phi < 90 degrees), linear model
dbeta           = 0*ones(Wp.N,Np);      % Scaled axial induction (0 < beta < 1), linear model
phi             = zeros(Wp.N,Np);   	% Wind angles in degrees

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end

%% Simulate Wind Farm towards Steady State

tic
index = 0;
RunWF_ss;  % index 1
toc

%% First forward simulation + Backwards adjoint

load((strcat('States/state',num2str(index),'_',num2str(Nr))))
% index = index + 1;
tic
RunWF % index 2
toc
tic
[grad,labda,J_partial]   = get_adjoint_Nc(Np,Nc,sol,Wp,index);
toc
grad 	= gradproj(grad,beta,beta_lim,eps);
dJmax   = max(max(abs(grad)))';

%% Gradient calculation

stop_ls         = 0;
grad            = get_adjoint_Nc(Np,Nc,sol,Wp,index);
dJmaxinit       = max(max(abs(grad)))';
alphai          = 1/(dJmaxinit);
alphai0         = alphai;
BETA(:,1:Nr)    = beta(:,1:Nr);
GRAD(:,1:Nr)    = grad(1:Nr,:)';
POWER(:,1:Nr)   = Power(:,1:Nr);
POWERTOT(:,1)   = sum(Power,2);
J(1)            = sum(POWERTOT(:,1));
P               = zeros(Nrmax,imax+1);
P2              = zeros(Nrmax,imax+2);
ALPHA           = [];
dP_normi        = dP_norm*ones(1,Nrmax);
constant        = 0;

index = 1;
    
for i = 1:Nrmax
    
    % Forward simulation 1 -> Np
    load((strcat('States/state',num2str(index),'_',num2str(Nr))))
    if constant == 0
        beta        = [beta(:,Nr+1:end) beta(:,end)*ones(1,Nr)];
        beta0       = beta(:,Nr);
    end
    constant    = 0;
    indexNr     = index;
    RunWF
    
    % Backward adjoint Np -> 1
    grad        = get_adjoint_Nc(Np,Nc,sol,Wp,index);
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
    
    while ls <= imax && ~(check1 && check2) && dP(i) > dP_normi(i)
        
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
        beta    = update_beta(beta_prev,grad,beta_lim,Nc,dbeta_max,alphai,beta0);

        load((strcat('States/state',num2str(indexNr),'_',num2str(Nr))))
        RunWF
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
%         beta        = beta(:,1)*ones(1,Np);
        beta(:,1:Nr)    = beta(:,1)*ones(1,Nr);
%         dP_normi(i+1)   = dP_norm*2;
        load((strcat('States/state',num2str(indexNr),'_',num2str(Nr))))
        RunWF
        P(i,2)      = sum(sum(Power));
        P2(i,3)     = P(i,2);
        constant    = 1;
        alphai      = alphai0;
        display(['i =',num2str(i),', constant Power, Cost J =',num2str(P(i,1),'%10.4e')]);
    elseif P(i,ls) < P(i,ls-1) && P(i,ls-1) > P(i,1)
        beta    = beta_prev2;
        Power   = Power_prev2;
        index   = index-1;
        stopls  = 1;
    elseif P(i,ls) < P(i,1)
        beta    = beta_prev;
        Power   = Power_prev;
        index   = indexNr;
    end
    
    % Save intermediate results
    BETA(:,1+Nr*(i-1):Nr*i)     = beta(:,1:Nr);
    GRAD(:,1+Nr*(i-1):Nr*i)     = grad(1:Nr,:)';
    POWER(:,1+Nr*(i-1):Nr*i)    = Power(:,1:Nr);
    POWERTOT(:,i)               = sum(Power,2);
    J(i)                        = sum(POWERTOT(:,i));
    
    if ShowGrad == 1
        figure;
        plot(grad);
        figure;
        plot(beta');
    end
    
    if ChangeDir == 1 && i == round(Nrmax/2)
        v_Inf   = u_Inf*sin(angleDir);
        u_Inf   = u_Inf*cos(angleDir);
        [B1,B2,Bm1,Bm2,bc] = Compute_B1_B2_bc(Wp,u_Inf);
        if Projection==1
            [Qsp, Bsp] = Solution_space(B1,B2,bc);
            solnew = Qsp\([vec(u(3:end-1,2:end-1)');vec(v(2:end-1,3:end-1)')]-Bsp); % Initial condition
        end
        if Animate==1;scrsz = get(0,'ScreenSize');figure('color',[0 166/255 214/255],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)],...
                'ToolBar','none','visible', 'on');end
    end
    
    if ChangeSpeed == 1 && i == round(Nrmax/2)
        u_Inf   = 10;
        [B1,B2,Bm1,Bm2,bc] = Compute_B1_B2_bc(Wp,u_Inf);
        if Projection==1
            [Qsp, Bsp] = Solution_space(B1,B2,bc);
            solnew = Qsp\([vec(u(3:end-1,2:end-1)');vec(v(2:end-1,3:end-1)')]-Bsp); % Initial condition
        end
        if Animate==1;scrsz = get(0,'ScreenSize');figure('color',[0 166/255 214/255],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)],...
                'ToolBar','none','visible', 'on');end
    end
    
end
    
%% Figures: opmaak verbeteren!

figure;hold on;plot(h:h:h*Nrmax*Nr,POWER')
plot(h:h:h*Nrmax*Nr,sum(POWER))
if ChangeSpeed == 1
    load(strcat('greedypower',num2str(Wp.N),'T_u8_u10'))
    plot(h:h:h*Nrmax*Nr,greedyP_u8_u10,'--')
elseif ChangeDir == 1
    load(strcat('greedypower',num2str(Wp.N),'T',num2str(angleDir*180/pi)))
    plot(h:h:h*Nrmax*Nr,greedyP,'--')
else
    load(strcat('greedypower',num2str(Wp.N),'T'))
    plot(h:h:h*Nrmax*Nr,greedyP_u8*ones(1,Nrmax*Nr),'--')
%     plot(h:h:h*Nrmax*Nr,[greedyP_u8*ones(1,round(Nrmax/2)*Nr) greedyP_u10*ones(1,round(Nrmax/2)*Nr)],'--')
end
grid on
xlabel('Time [s]')
ylabel('Power [W]')
legend('Turbine 1','Turbine 2','Turbine 3','Total power','Greedy control power')

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

delete('States/*.mat')
