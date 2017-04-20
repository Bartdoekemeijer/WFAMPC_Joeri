%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'WFAMPC_initialize.m'
%  This script loads the model and MPC settings. It also prepares
%  the meshing and all variables required for simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WFSim settings in WFAMPC
Animate       = 0;                   % Show 2D flow fields
conv_eps      = 1e-6;                % Convergence threshold
max_it_dyn    = 1;                   % Maximum number of iterations for k > 1

options.Projection     = 0;                      % Use projection (true/false)
options.Linearversion  = 1;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 1;                      % Calculate linear solution of WFSim
options.Derivatives    = 1;                      % Compute derivatives
options.exportPressures= ~options.Projection;    % Calculate pressure fields

[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,0);

% MPC settings in WFAMPC
ShowGrad     = 0;                   % Show plot of gradient for each simulation
ShowBeta     = 0;                   % Show plot of beta for each simulation

method          = 'grad_ratio';     % Control method used
beta_lim        = [0.1;0.9];        % Limits of input beta
dbeta_max       = 0.1;              % Limit on change of beta per time step
Np              = 201;              % Prediction Horizon (in time steps)
Nc              = 201;              % Control horizon (in time steps)
Nr              = 10;               % Receding horizon (in time steps)
gamma           = 5e-7;             % For Steepest-descent method
imax            = 10;               % Max number of iterations within timestep
iter_ls         = 1;
pen             = 1e4*ones(1,Nc);   % penalty on dbeta
labda_end       = 100;              % end-point penalty
eps_grad        = 0.01;
dP_norm         = 0.002;
ChangeDir       = 0;
angleDir        = 30*(pi/180);
ChangeSpeed     = 0;
Nrmax           = 30;               % Number of receding horizons
counter         = 0;
