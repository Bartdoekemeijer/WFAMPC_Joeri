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
options.AMPC.ShowGrad     = 0;                   % Show plot of gradient for each simulation
options.AMPC.ShowBeta     = 0;                   % Show plot of beta for each simulation
options.AMPC.method       = 'grad_ratio';        % Control method used

options.AMPC.beta_lim        = [0.1;0.9];        % Limits of input beta
options.AMPC.dbeta_max       = 0.1;              % Limit on change of beta per time step
options.AMPC.Np              = 200;              % Prediction Horizon (in time steps)
options.AMPC.Nc              = 200;              % Control horizon (in time steps)
options.AMPC.Nr              = 10;               % Receding horizon (in time steps)
options.AMPC.gamma           = 5e-7;             % For Steepest-descent method
options.AMPC.imax            = 10;               % Max number of iterations within timestep
options.AMPC.iter_ls         = 1;
options.AMPC.pen             = 1e4*ones(1,options.AMPC.Nc);   % penalty on dbeta
options.AMPC.labda_end       = 100;              % end-point penalty
options.AMPC.eps_grad        = 0.01;
options.AMPC.dP_norm         = 0.002;
options.AMPC.Nrmax           = 30;               % Number of receding horizons
options.AMPC.counter         = 0;
