function Wp = meshing(plotting,Wp,MeshingMethod)
% ddx   = \Delta x_{I,I+1}
% ddy   = \Delta y_{J,J+1}
% ddx2  = \Delta x_{i,i+1}
% ddy2  = \Delta y_{j,j+1}
% dx    = \Delta x_{I-1,I}
% dy    = \Delta y_{J-1,J}
% dx2   = \Delta x_{i-1,i}
% dy2   = \Delta y_{j-1,j}
% ldxx  = I
% ldyy  = J
% ldxx2 = i
% ldyy2 = j

switch lower(Wp.name)
    
    case {'benchmark1_50x25'}
        turbine.MaxCp  = 0.4866;                 % Maximup coefficient
        turbine.Drotor = 90;                     % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   Ly/2;                    % X-coordinate of rotor (center)
        Crx         =   500+2*7*turbine.Drotor;  % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                      % Number of grid points (x-direction)
        Ny          =   25;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid poin
    
        case {'benchmark2_50x25'}
        turbine.MaxCp  = 0.4866;                 % Maximup coefficient
        turbine.Drotor = 90;                     % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [400; 400+7*turbine.Drotor]; % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                      % Number of grid points (x-direction)
        Ny          =   25;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid poin
        
    case {'benchmark2_50x25na'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2+200; Ly/2-200];            % X-coordinate of rotor (center
        Crx         =   [400; 400+7*turbine.Drotor]; % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   100;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
        
    case {'benchmark3_100x50'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [500; 500+7*turbine.Drotor;  500+2*7*turbine.Drotor];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   100;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
        
    case {'benchmark3_50x25'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [500; 500+7*turbine.Drotor;  500+2*7*turbine.Drotor];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                     % Number of grid points (x-direction)
        Ny          =   25;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
        
    case {'benchmark3_20x10'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [500; 500+7*turbine.Drotor;  500+2*7*turbine.Drotor];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   20;                     % Number of grid points (x-direction)
        Ny          =   10;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
        m           =   2;
        lmu         =   2.5;
    case {'benchmark4_50x25'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2; Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [500; 500+7*turbine.Drotor; 500+2*7*turbine.Drotor; 500+3*7*turbine.Drotor];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                     % Number of grid points (x-direction)
        Ny          =   25;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
        
    case {'benchmark6_50x25'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2-3.5*turbine.Drotor; Ly/2+3.5*turbine.Drotor; Ly/2-3.5*turbine.Drotor;...
                            Ly/2+3.5*turbine.Drotor; Ly/2-3.5*turbine.Drotor; Ly/2+3.5*turbine.Drotor];    % X-coordinate of rotor (center)
        Crx         =   [500;500; 500+7*turbine.Drotor;500+7*turbine.Drotor;  500+2*7*turbine.Drotor;500+2*7*turbine.Drotor];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                     % Number of grid points (x-direction)
        Ny          =   25;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
       
        
    case {'benchmark6_50x50'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;                    % Length of the grid in x (N-S direction)
        Ly          =   3000;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2-3.5*turbine.Drotor; Ly/2+3.5*turbine.Drotor; Ly/2-3.5*turbine.Drotor;...
                            Ly/2+3.5*turbine.Drotor; Ly/2-3.5*turbine.Drotor; Ly/2+3.5*turbine.Drotor];    % X-coordinate of rotor (center)
        Crx         =   [500;500; 500+7*turbine.Drotor;500+7*turbine.Drotor;  500+2*7*turbine.Drotor;500+2*7*turbine.Drotor];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
        
    case {'benchmark4_100x50'}
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;            % Length of the grid in x (N-S direction)
        Ly          =   1250;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry         =   [Ly/3; 2*Ly/3; Ly/3; 2*Ly/3];
        Crx         =   [1000; 1000; 1000+8*turbine.Drotor; 1000+8*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs        =   Crx-360;            % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   100;   % Number of grid points (x-direction)
        Ny          =   50;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)\
        
    case {'benchmark4_50x25'}
        turbine.Drotor = 90;            % Rotor diameter
        Lx          =   3000;            % Length of the grid in x (N-S direction)
        Ly          =   1250;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry         =   [Ly/3; 2*Ly/3; Ly/3; 2*Ly/3];
        Crx         =   [1000; 1000; 1000+8*turbine.Drotor; 1000+8*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs        =   Crx-360;            % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   100;   % Number of grid points (x-direction)
        Ny          =   50;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    
case {'benchmark2_50x25na'}
        turbine.MaxCp  = 0.4866;                  % Maximup coefficient
        turbine.Drotor = 90;                      % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1250;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2+200; Ly/2-200];            % X-coordinate of rotor (center
        Crx         =   [400; 400+7*turbine.Drotor]; % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   100;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
    case {'amalia'}
        load centers_Amalia
        load V80_data
        Lx      =   9200;            % Length of the grid in x (N-S direction)
        Ly      =   7500;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   80;             % Rotor diameter
        Crx     =   1820+Centers_turbine(:,2);%[55; 80];             % X-coordinate of rotor (center)
        Crxs     =   Crx-80;%[55; 80];             % X-coordinate of rotor (center)
        %        xlines   =   floor(Crxs/Dx);  % Y grid number of the turbine
        Cry     =   1820+Centers_turbine(:,1);%[20; 200];             % Y-coordinate of rotor (center)
        Qx          =   1;   % change of grid size (x-direction)
        Qy          =   1;   % change of grid size (y-direction)
        Nx          =   300;   % Number of grid points (x-direction)
        Ny          =   200;   % Number of grid points (y-direction)
        sigmax      =   100;   % Number of grid points (x-direction)
        sigmay      =   100;   % Number of grid points (y-direction)
    case {'benchmark6_ieee'}
        load V90_data
        Lx      =   5000;            % Length of the grid in x (N-S direction)
        Ly      =   3000;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry     =   [1000; 1000+8*turbine.Drotor; 1000; 1000+8*turbine.Drotor; 1000; 1000+8*turbine.Drotor];
        Crx     =   [1800; 1800; 1800+10*turbine.Drotor; 1800+10*turbine.Drotor; 1800+20*turbine.Drotor; 1800+20*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs     =   Crx-360;            % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   100;   % Number of grid points (x-direction)
        Ny          =   50;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)\
    case {'benchmark6_yaw'}
        load V90_data
        Lx      =   3800;            % Length of the grid in x (N-S direction)
        Ly      =   2700;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry     =   [1080; 1080+6*turbine.Drotor; 1080-0*turbine.Drotor; 1080+6*turbine.Drotor+0*turbine.Drotor; 1080; 1080+6*turbine.Drotor];
        Crx     =   [1260; 1260; 1260+7*turbine.Drotor; 1260+7*turbine.Drotor; 1260+14*turbine.Drotor; 1260+14*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs     =   Crx-360;            % X-coordinate of rotor (center)
        Qx          =   1;%40;   % change of grid size (x-direction)
        Qy          =   1;%20;   % change of grid size (y-direction)
        Nx          =   150;   % Number of grid points (x-direction)
        Ny          =   150;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_yaw'}
        load V90_data
        Lx      =   2000;            % Length of the grid in x (N-S direction)
        Ly      =   600;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry     =   [300;300-45];             % X-coordinate of rotor (center)
        Crx     =   [500; 500+6*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs     =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   100;   % Number of grid points (x-direction)
        Ny          =   50;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_sowfa'}
        load V90_data
        Lx      =   4500;            % Length of the grid in x (N-S direction)
        Ly      =   2000;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   126;             % Rotor diameter
        Cry     =   [1000; 1000];             % X-coordinate of rotor (center)
        Crx     =   [1910; 2630];             % Y-coordinate of rotor (center)
        Crxs     =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   450;   % Number of grid points (x-direction)
        Ny          =   200;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_sowfa2'}
        load V90_data
        Lx      =   3500;            % Length of the grid in x (N-S direction)
        Ly      =   2000;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   126;             % Rotor diameter
        Cry     =   [1000; 1000];             % X-coordinate of rotor (center)
        Crx     =   [1910; 2630];             % Y-coordinate of rotor (center)
        Crxs     =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   175;   % Number of grid points (x-direction)
        Ny          =   100;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_sowfa3'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 126.3992;            % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1000;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [400;1031.9960];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   100;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)
    case {'benchmark3_sowfa3'} %[0.31; 0.35; 0.49]
        load V90_data
        Lx          =   3000;%6000;            % Length of the grid in x (N-S direction)
        Ly          =   1500;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   126.3992;             % Rotor diameter
        Cry         =   [Ly/2; Ly/2; Ly/2];             % X-coordinate of rotor (center)
        Crx         =   [Lx/2-5*turbine.Drotor; Lx/2; Lx/2+5*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =   1;   % change of grid size (x-direction)
        Qy          =   1;   % change of grid size (y-direction)
        Nx          =   100;%150;  % Number of grid points (x-direction)
        Ny          =   50;%120;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_sowfa4'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 126.3992;            % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1000;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2+250; Ly/2-250];            % X-coordinate of rotor (center)
        Crx         =   [400; 1031.9960];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   100;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)    
    case {'benchmark1_sowfa3'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 126.3992;            % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1000;                    % Length of the grid in y (O-W direction)
        Cry         =   Ly/2;            % X-coordinate of rotor (center)
        Crx         =   1000;         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   100;                     % Number of grid points (x-direction)
        Ny          =   50;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)  
    case {'benchmark2_ieee'}
        load V90_data
        Lx      =   6000;            % Length of the grid in x (N-S direction)
        Ly      =   2000;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry     =   [1000; 1000];             % X-coordinate of rotor (center)
        Crx     =   [1800; 3800];             % Y-coordinate of rotor (center)
        Crxs     =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =  1;% 40;   % change of grid size (x-direction)
        Qy          =   1;%p;20;   % change of grid size (y-direction)
        Nx          =   100;   % Number of grid points (x-direction)
        Ny          =   50;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_coarse'}
        load V90_data
        Lx      =   7000;            % Length of the grid in x (N-S direction)
        Ly      =   2400;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry         =   [1200; 1200];             % X-coordinate of rotor (center)
        Crx         =   [3800; 4800];             % Y-coordinate of rotor (center)
        Crxs        =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =  40;   % change of grid size (x-direction)
        Qy          =  20;   % change of grid size (y-direction)
        Nx          =   60;   % Number of grid points (x-direction)
        Ny          =   30;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark2_coarser'}
        load V90_data
        Lx      =   4000;            % Length of the grid in x (N-S direction)
        Ly      =   1500;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Cry         =   [750; 750];             % X-coordinate of rotor (center)
        Crx         =   [1500; 1500+7*turbine.Drotor];             % Y-coordinate of rotor (center)
        Crxs        =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        Qx          =  40;   % change of grid size (x-direction)
        Qy          =  20;   % change of grid size (y-direction)
        Nx          =   20;   % Number of grid points (x-direction)
        Ny          =   10;   % Number of grid points (y-direction)
        sigmax      =   40;   % Number of grid points (x-direction)
        sigmay      =   80;   % Number of grid points (y-direction)
    case {'benchmark1_ieee'}
        load V90_data
        Lx          =   4400;            % Length of the grid in x (N-S direction)
        Ly          =   2200;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Crx         =   1800;             % X-coordinate of rotor (center)
        Crxs        =   Crx-0.2*turbine.Drotor;%[55; 80];             % X-coordinate of rotor (center)
        Cry         =   1100;             % Y-coordinate of rotor (center)
        Qx          =   20;   % change of grid size (x-direction)
        Qy          =   40;   % change of grid size (y-direction)
        Nx          =   120;   % Number of grid points (x-direction)
        Ny          =   100;   % Number of grid points (y-direction)
        sigmax      =   350;   % Number of grid points (x-direction)
        sigmay      =   400;   % Number of grid points (y-direction)
    case {'benchmark1_close'}
        load V90_data
        Lx          =   1500;            % Length of the grid in x (N-S direction)
        Ly          =   500;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Crx         =   400;             % X-coordinate of rotor (center)
        Crxs        =   Crx-0.2*turbine.Drotor;%[55; 80];             % X-coordinate of rotor (center)
        Cry         =   250;             % Y-coordinate of rotor (center)
        Qx          =   20;   % change of grid size (x-direction)
        Qy          =   40;   % change of grid size (y-direction)
        Nx          =   120;   % Number of grid points (x-direction)
        Ny          =   100;   % Number of grid points (y-direction)
        sigmax      =   350;   % Number of grid points (x-direction)
        sigmay      =   400;   % Number of grid points (y-direction)
    case {'benchmark1_coarse'}
        load V90_data
        Lx      =   1500;            % Length of the grid in x (N-S direction)
        Ly      =   500;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   90;             % Rotor diameter
        Crx         =   400;             % X-coordinate of rotor (center)
        Crxs        =   Crx-0.2*turbine.Drotor;%[55; 80];             % X-coordinate of rotor (center)
        Cry         =   250;             % Y-coordinate of rotor (center)
        Qx          =  20;   % change of grid size (x-direction)
        Qy          =  40;   % change of grid size (y-direction)
        Nx          =   120;   % Number of grid points (x-direction)
        Ny          =   100;   % Number of grid points (y-direction)
        sigmax      =   350;   % Number of grid points (x-direction)
        sigmay      =   400;   % Number of grid points (y-direction)
    case {'benchmark1_minimal'}
        load V90_data
        Lx          =   5;            % Length of the grid in x (N-S direction)
        Ly          =   5;            % Length of the grid in y (O-W direction)
        turbine.Drotor  =   1;        % Rotor diameter
        Cry         =   4;            % X-coordinate of rotor (center) ?
        Crx         =   3;            % Y-coordinate of rotor (center) ?
        Crxs        =   Crx-0.4*turbine.Drotor;%[55; 80]; % X-coordinate of rotor (center)
        Qx          =   1;   % change of grid size (x-direction)
        Qy          =   1;   % change of grid size (y-direction)
        Nx          =   5;   % Number of grid points (x-direction)
        Ny          =   5;   % Number of grid points (y-direction)
        sigmax      =   5;   % Number of grid points (x-direction) exp grid
        sigmay      =   6;   % Number of grid points (y-direction) exp grid   
        turbine.Drotor=1;
        
    case {'benchmark1_simple'}
        load V90_data
        ss          =   20;
        s           =   1;
        Lx          =   ss;                 % Length of the grid in x (N-S direction)
        Ly          =   ss;                 % Length of the grid in y (O-W direction)
        turbine.Drotor  =   2;              % Rotor diameter
        Cry         =   ss/2;               % X-coordinate of rotor (center) ?
        Crx         =   ss/2;               % Y-coordinate of rotor (center) ?
        Crxs        =   Crx-0.4*turbine.Drotor;%[55; 80]; % X-coordinate of rotor (center)
        Qx          =   s;      % change of grid size (x-direction)
        Qy          =   s;      % change of grid size (y-direction)
        Nx          =   ss;     % Number of grid points (x-direction)
        Ny          =   ss;     % Number of grid points (y-direction)
        sigmax      =   ss;     % Number of grid points (x-direction) exp grid
        sigmay      =   ss;     % Number of grid points (y-direction) exp grid
        
    otherwise
        disp('Wind Farm not in list')
end
N=size(Crx,1);

%% For X
% initial grid
Xgrid=0:Lx/5000:Lx;
Xres=zeros(1,length(Xgrid));

for j=1:1:N
    Xres=exp(-(Xgrid-Crx(j)).^2/2/sigmax^2)+Xres;
end
switch lower(MeshingMethod)
    case {'lin'}
        Xres=(sign(Xres+1e-12)+1)./2;
    otherwise
        disp('exponential meshing')
end

XRES=((1+Xres*Qx));

q1 = 0;
for i=1:1:length(XRES)
    q1(i+1)=q1(i)+XRES(i);
end

%plot(q1./q1(end)*Lx,(0:1:length(XRES))/length((XRES))*Lx)
Dist_dxx = interp1(q1./q1(end)*Lx,(0:1:length(XRES))/length(XRES).*Lx,0:Lx/(Nx-1):Lx);
%Dist_dxx = linspace(0,Lx + Lx/(Nx-1), Nx+1);
dvdxx    = diff(Dist_dxx);
dxx      = [dvdxx dvdxx(end)];

dxx2(1)=dxx(1)/2;
for i=1:1:length(dxx)
    if i==length(dxx); else
        dxx2(i+1)=dxx(i)/2+dxx(i+1)/2;
    end
end
%% For Y
% initial grid
Ygrid=0:Ly/5000:Ly;
Yres=zeros(1,length(Ygrid));

for j=1:1:N
    Yres=exp(-(Ygrid-Cry(j)).^2/2/sigmay^2)+Yres;
end
switch lower(MeshingMethod)
    case {'lin'}
        Yres=(sign(Yres+1e-12)+1)./2;
    otherwise
        disp(' ')
end
YRES=((1+Yres*Qy));

q1=[0];
for i=1:1:length(YRES)
    q1(i+1)=q1(i)+YRES(i);
end

%plot(q1,(0:1:length(XRES))/length((XRES)*Lx)
Dist_dyy = interp1(q1./q1(end)*Ly,(0:1:length(YRES))/length(YRES).*Ly,0:Ly/(Ny-1):Ly);
%Dist_dyy = linspace(0,Ly + Ly/(Ny-1), Ny+1);
dvdyy    = diff(Dist_dyy);
dyy      = [dvdyy dvdyy(end)];

dyy2(1)=dyy(1)/2;
for i=1:1:length(dyy)
    if i==length(dyy); else
        dyy2(i+1)=dyy(i)/2+dyy(i+1)/2;
    end
end

%%
dxx=dxx'*ones(1,Ny);
dxx2=dxx2'*ones(1,Ny);
dyy=ones(Nx,1)*dyy;
dyy2=ones(Nx,1)*dyy2;

% ldxx=Dist_dxx(1:end-1)'*ones(1,Ny);
% ldyy=ones(Nx,1)*Dist_dyy(1:end-1);
%
% Dist_dxx2=Dist_dxx(1:end-1)+diff(Dist_dxx)-0.5*Dist_dxx(1,2);
% Dist_dyy2=Dist_dyy(1:end-1)+diff(Dist_dyy)-0.5*Dist_dyy(1,2);

ldxx=Dist_dxx'*ones(1,Ny);
ldyy=ones(Nx,1)*Dist_dyy;

dvdxx = diff(Dist_dxx);
dvdyy = diff(Dist_dyy);

Dist_dxx2=Dist_dxx+[dvdxx dvdxx(end)]-0.5*Dist_dxx(1,2);
Dist_dyy2=Dist_dyy+[dvdyy dvdyy(end)]-0.5*Dist_dyy(1,2);

ldxx2=Dist_dxx2(1:end)'*ones(1,Ny);
ldyy2=ones(Nx,1)*Dist_dyy2(1:end);

if plotting==1;
    figure(20)
    axis([0 Ly 0 Lx])
    vline(Dist_dyy,'k')
    hline(Dist_dxx,'k')
    xlabel('y');ylabel('x');
end

for i=1:1:N
    [MINNIE,xline(i,:)]   =   min((Crx(i)-ldxx2(:,1)).^2);%floor(Crx/Dx);  % Y grid number of the turbine
end
for i=1:1:N
    [MINNIE,xlines(i,:)]   =   min((Crxs(i)-(ldxx2(:,1))).^2);%floor(Crx/Dx);  % Y grid number of the turbine
end

for i=1:1:N
    [Minnie,Q1]=min(((Cry(i)-0.5*turbine.Drotor-ldyy2(1,:)')).^2);
    [Minnie,Q2]=min(((Cry(i)+0.5*turbine.Drotor-ldyy2(1,:)')).^2);
    
    yline{i}   =   floor(Q1):1:floor(Q2);
    ylinev{i}   =   floor(Q1):1:floor(Q2)+1; %added JW
  %  ylineu{i}   =   floor(Q1)-1:1:floor(Q2)+1; %added JW
end

Wp.dxx     = dxx;
Wp.dyy     = dyy;
Wp.dxx2    = dxx2;
Wp.dyy2    = dyy2;
Wp.ldxx    = ldxx;
Wp.ldyy    = ldyy;
Wp.ldxx2   = ldxx2;
Wp.ldyy2   = ldyy2;
Wp.Nx      = Nx;
Wp.Ny      = Ny;
Wp.Lx      = Lx;
Wp.Ly      = Ly;
Wp.xline   = xline;
Wp.yline   = yline;
Wp.ylinev   = ylinev;
%  Wp.ylineu   = ylineu;
Wp.xlines  = xlines;
Wp.N       = N;
Wp.turbine = turbine;
Wp.Crx     = Crx;
Wp.Cry     = Cry;
Wp.Crxs    = Crxs;
Wp.lmu     = lmu;
Wp.m       = m;

% meshgrid(Wp.ldyy(1,1:end),Wp.ldxx(1:end,1))

