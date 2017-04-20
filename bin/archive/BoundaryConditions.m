function [ax,ay,dax,day,bx,by,bbx,bby,dbcdx] = BoundaryConditions(Wp,ax,ay,dax,day,u,v)

Nx = Wp.Nx;
Ny = Wp.Ny;

bbx = sparse((Nx-3)*(Ny-2),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));
bby = sparse((Nx-2)*(Ny-3),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));

% Zero gradient outflow
% for u-direction
ax.aP((Nx-1),(2:Ny-1))  = ax.aP((Nx-1),(2:Ny-1)) - ax.aE((Nx-1),(2:Ny-1)); %NORTH
ax.aP((3:Nx-1),(Ny-1))  = ax.aP((3:Nx-1),(Ny-1)) - ax.aN((3:Nx-1),(Ny-1)); %EAST
ax.aP((3:Nx-1),(2))     = ax.aP((3:Nx-1),2)      - ax.aS((3:Nx-1),2);

dax.P((Nx-1),(2:Ny-1))  = dax.P((Nx-1),(2:Ny-1)) - dax.E((Nx-1),(2:Ny-1)); %NORTH
dax.P((3:Nx-1),(Ny-1))  = dax.P((3:Nx-1),(Ny-1)) - dax.N((3:Nx-1),(Ny-1)); 
dax.P((3:Nx-1),(2))     = dax.P((3:Nx-1),2)      - dax.S((3:Nx-1),2);

%% inflow conditions (no zero gradient this seems to work
% dax.NW(:,2)        = dax.NW(:,2)+  dax.aPN(:,2).*u(:,2); % a bit random (i can not explain)
% dax.NE(1:Nx,2)     = dax.NE(1:Nx,2)  + dax.aPN(1:Nx,2).*u(1:Nx,2); % a bit random (i can not explain)
% 
% dax.SW(:,end-1)    = dax.SW(:,end-1)-  dax.aS(1:Nx,end-1).*u(1:Nx,end-1) ; % a bit random (i can not explain)
% dax.SE(1:Nx,end-1) = dax.SE(1:Nx,end-1) - dax.aS(1:Nx,end-1).*u(1:Nx,end-1); % a bit random (i can not explain)
% 
% 
% dax.P(2:3,:)       = dax.P(2:3,:)+dax.aW(2:3,1:Ny).*u(1:2,1:Ny); % y momentum west side
% 
% 
% 
% day.SW(end-1,:)    = day.SW(end-1,:)-  day.aW(end-1,:).*v(end-1,:) ; % a bit random (i can not explain)
% day.NW(end-1,:)    = day.NW(end-1,:)-  day.aW(end-1,:).*v(end-1,:) ; % a bit random (i can not explain)

%dax.NW(2:3,1:Ny-1) = dax.NW(2:3,1:Ny-1)-dax.aS(2:3,1:Ny-1).*u(2:3,2:Ny) ; %JW

% Sjoerd found the following

dax.NW(:,2)        = dax.NW(:,2)    -  dax.SW(:,2); 
dax.NE(:,2)        = dax.NE(:,2)    -  dax.SE(:,2); 

dax.SW(Nx-1,:)     = dax.SW(Nx-1,:) -  dax.NW(Nx-1,:); 
dax.SE(Nx-1,:)     = dax.SE(Nx-1,:) -  dax.NE(Nx-1,:);

day.SE(Nx-1,:)     = day.SE(Nx-1,:) -  day.SW(Nx-1,:);
day.NE(Nx-1,:)     = day.NE(Nx-1,:) -  day.SW(Nx-1,:);


%% for v-direction
ay.aP((Nx-1),(1:Ny))  = ay.aP((Nx-1),(1:Ny)) - ay.aE((Nx-1),(1:Ny));
ay.aP((1:Nx),(Ny-1))  = ay.aP((1:Nx),(Ny-1)) - ay.aN((1:Nx),(Ny-1));
ay.aP((1:Nx),(3))     = ay.aP((1:Nx),3)      - ay.aS((1:Nx),3); % changed to 3 3 2 instead of 2 2 1

day.P((Nx-1),(1:Ny))  = day.P((Nx-1),(1:Ny)) - day.E((Nx-1),(1:Ny));
day.P((1:Nx),(Ny-1))  = day.P((1:Nx),(Ny-1)) - day.N((1:Nx),(Ny-1));
day.P((1:Nx),(3))     = day.P((1:Nx),3)      - day.S((1:Nx),3);% changed to 3 3 2 instead of2 2 1


% Inflow boundary for non linear model
bx      = kron([1;zeros(Nx-4,1)],(ax.aW(3,2:end-1).*u(2,2:end-1))'); %changed to 3: 2 instead of 2:2
by      = [v(1,3:Ny-1)'.*ay.aW(2,3:Ny-1)';zeros((Nx-3)*(Ny-3),1)]; %changed to 2:3 inst

% Inflow boundary for linear model
bbx(1:Ny-2,1:Ny-2) = diag((dax.aW(3,2:end-1).*u(2,2:end-1))');
bby(1:Ny-2,1:Ny-2) = diag((day.aW(2,2:end-1).*v(1,2:end-1))');

dbcdx = [bbx;bby;sparse((Nx-2)*(Ny-2),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2))];

end