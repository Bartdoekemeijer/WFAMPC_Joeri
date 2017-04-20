function [ax,dax,ay,day,cx,cy,ccx,ccy,dcdx] = Dynamical(Wp,ax,dax,ay,day,u,v,dt)
global Rho

% Fully implicit (page 248 Versteeg) See also page 257
ax.aP   = ax.aP + Rho.*Wp.dxx.*Wp.dyy2./dt;       % Rho.*dxx.*dyy2./dt = a_P^0
ay.aP   = ay.aP + Rho.*Wp.dxx2.*Wp.dyy./dt;

dax.P   = dax.P + Rho.*Wp.dxx.*Wp.dyy2./dt;       
day.P   = day.P + Rho.*Wp.dxx2.*Wp.dyy./dt;

% Use u and v from the previous time step
cx     = vec(Rho.*Wp.dxx(3:end-1,2:end-1)'.*Wp.dyy2(3:end-1,2:end-1)'./dt).*vec(u(3:end-1,2:end-1)');
cy     = vec(Rho.*Wp.dxx2(2:end-1,3:end-1)'.*Wp.dyy(2:end-1,3:end-1)'./dt).*vec(v(2:end-1,3:end-1)');
ccx    = vec(Rho.*Wp.dxx(3:end-1,2:end-1)'.*Wp.dyy2(3:end-1,2:end-1)'./dt);
ccy    = vec(Rho.*Wp.dxx2(2:end-1,3:end-1)'.*Wp.dyy(2:end-1,3:end-1)'./dt);

dcdx   = blkdiag(diag(ccx),diag(ccy),sparse((Wp.Ny-2)*(Wp.Nx-2),(Wp.Ny-2)*(Wp.Nx-2)));



