function sol = field2sol(u,v,p,Wp)

clear sol pex
sol(1:(Wp.Nx-3)*(Wp.Ny-2))=vec(u(3:Wp.Nx-1,2:Wp.Ny-1)');
sol((Wp.Nx-3)*(Wp.Ny-2)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3))=vec(v(2:Wp.Nx-1,3:Wp.Ny-1)');
pex=vec(p(2:end-1,2:end-1)');
sol((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2)-2)=pex(1:end-2);
sol=sol';

end