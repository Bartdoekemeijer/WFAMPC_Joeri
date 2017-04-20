if  Projection==1;
    
    Ct                      = blkdiag(spdiags(ccx,0,length(ccx),length(ccx)),spdiags(ccy,0,length(ccy),length(ccy)));
    Et                      = Qsp'*blkdiag(Ax,Ay)*Qsp;
    
    if (it==1 && k==1);pRCM = symrcm(Et);end
    
    At                      = Qsp'*Ct*Qsp;
    St                      = Qsp'*[bx;by] - Qsp'*blkdiag(Ax,Ay)*Bsp + Qsp'*Ct*Bsp;
    Bt                      = Qsp'*[Sm.xx;Sm.yy];
    Ft                      = At*solnew + Bt*beta(:,k) + St;
    
    solnew(pRCM,1)          = Et(pRCM,pRCM)\Ft(pRCM);
    sol                     = Qsp*solnew + Bsp;
    uu(3:end-1,2:end-1)     = reshape(sol(1:(Wp.Nx-3)*(Wp.Ny-2)),Wp.Ny-2,Wp.Nx-3)';
    vv(2:end-1,3:end-1)     = reshape(sol((Wp.Nx-3)*(Wp.Ny-2)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)),Wp.Ny-3,Wp.Nx-2)';
    
    % Linear version
        
    if k==1; sollnew = 0*solnew;ul=u;vl=v;pl=p;end
    
    Atl = dcdx+dbcdx+dSm.dx-Al;
    Atl = Qsp'*Atl(1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3),1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3))*Qsp;
    Etl = Et;
    Btl = Qsp'*[Sm.dxx;Sm.dyy];
    Ftl = Atl*sollnew + Btl*dbeta(:,k);
    
    sollnew(pRCM,1)          = Etl(pRCM,pRCM)\Ftl(pRCM);
    soll                     = Qsp*sollnew;
    
    du(3:end-1,2:end-1)     = reshape(soll(1:(Wp.Nx-3)*(Wp.Ny-2)),Wp.Ny-2,Wp.Nx-3)';
    dv(2:end-1,3:end-1)     = reshape(soll((Wp.Nx-3)*(Wp.Ny-2)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)),Wp.Ny-3,Wp.Nx-2)';
    
    ul(3:end-1,2:end-1)     = ul(3:end-1,2:end-1)+du(3:end-1,2:end-1);
    vl(2:end-1,3:end-1)     = vl(2:end-1,3:end-1)+dv(2:end-1,3:end-1);
    
    
else
    % Remove zero columns and rows from A and b
    A(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    b(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    A(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = [];
    
    A(:,end) = [];A(end,:) = [];b(end)=[];
    
    if (it==1 && k==1);pRCM = symrcm(A);end
    
    sol(pRCM,1)             = A(pRCM,pRCM)\b(pRCM);
    
    uu(3:end-1,2:end-1)     = reshape(sol(1:(Wp.Nx-3)*(Wp.Ny-2)),Wp.Ny-2,Wp.Nx-3)';
    vv(2:end-1,3:end-1)     = reshape(sol((Wp.Nx-3)*(Wp.Ny-2)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)),Wp.Ny-3,Wp.Nx-2)';
    %pp(2:end-1,2:end-1)     = reshape([sol((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+1:end);0],Wp.Ny-2,Wp.Nx-2)';
    pp(2:end-1,2:end-1)     = reshape([sol((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+1:end);0;0],Wp.Ny-2,Wp.Nx-2)';
    pp(isinf(pp))           = 0;
    
    
    % Linear version
    All        = (dcdx+dbcdx+dSm.dx-Al);
    
    All(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    bl(size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1,:) = [];
    All(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Wp.Ny-2)+1) = [];
    
    All(:,end) = [];All(end,:) = [];bl(end) = [];
    
    
    if k==1; soll = zeros(size(All,2),1);ul=u;vl=v;pl=p;end
    
    bll                     = All*soll+bl;
    soll(pRCM,1)            = A(pRCM,pRCM)\bll(pRCM);
    
    du(3:end-1,2:end-1)     = reshape(soll(1:(Wp.Nx-3)*(Wp.Ny-2)),Wp.Ny-2,Wp.Nx-3)';
    dv(2:end-1,3:end-1)     = reshape(soll((Wp.Nx-3)*(Wp.Ny-2)+1:(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)),Wp.Ny-3,Wp.Nx-2)';
    %dp(2:end-1,2:end-1)     = reshape([soll((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+1:end);0],Wp.Ny-2,Wp.Nx-2)';
    dp(2:end-1,2:end-1)     = reshape([soll((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+1:end);0;0],Wp.Ny-2,Wp.Nx-2)';
    dp(isinf(dp))           = 0;
    
    ul(3:end-1,2:end-1)     = ul(3:end-1,2:end-1)+du(3:end-1,2:end-1);
    vl(2:end-1,3:end-1)     = vl(2:end-1,3:end-1)+dv(2:end-1,3:end-1);
    %pl(2:end-1,2:end-1)     = pl(2:end-1,2:end-1)+dp(2:end-1,2:end-1);
    pl(2:end-1,2:end-1)     = reshape([soll((Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+1:end);0;0],Wp.Ny-2,Wp.Nx-2)';
    pl(isinf(pl))           = 0;
    
end