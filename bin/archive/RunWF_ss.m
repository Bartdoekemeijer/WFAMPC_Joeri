index = index + 1;

for k=1:Np
    
    it       = 0;
    eps      = 1e19;
    epss     = 1e20;
    
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1; max_it=max_it_dyn; end
        
        % Spatial discretization
        [ax,ay,dax,day] = SpatialDiscr_Hybrid2(Wp,u,v,Turb);
        
        %% Define Ax=b
        
        % Define b
        % 1) Dynamical term
        if k==1; dt = Inf; else dt = h; end % First converge to a steady state solution before running time simultation
                    % Inf for static flow field
        [ax,dax,ay,day,cx,cy,ccx,ccy,dcdx] = Dynamical(Wp,ax,dax,ay,day,u,v,dt);
        
        % 2) Actuator:
        [Sm,dSm,Ueffect(:,k),Ur(:,k),a(:,k),Power(:,k),CT(:,k),CP(:,k),phi(:,k),Fxmean(:,k),Fymean(:,k),Umean(:,k),Vmean(:,k),Umatrix]...
            = Actuator2(Wp,u,v,Phi(:,k),beta(:,k),dbeta(:,k),dPhi(:,k));
        
        % 3) Zero gradient boundary conditions momentum equations
        [ax,ay,dax,day,bx,by,bbx,bby,dbcdx] = BoundaryConditions(Wp,ax,ay,dax,day,u,v);
        
        b       = [bx+cx+c*vec(Sm.x');
            by+cy+c*vec(Sm.y');
            bc            ];
        
        bl      = [c*vec(Sm.dx');
            c*vec(Sm.dy');
            0.*bc            ];
        
        % Define A
        [Ay,~]   =  MakingSparseMatrix(Wp,ay,2,3,1);
        [Ax,~]   =  MakingSparseMatrix(Wp,ax,3,2,1);
        A        = [blkdiag(Ax,Ay) [B1;B2]; [B1;B2]' sparse((Wp.Nx-2)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-2))];
        
        % Define linearised A
        Ayl         = MakingSparseMatrixl(Wp,day,2,3,1);
        Axl         = MakingSparseMatrixl(Wp,dax,3,2,1);
        [Axlo,Aylo] = MakingSparseMatrixlo(Wp,dax,day);
        Al          = [Axl Axlo B1;Aylo Ayl B2;B1' B2' sparse((Wp.Nx-2)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-2))]-A;
        
        
        derivatives = derivatives_calc(A,Al,Ax,Ay,Wp,B1,b,bl,dSm,ccx,ccy,dbcdx);          % Added JW function;
        
        %% Compute solution
        % sol = [u_(i,J) u_(i,J+1)...u_(i,Ny-1) u_(i+1,J) u_(i+1,J+1) ...
        % u_(Nx-1,Ny-1) v_(I,j) v_(I,j+1)...v_(I,Ny-1) v_(I+1,j)...v_(Nx-1,Ny-1) p]
        
        ComputeSolution
        
        % Check if solution converged
        Normv{k}(it) = norm(vec(v(2:end-1,3:end-1)-vv(2:end-1,3:end-1)));
        Normu{k}(it) = norm(vec(u(3:end-1,2:end-1)-uu(3:end-1,2:end-1)));
        eps          = sqrt((Normv{k}(it)+Normu{k}(it)))/((Wp.Ny-2)*(Wp.Nx-2))/2;
        
        display(['k ',num2str(time(k),'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv{k}(it),'%10.2e'), ', Nu=', num2str(Normu{k}(it),'%10.2e'), ', TN=',num2str(eps,'%10.2e'),', Np=','Mean effective=',num2str(mean(Ueffect(1,k)),'%-1000.2f')]) ;
        
        if k==1; alpha      = min(1-.8^it,1); else alpha=1; end;
        % alpha naar 0.9 veranderen? Kijken of dat betere ss oplevert
        
        u(3:end-1,2:end-1)  = (1-alpha)*u(3:end-1,2:end-1)+(alpha)*uu(3:end-1,2:end-1);
        v(2:end-1,3:end-1)  = (1-alpha)*v(2:end-1,3:end-1)+(alpha)*vv(2:end-1,3:end-1);
        p(2:end-1,2:end-1)  = (1-alpha)*p(2:end-1,2:end-1)+(alpha)*pp(2:end-1,2:end-1);
        
        if k==1; ul=u;vl=v;pl=p; end;
        
        % Update velocities for next iteration and boundary conditions
        [u,v,p]    = Updateboundaries(Wp,u,v,p);
        [ul,vl,pl] = Updateboundaries(Wp,ul,vl,pl);
        [du,dv,dp] = Updateboundaries(Wp,du,dv,dp);
        
    end
        
    sol = field2sol(u,v,p,Wp);
    save((strcat('States/state',num2str(index),'_',num2str(k))),'derivatives','sol','u','v','p','du','dv','dp','uu','vv','pp');
    
end
