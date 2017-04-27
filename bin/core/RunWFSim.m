index                   = index + 1;

for k=1:options.AMPC.Np
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1
            max_it = max_it_dyn;
        end
        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);                      % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                      % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);            % Map solution to field
        
    end
    
    derivatives = sys.derivatives;
    x           = sol.x;
    save((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(k))),'sys','x');    
end;