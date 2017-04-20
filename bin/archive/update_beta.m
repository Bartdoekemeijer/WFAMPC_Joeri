function beta = update_beta(beta,grad,beta_lim,Nc,dbeta_max,alphai,beta0)

% Find new beta
dJ      = grad';
s       = sign(dJ);
betas   = beta(:,1:Nc);
betas   = betas-((beta_lim(2)-betas).*(s==-1)+...
            (betas-beta_lim(1)).*(s==1)).*dJ*alphai;
        
% s.t. constraint
if nargin > 5
%     betas   = [betas(:,1) betas];
    betas   = [beta0 betas];
    for i = 1:Nc
        dbeta           = diff(betas(:,i:i+1),1,2);
        betas(:,i+1)    = (abs(dbeta)>dbeta_max).*(betas(:,i)+...
                            dbeta_max.*sign(dbeta))+...
                            (abs(dbeta)<=dbeta_max).*betas(:,i+1);
    end
    betas = betas(:,2:end);
end

beta    = [betas betas(:,end)*ones(1,length(beta)-Nc)];
beta    = betalim2(beta,beta_lim,dbeta_max);

