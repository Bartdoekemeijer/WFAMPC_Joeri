function beta = betalim2(beta,beta_lim,dbeta_max)

[a,b]       = size(beta);
dbeta       = diff(beta')';

for ai = 1:a
    for bi = 1:b-1
        if beta(ai,bi) > beta_lim(2)
            beta(ai,bi) = beta_lim(2);
        elseif beta(ai,bi) < beta_lim(1)
            beta(ai,bi) = beta_lim(1);
        end
        
        if nargin > 2
            % max dbeta
            if dbeta(ai,bi) > dbeta_max
                beta(ai,bi+1) = beta(ai,bi)+dbeta_max;
            elseif dbeta(ai,bi) < -dbeta_max
                beta(ai,bi+1) = beta(ai,bi)-dbeta_max;
            end
        end
            
    end
    if beta(ai,b) > beta_lim(2)
        beta(ai,b) = beta_lim(2);
    elseif beta(ai,b) < beta_lim(1)
        beta(ai,b) = beta_lim(1);
    end
end

end