function grad = gradproj(grad,beta,beta_lim,eps)

betaup      = (beta' > beta_lim(2) - eps) & (grad < 0);
betadown    = (beta' < beta_lim(1) + eps) & (grad > 0);

grad        = grad.*(1-betaup).*(1-betadown);
